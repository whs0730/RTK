#include "rtk_stream.h"

#include <cmath>
#include <cstring>
#include <deque>
#include <vector>

// 接收端状态持续更新，队列中的station_state_t则是观测报文到达时的值副本。
// SPP先验另存于接收端，只在队列历元真正完成SPP后更新，避免接收和处理进度混淆。
struct rtk_stream_receiver_t
{
    TcpOem4Stream connection;
    raw_t raw;
    vector<uint8_t> buffer;
    size_t remaining_skip_bytes = 0; // 不处理的报文也按声明总长度跳过，允许跨多次recv完成
    deque<unique_ptr<station_state_t>> epochs;
    gtime_t latest_observation{};
    bool has_observation = false;
    sol_t last_spp;
    bool has_last_spp = false;
};

static void DiscardBytes(vector<uint8_t>& buffer, size_t count)
{
    buffer.erase(buffer.begin(), buffer.begin() + count);
}

// 完整帧及长度已由调用方检查；这里只读OEM头的周/毫秒，先于有状态的观测解码去重。
// 最新时刻只由成功解码的新历元推进，重复或倒序报文不应改变锁定时间与半周历史。
static bool IsObservationNotNew(rtk_stream_receiver_t& receiver,
    rtk_stream_receiver_statistics_t& statistics)
{
    if (!receiver.has_observation) return false;
    uint32_t milliseconds = 0;
    memcpy(&milliseconds, receiver.buffer.data() + 16, sizeof(milliseconds));
    const int week = U2(receiver.buffer.data() + 14);
    if (week == 0) return false;
    const gtime_t time = gpst2time(adjgpsweek(week), milliseconds * 0.001);
    const double dt = timediff(time, receiver.latest_observation);
    if (dt > 1e-9) return false;
    if (fabs(dt) <= 1e-9) ++statistics.duplicate_epochs;
    else ++statistics.discarded_epochs;
    return true;
}

// TCP读取边界不等于OEM报文边界：缓存未完成帧，循环处理一次recv中的多帧。
// 两站不用input_oem4s的静态组帧缓存；每站独立组帧后再调用共用报文解码器。
static void DecodeReceivedFrames(rtk_stream_receiver_t& receiver,
    rtk_stream_receiver_statistics_t& statistics)
{
    auto& buffer = receiver.buffer;
    while (!buffer.empty())
    {
        if (receiver.remaining_skip_bytes > 0)
        {
            const size_t count = buffer.size() < receiver.remaining_skip_bytes ? buffer.size() : receiver.remaining_skip_bytes;
            DiscardBytes(buffer, count);
            receiver.remaining_skip_bytes -= count;
            continue;
        }
        if (buffer.size() < 3) return;
        size_t start = 0;
        while (start + 2 < buffer.size() &&
            !(buffer[start] == OEM4SYNC1 && buffer[start+1] == OEM4SYNC2 && buffer[start+2] == OEM4SYNC3))
        {
            ++start;
        }
        if (start + 2 == buffer.size())
        {
            // 本批未找到完整同步字，保留末两字节以识别跨recv边界的三字节同步字。
            DiscardBytes(buffer, buffer.size() - 2);
            return;
        }
        if (start > 0) DiscardBytes(buffer, start);
        if (buffer.size() < 10) return;
        if (buffer[3] != OEM4HLEN)
        {
            // 本模块只接收28字节长头格式；头长不符时移一字节重新搜索同步字。
            ++statistics.decode_errors;
            DiscardBytes(buffer, 1);
            continue;
        }
        const size_t frame_size = static_cast<size_t>(U2(buffer.data() + 8)) + OEM4HLEN + 4;
        const int type = U2(buffer.data() + 4);
        if (type != ID_RANGECMP && type != ID_GPSEPHEM && type != ID_BDSEPHEMERIS)
        {
            // RANGE必须在解码前跳过，避免同秒重复交付并修改锁定/LLI历史。
            receiver.remaining_skip_bytes = frame_size;
            continue;
        }
        if (frame_size > MAXRAWLEN)
        {
            ++statistics.decode_errors;
            receiver.remaining_skip_bytes = frame_size;
            continue;
        }
        if (buffer.size() < frame_size) return;
        const size_t payload_size = frame_size - OEM4HLEN - 4;
        // 现有星历解码器按固定偏移读取字段，短报文不能交给它更新导航缓存。
        const size_t minimum_payload = type == ID_GPSEPHEM ? 204 : type == ID_BDSEPHEMERIS ? 196 : 4;
        bool payload_valid = payload_size >= minimum_payload;
        if (payload_valid && type == ID_RANGECMP)
        {
            // RANGECMP正文为4字节记录数加每条24字节记录，声明数量不能越过正文末尾。
            uint32_t observations = 0;
            memcpy(&observations, buffer.data() + OEM4HLEN, sizeof(observations));
            payload_valid = observations <= (payload_size - 4) / 24;
        }
        if (!payload_valid)
        {
            ++statistics.decode_errors;
            DiscardBytes(buffer, frame_size);
            continue;
        }
        if (type == ID_RANGECMP && IsObservationNotNew(receiver, statistics))
        {
            DiscardBytes(buffer, frame_size);
            continue;
        }
        memcpy(receiver.raw.buff, buffer.data(), frame_size);
        // decode_oem4的len不含末尾4字节CRC，但buff保留CRC供该解码器校验。
        receiver.raw.len = static_cast<int>(frame_size) - 4;
        receiver.raw.nbyte = 0;
        DiscardBytes(buffer, frame_size);
        const int result = decode_oem4(&receiver.raw);
        if (result < 0)
        {
            ++statistics.decode_errors;
            continue;
        }
        if (result != 1 || receiver.raw.obs.n <= 0) continue;

        // 与离线RTK一致，每个完整RANGECMP报文交付一个观测历元。
        // 导航数据一同快照，避免稍后接收的星历改变队列中旧历元的计算条件。
        auto epoch = make_unique<station_state_t>();
        epoch->raw = receiver.raw;
        receiver.latest_observation = receiver.raw.obs.data[0].time;
        receiver.has_observation = true;
        ++statistics.observation_epochs;
        receiver.epochs.push_back(std::move(epoch));
        if (receiver.epochs.size() > 5)
        {
            receiver.epochs.pop_front();
            ++statistics.discarded_epochs;
        }
    }
}

RtkStream::RtkStream() = default;

RtkStream::~RtkStream()
{
    Close();
}

bool RtkStream::Open(const rtk_stream_options_t& options)
{
    Close();
    error_message_.clear();
    statistics_ = rtk_stream_statistics_t{};
    const auto& config = options.preparation;
    if (options.base_ip.empty() || options.rover_ip.empty() || options.base_port == 0 || options.rover_port == 0 ||
        options.receive_timeout_ms <= 0 || !isfinite(config.sync_tolerance) || config.sync_tolerance <= 0.0 ||
        !isfinite(config.elevation_mask) || config.elevation_mask < 0.0 || config.elevation_mask >= 90.0 ||
        !isfinite(config.snr_mask) || config.snr_mask < 0.0)
    {
        error_message_ = "Invalid RTK stream addresses, timeouts or preparation options";
        return false;
    }
    options_ = options;
    base_first_ = true;
    base_ = make_unique<rtk_stream_receiver_t>();
    rover_ = make_unique<rtk_stream_receiver_t>();
    base_->buffer.reserve(2 * MAXRAWLEN);
    rover_->buffer.reserve(2 * MAXRAWLEN);
    if (!base_->connection.Open(options.base_ip.c_str(), options.base_port, options.receive_timeout_ms))
    {
        error_message_ = "Base stream: " + base_->connection.LastError();
        Close();
        return false;
    }
    if (!rover_->connection.Open(options.rover_ip.c_str(), options.rover_port, options.receive_timeout_ms))
    {
        error_message_ = "Rover stream: " + rover_->connection.LastError();
        Close();
        return false;
    }
    open_ = true;
    return true;
}

void RtkStream::Close()
{
    open_ = false;
    base_.reset();
    rover_.reset();
}

bool RtkStream::IsOpen() const
{
    return open_;
}

bool RtkStream::Receive(rtk_stream_receiver_t& receiver,
    rtk_stream_receiver_statistics_t& statistics, const char* name)
{
    uint8_t data[MAXRAWLEN];
    // 双站轮询不沿用单站接收接口默认的980 ms休眠。
    const int count = receiver.connection.RecvPacket(data, MAXRAWLEN, 0);
    if (count < 0)
    {
        error_message_ = string(name) + " stream: " + receiver.connection.LastError();
        return false;
    }
    if (count > 0)
    {
        statistics.received_bytes += static_cast<uint64_t>(count);
        receiver.buffer.insert(receiver.buffer.end(), data, data + count);
        DecodeReceivedFrames(receiver, statistics);
    }
    return true;
}

bool RtkStream::PrepareCachedEpoch(rtk_epoch_t& epoch)
{
    while (!base_->epochs.empty() && !rover_->epochs.empty())
    {
        const auto& b = base_->epochs.front()->raw.obs.data[0].time;
        const auto& r = rover_->epochs.front()->raw.obs.data[0].time;
        const double dt = timediff(b, r);
        if (fabs(dt) >= options_.preparation.sync_tolerance)
        {
            // 两队列均按观测时间递增，较早队头不可能匹配另一站更晚的历元，可直接丢弃。
            // 严格同步门限与文件模式相同；不把到达时间当作卫星观测时间。
            if (dt < 0.0)
            {
                base_->epochs.pop_front();
                ++statistics_.base.discarded_epochs;
            }
            else
            {
                rover_->epochs.pop_front();
                ++statistics_.rover.discarded_epochs;
            }
            continue;
        }
        auto base = std::move(base_->epochs.front());
        auto rover = std::move(rover_->epochs.front());
        base_->epochs.pop_front();
        rover_->epochs.pop_front();
        ++statistics_.synchronized_epochs;
        // 每个配对只处理一次。两站SPP独立更新先验，一站失败不抹掉另一站的成功先验。
        base->last_spp_sol = base_->last_spp;
        base->has_last_spp = base_->has_last_spp;
        rover->last_spp_sol = rover_->last_spp;
        rover->has_last_spp = rover_->has_last_spp;
        const bool base_ok = ProcessStationEpoch(*base);
        const bool rover_ok = ProcessStationEpoch(*rover);
        if (base_ok)
        {
            base_->last_spp = base->last_spp_sol;
            base_->has_last_spp = true;
        }
        if (rover_ok)
        {
            rover_->last_spp = rover->last_spp_sol;
            rover_->has_last_spp = true;
        }
        if (!base_ok || !rover_ok)
        {
            ++statistics_.spp_failed_epochs;
            continue;
        }
        if (!BuildRtkEpoch(*base, *rover, options_.preparation, epoch))
        {
            ++statistics_.preparation_failed_epochs;
            continue;
        }
        ++statistics_.prepared_epochs;
        return true;
    }
    epoch = rtk_epoch_t{};
    return false;
}

RtkStreamReadResult RtkStream::ReadPreparedEpoch(rtk_epoch_t& epoch)
{
    epoch = rtk_epoch_t{};
    if (!open_)
    {
        if (error_message_.empty()) error_message_ = "RTK stream is not open";
        return RtkStreamReadResult::Error;
    }
    if (PrepareCachedEpoch(epoch)) return RtkStreamReadResult::EpochReady;
    // 交替先轮询的测站，避免总让同一路先等待超时；一次调用最多交付一个准备历元。
    bool ok = false;
    if (base_first_)
    {
        ok = Receive(*base_, statistics_.base, "Base") && Receive(*rover_, statistics_.rover, "Rover");
    }
    else
    {
        ok = Receive(*rover_, statistics_.rover, "Rover") && Receive(*base_, statistics_.base, "Base");
    }
    base_first_ = !base_first_;
    if (!ok)
    {
        // 连接丢失后不带旧队列继续交付；重连后的解码历史和EKF都需要重新开始。
        Close();
        return RtkStreamReadResult::Error;
    }
    return PrepareCachedEpoch(epoch) ? RtkStreamReadResult::EpochReady : RtkStreamReadResult::Waiting;
}

const string& RtkStream::LastError() const
{
    return error_message_;
}

const rtk_stream_statistics_t& RtkStream::Statistics() const
{
    return statistics_;
}
