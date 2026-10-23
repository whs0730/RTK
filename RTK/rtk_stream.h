#pragma once

#include "stream_decode.h"
#include "rtk_prepare.h"
#include <cstdint>
#include <memory>

struct rtk_stream_options_t
{
    // 课件RTS1的两路OEM二进制TCP流，调用方可覆盖地址。
    string base_ip = "47.114.134.129";
    unsigned short base_port = 7190;
    string rover_ip = "8.148.22.229";
    unsigned short rover_port = 7201;
    int receive_timeout_ms = 100; // 单路socket接收超时，ms；不是历元同步门限或断流判据
    rtk_config_t preparation;
};

struct rtk_stream_receiver_statistics_t
{
    uint64_t received_bytes = 0;
    uint64_t observation_epochs = 0; // 完成解码并入队的新历元，不含重复、倒序和坏报文
    uint64_t duplicate_epochs = 0;
    uint64_t discarded_epochs = 0;  // 倒序、缓存溢出或已无法配对的历元
    uint64_t decode_errors = 0;
};

struct rtk_stream_statistics_t
{
    rtk_stream_receiver_statistics_t base;
    rtk_stream_receiver_statistics_t rover;
    uint64_t synchronized_epochs = 0; // 两站配对成功即计数，后续SPP/共视筛选仍可能失败
    uint64_t spp_failed_epochs = 0;
    uint64_t preparation_failed_epochs = 0;
    uint64_t prepared_epochs = 0;
};

enum class RtkStreamReadResult
{
    Error = -1,     // 未打开或接收失败；接收失败会关闭两路连接，不自动重连
    Waiting = 0,    // 本轮没有可交付历元，允许继续调用；不是EOF或解算失败
    EpochReady = 1  // 两站同步、SPP和共视准备均成功，尚未构建双差模型
};

struct rtk_stream_receiver_t;

// 两站各自解码RANGE/RANGECMP并缓存最近5个观测历元，配对后复用SPP与共同星历准备流程。
// 不设置函数模型的已知基站坐标；接入实时主流程时应采用RTS1资料中的基站坐标。
// 单实例由同一线程调用；队列仅用于吸收两路到达时延差，不外推、插值或补齐缺失观测。
class RtkStream
{
public:
    RtkStream();
    ~RtkStream();
    RtkStream(const RtkStream&) = delete;
    RtkStream& operator=(const RtkStream&) = delete;

    // 失败时关闭两路连接。重新Open会清空解码、缓存和SPP先验；调用方也应重置EKF。
    bool Open(const rtk_stream_options_t& options = rtk_stream_options_t{});
    // 释放连接与待配对历元，保留统计和错误供主流程在退出后输出。
    void Close();
    bool IsOpen() const;

    // 单线程交替轮询；无完整同步历元时返回Waiting，连接中断时返回Error。
    // 每次调用清空epoch，只有EpochReady时可使用；接收等待上限由两站超时共同决定。
    // 优先消费已有队列；SPP或共视筛选失败的配对会被消费并计数，不作为Error返回。
    RtkStreamReadResult ReadPreparedEpoch(rtk_epoch_t& epoch);
    const string& LastError() const;
    const rtk_stream_statistics_t& Statistics() const;

private:
    bool Receive(rtk_stream_receiver_t& receiver, rtk_stream_receiver_statistics_t& statistics, const char* name);
    bool PrepareCachedEpoch(rtk_epoch_t& epoch);

    unique_ptr<rtk_stream_receiver_t> base_;
    unique_ptr<rtk_stream_receiver_t> rover_;
    rtk_stream_options_t options_;
    rtk_stream_statistics_t statistics_;
    string error_message_;
    bool open_ = false;
    bool base_first_ = true;
};
