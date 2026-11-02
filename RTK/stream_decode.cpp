#include "stream_decode.h"
#include <chrono>
#include <cstring>
#include <thread>

TcpOem4Stream::~TcpOem4Stream()
{
    Close();
}

bool TcpOem4Stream::Open(const char* ip, unsigned short port, int recv_timeout_ms)
{
    Close();
#ifdef _WIN32
    WSADATA data;
    const int startup_error = WSAStartup(MAKEWORD(2, 2), &data);
    if (startup_error != 0)
    {
        last_error_ = "WSAStartup failed: " + to_string(startup_error);
        return false;
    }
    wsa_started_ = true;
    // 同步轮询不使用重叠 I/O，保留 Winsock 的 flags=0 创建方式。
    sock_ = WSASocketW(AF_INET, SOCK_STREAM, IPPROTO_TCP, nullptr, 0, 0);
#else
    sock_ = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
#endif
    if (sock_ == socket_platform::invalid)
    {
        SetLastSocketError("socket failed");
        Close();
        return false;
    }

    sockaddr_in address{};
    address.sin_family = AF_INET;
    address.sin_port = htons(port);
    if (inet_pton(AF_INET, ip, &address.sin_addr) != 1)
    {
        last_error_ = string("Invalid IPv4 address: ") + ip;
        Close();
        return false;
    }
    if (connect(sock_, reinterpret_cast<sockaddr*>(&address), sizeof(address)) != 0)
    {
        SetLastSocketError("connect failed");
        Close();
        return false;
    }
    if (!socket_platform::SetNonblocking(sock_))
    {
        SetLastSocketError("set nonblocking failed");
        Close();
        return false;
    }
    receive_timeout_ms_ = recv_timeout_ms;
    last_error_.clear();
    return true;
}

void TcpOem4Stream::Close()
{
    if (sock_ != socket_platform::invalid)
    {
        socket_platform::Close(sock_);
        sock_ = socket_platform::invalid;
    }
#ifdef _WIN32
    if (wsa_started_)
    {
        WSACleanup();
        wsa_started_ = false;
    }
#endif
}

bool TcpOem4Stream::ReadByte(uint8_t* data)
{
    return ReadExact(data, 1);
}

bool TcpOem4Stream::ReadExact(uint8_t* data, int len)
{
    if (data == nullptr || len <= 0 || sock_ == socket_platform::invalid)
    {
        last_error_ = "stream is not open";
        return false;
    }
    int got = 0;
    while (got < len)
    {
        const int n = RecvPacket(data + got, len - got, 0);
        if (n < 0) return false;
        got += n;
    }
    return true;
}

int TcpOem4Stream::RecvPacket(uint8_t* data, int max_len, int sleep_ms)
{
    if (data == nullptr || max_len <= 0 || sock_ == socket_platform::invalid)
    {
        last_error_ = "stream is not open";
        return -1;
    }
    if (sleep_ms > 0)
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
    }

    fd_set readable;
    FD_ZERO(&readable);
    FD_SET(sock_, &readable);
    timeval timeout{ receive_timeout_ms_ / 1000, (receive_timeout_ms_ % 1000) * 1000 };
#ifdef _WIN32
    const int nfds = 0;
#else
    // POSIX select 需要最大文件描述符加一；Winsock 忽略该参数。
    const int nfds = sock_ + 1;
#endif
    const int ready = select(nfds, &readable, nullptr, nullptr, &timeout);
    if (ready == 0) return 0;
    if (ready < 0)
    {
        const int error = socket_platform::LastError();
        if (socket_platform::Interrupted(error)) return 0;
        last_error_ = "select failed: " + to_string(error);
        return -1;
    }

    const int n = static_cast<int>(recv(sock_, reinterpret_cast<char*>(data), max_len, 0));
    if (n > 0) return n;
    if (n == 0)
    {
        last_error_ = "remote host closed the connection";
        return -1;
    }
    const int error = socket_platform::LastError();
    // 可读状态变化后 recv 仍可能暂无数据，留到下一轮；真正的断开仍返回错误。
    if (socket_platform::Interrupted(error) || socket_platform::WouldBlock(error)) return 0;
    last_error_ = "recv failed: " + to_string(error);
    return -1;
}

const std::string& TcpOem4Stream::LastError() const
{
    return last_error_;
}

void TcpOem4Stream::SetLastSocketError(const char* prefix)
{
    last_error_ = std::string(prefix) + ": " + to_string(socket_platform::LastError());
}

// 尝试从缓冲区 Buff 中解码 OEM4 报文
// 参数：raw       - 解码后的 Oem4 数据结构
//       Buff      - 存放接收数据的缓冲区
//       lenD      - 缓冲区中有效数据长度
//       ReadFlag  - 输出解码结果标志：
//                    0 正常，负值表示错误，-1 超长报文等
// 返回：true 表示取出了一个报文（成功或错误都算取走），false 需要更多数据

static bool TryDecodeOem4FromStreamBuff(raw_t* raw, uint8_t* Buff, int& lenD, int& ReadFlag)
{
    if (raw == nullptr) {                                   
        ReadFlag = -2;                                      // 标志错误码 -2
        return true;                                        
    }

    while (lenD >= 3) {                                     // 至少需要 3 字节才能搜索同步头
        int sync_pos = -1;

        // 在 Buff 中搜索 3 字节同步头
        for (int i = 0; i <= lenD - 3; i++) {
            if (Buff[i] == OEM4SYNC1 &&
                Buff[i + 1] == OEM4SYNC2 &&
                Buff[i + 2] == OEM4SYNC3) {
                sync_pos = i;                                // 同步头位置
                break;
            }
        }

        // 未找到同步头，保留最后两字节（可能包含部分同步头），丢弃其余
        if (sync_pos < 0) {
            Buff[0] = Buff[lenD - 2];                       // 将倒数第二字节移到开头
            Buff[1] = Buff[lenD - 1];                       // 将倒数第一字节移到开头
            lenD = 2;                                       // 更新有效长度
            return false;                                   // 需要更多数据
        }

        // 同步头前面有无效数据，丢弃它们
        if (sync_pos > 0) {
            lenD -= sync_pos;                               // 新长度
            memmove(Buff, Buff + sync_pos, lenD);           // 将同步头及之后数据移到缓冲区头部
        }

        // 至少需要 10 字节头部
        if (lenD < 10) {
            return false;                                   // 等待更多数据
        }

        // 计算报文长度：len = 消息长度（从第9字节开始的2字节） + 头部固定长度
        raw->len = U2(Buff + 8) + OEM4HLEN;                 
        if (raw->len > MAXRAWLEN - 4) {                     
            lenD--;                                         // 丢弃当前第一个字节，尝试重新同步
            memmove(Buff, Buff + 1, lenD);
            raw->nbyte = 0;
            ReadFlag = -1;                                  // 返回错误 -1（报文过长）
            return true;                                    // 虽然错误，但已取走该报文框架
        }

        int frame_len = raw->len + 4;                       // 总帧长 = 报文长度 + CRC(4字节)
        if (lenD < frame_len) {                             // 缓冲区数据不够一帧
            return false;                                   // 等待更多数据
        }

        // 拷贝完整帧到 raw->buff
        memcpy(raw->buff, Buff, frame_len);
        raw->nbyte = 0;

        // 移除已取走的帧数据
        lenD -= frame_len;
        if (lenD > 0) {
            memmove(Buff, Buff + frame_len, lenD);          // 将剩余数据前移
        }

        ReadFlag = decode_oem4(raw);                        // 调用实际解码函数，返回解码结果
        return true;                                        // 成功取出并解码一帧
    }

    return false;                                           // 缓冲区不足 3 字节，需更多数据
}


// 从 TCP 流中不断接收数据，并尝试解码出 OEM4 报文
// 参数：raw     - 输出解码后的数据
//       stream  - TCP 流对象
// 返回：>0 解码成功，0 无报文或等待数据，-2 错误

int input_oem4s(raw_t* raw, TcpOem4Stream* stream)
{
    static uint8_t buff[OEM4RECVLEN];                        // 临时接收缓冲区
    static uint8_t Buff[MAXRAWLEN + OEM4RECVLEN];            // 未完成帧加一次recv
    static int lenD = 0;                                    // 拼接缓冲区有效数据长度

    if (raw == nullptr || stream == nullptr) {               
        return -2;
    }

    int lenR = 0;
    int ReadFlag = 0;

    // 首先尝试从已有的拼接缓冲区中解码
    if (TryDecodeOem4FromStreamBuff(raw, Buff, lenD, ReadFlag)) {
        return ReadFlag;                                     // 取出一帧，返回解码标志
    }

    // 缓冲区数据不够一帧，从网络接收新数据
    lenR = stream->RecvPacket(buff, OEM4RECVLEN);             // 接收一次，通常是一次轮询
    if (lenR < 0) {
        return -2;                                           // 接收错误
    }

    if (lenR <= 0) {
        return 0;                                            // 无数据可读
    }

    // 防止拼接缓冲区溢出
    if ((lenD + lenR) > MAXRAWLEN + OEM4RECVLEN) {
        lenD = 0;                                            // 溢出则清空，重新同步
    }

    // 将新接收的数据追加到拼接缓冲区
    memcpy(Buff + lenD, buff, lenR);
    lenD += lenR;

    // 再次尝试解码
    if (TryDecodeOem4FromStreamBuff(raw, Buff, lenD, ReadFlag)) {
        return ReadFlag;
    }

    return 0;                                                // 仍然不够一帧，等待下次调用
}
