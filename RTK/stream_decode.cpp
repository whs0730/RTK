#include "stream_decode.h"               
#include<string>
#include <cstring>                       // memcpy, memmove 等
#include <windows.h>                     // Windows API

#pragma warning(disable:4996)            // 关闭 VS 对某些旧函数（如 inet_addr）的安全警告

// 打开 TCP 连接
static bool OpenSocket(SOCKET& sock, const char IP[], const unsigned short Port)
{
    WSADATA wsaData;                         // 保存 Winsock 初始化信息
    SOCKADDR_IN addrSrv;                     // 保存 IPv4 地址和端口

    if (!WSAStartup(MAKEWORD(2, 2), &wsaData))
    {
        // 创建 TCP 套接字
        // 本接口只做同步轮询，不发起重叠I/O。
        if ((sock = WSASocketW(AF_INET, SOCK_STREAM, IPPROTO_TCP, nullptr, 0, 0)) != INVALID_SOCKET)
        {
            addrSrv.sin_addr.S_un.S_addr = inet_addr(IP);   // 把字符串 IP 转成网络字节序地址
            addrSrv.sin_family = AF_INET;                  // IPv4
            addrSrv.sin_port = htons(Port);                // 把端口从主机字节序转成网络字节序

            // 连接服务器
            if (connect(sock, (SOCKADDR*)&addrSrv, sizeof(SOCKADDR)) == SOCKET_ERROR)
            {
                closesocket(sock);                         // 连接失败，关闭套接字
                sock = INVALID_SOCKET;                     // 将句柄置为无效
                WSACleanup();                              // 清理 Winsock
                return false;
            }
            return true;                                  
        }
        WSACleanup();                                      // socket() 失败，清理 Winsock
    }
    return false;                                        
}

// 关闭 TCP 连接并清理 Winsock

static void CloseSocket(SOCKET& sock)
{
    closesocket(sock);                                     // 关闭套接字
    WSACleanup();                                          // 清理 Winsock 资源
}

// 析构函数：确保关闭连接

TcpOem4Stream::~TcpOem4Stream()
{
    Close();                                               // 关闭连接，释放资源
}

// 打开指定 IP 和端口的 TCP 流，并设置接收超时
bool TcpOem4Stream::Open(const char* ip, unsigned short port, int recv_timeout_ms)
{
    Close();                                               // 先关闭可能已存在的连接

    if (!OpenSocket(sock_, ip, port)) {                    // 尝试打开 TCP 连接
        SetLastSocketError("OpenSocket failed");           // 记录错误信息
        Close();                                           // 清理资源
        return false;
    }
    wsa_started_ = true;                                   // 标记 Winsock 已启动

    u_long nonblocking = 1;
    if (ioctlsocket(sock_, FIONBIO, &nonblocking) == SOCKET_ERROR)
    {
        SetLastSocketError("ioctlsocket(FIONBIO) failed");
        Close();
        return false;
    }
    receive_timeout_ms_ = recv_timeout_ms;

    last_error_.clear();                                   // 清除上一次错误
    return true;
}

// 关闭 TCP 连接

void TcpOem4Stream::Close()
{
    if (sock_ != INVALID_SOCKET || wsa_started_) {         // 如果套接字有效或曾经初始化过 Winsock
        CloseSocket(sock_);                                // 关闭套接字并清理 Winsock
        sock_ = INVALID_SOCKET;                            // 置为无效套接字
        wsa_started_ = false;                              // 标记 Winsock 未启动
    }
}


// 读取一个字节

bool TcpOem4Stream::ReadByte(uint8_t* data)
{
    return ReadExact(data, 1);                             // 委托给 ReadExact，读取 1 字节
}

// 精确读取指定长度的数据（循环接收直到读满 len 字节）

bool TcpOem4Stream::ReadExact(uint8_t* data, int len)
{
    if (data == nullptr || len <= 0 || sock_ == INVALID_SOCKET) {
        last_error_ = "stream is not open";
        return false;
    }

    int got = 0;
    while (got < len) {
        const int n = RecvPacket(data + got, len - got, 0);
        if (n < 0) return false;
        got += n;
    }
    return true;
}


// 接收一个数据包（非精确读取，一次 recv 调用）
// 参数：data    - 接收缓冲区
//       max_len - 缓冲区最大长度
//       sleep_ms- 接收前等待时间（毫秒），>0 时先 Sleep
int TcpOem4Stream::RecvPacket(uint8_t* data, int max_len, int sleep_ms)
{
    // 参数检查
    if (data == nullptr || max_len <= 0 || sock_ == INVALID_SOCKET) {
        last_error_ = "stream is not open";                
        return -1;
    }

    if (sleep_ms > 0) {
        Sleep(sleep_ms);                                    // 接收前等待
    }

    fd_set readable;
    FD_ZERO(&readable);
    FD_SET(sock_, &readable);
    timeval timeout{ receive_timeout_ms_ / 1000, (receive_timeout_ms_ % 1000) * 1000 };
    // 用select等待数据；不让日常轮询触发阻塞recv的SO_RCVTIMEO。
    const int ready = select(0, &readable, nullptr, nullptr, &timeout);
    if (ready == 0) return 0;
    if (ready == SOCKET_ERROR)
    {
        const int error = WSAGetLastError();
        if (error == WSAEINTR) return 0;
        last_error_ = "select failed: " + to_string(error);
        return -1;
    }

    int n = recv(sock_, reinterpret_cast<char*>(data), max_len, 0);
    if (n > 0) {
        return n;                                           // 成功接收到数据，返回字节数
    }

    if (n == 0) {
        last_error_ = "remote host closed the connection";  // 连接关闭
        return -1;
    }

    int err = WSAGetLastError();                            // 获取错误
    // 可读状态发生变化时非阻塞recv仍可能暂无数据，留到下次轮询。
    if (err == WSAEINTR || err == WSAEWOULDBLOCK) {
        return 0;
    }

    last_error_ = "recv failed: " + to_string(err);    // 记录错误
    return -1;
}
// 返回最后一次错误信息
const std::string& TcpOem4Stream::LastError() const
{
    return last_error_;
}
// 设置最后一次错误信息，附加前缀和当前 Winsock 错误码
void TcpOem4Stream::SetLastSocketError(const char* prefix)
{
    last_error_ = std::string(prefix) + ": " + to_string(WSAGetLastError());
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
