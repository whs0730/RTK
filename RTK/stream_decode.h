#pragma once
//Windows 头文件
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif

#include "decode.h"

#include <winsock2.h>//Windows socket 网络通信
#include <cstdint>
#include <string>

#pragma comment(lib, "Ws2_32.lib")//链接 Winsock 库
//TCP实时流读取器类
class TcpOem4Stream {
public:
    TcpOem4Stream() = default;
    ~TcpOem4Stream();
    //用IP 和端口连接服务器，默认接收超时 5000 ms
    bool Open(const char* ip, unsigned short port, int recv_timeout_ms = 5000); 
    void Close();//关闭 socket
    //读取数据
    bool ReadByte(uint8_t* data);//读 1 个字节 
    bool ReadExact(uint8_t* data, int len);//必须读满指定长度 len 才返回 true
    //一次尽量读取一批数据，最多 max_len 字节
    int RecvPacket(uint8_t* data, int max_len, int sleep_ms = 980);
    const string& LastError() const;//返回最近一次错误信息

private:
    void SetLastSocketError(const char* prefix);

    SOCKET sock_ = INVALID_SOCKET;//保存 socket 句柄。初值是无效 socket
    bool wsa_started_ = false;//记录 Winsock 是否已经启动
    std::string last_error_;//保存错误信息
};
//实时解码入口
int input_oem4s(raw_t* raw, TcpOem4Stream* stream);
