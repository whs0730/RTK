#pragma once

#include "decode.h"

#include "socket_platform.h"
#include <cstdint>
#include <string>

//TCP实时流读取器类
class TcpOem4Stream {
public:
    TcpOem4Stream() = default;
    ~TcpOem4Stream();
    // 用IP和端口连接服务器，默认等待可读数据5000 ms。
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

    socket_platform::Handle sock_ = socket_platform::invalid;
#ifdef _WIN32
    bool wsa_started_ = false;
#endif
    int receive_timeout_ms_ = 5000;
    std::string last_error_;//保存错误信息
};
//实时解码入口
int input_oem4s(raw_t* raw, TcpOem4Stream* stream);
