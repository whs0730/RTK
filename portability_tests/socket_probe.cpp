#include "stream_decode.h"
#include <array>
#include <chrono>
#include <cstdlib>
#include <iostream>

static void Require(bool condition, const std::string& message)
{
    if (!condition)
    {
        std::cerr << message << '\n';
        std::exit(EXIT_FAILURE);
    }
}

int main(int argc, char** argv)
{
    Require(argc == 2, "expected loopback port");
    TcpOem4Stream stream;
    const auto port = static_cast<unsigned short>(std::stoi(argv[1]));
    Require(stream.Open("127.0.0.1", port, 50), stream.LastError());
    std::array<uint8_t, 5> bytes{};
    const auto start = std::chrono::steady_clock::now();
    Require(stream.RecvPacket(bytes.data(), bytes.size(), 0) == 0, "quiet socket must time out");
    Require(std::chrono::steady_clock::now() - start >= std::chrono::milliseconds(30),
        "select returned before the receive timeout");
    Require(stream.ReadExact(bytes.data(), bytes.size()), stream.LastError());
    Require(bytes == std::array<uint8_t, 5>{0, 1, 127, 128, 255}, "fragmented payload mismatch");
    Require(stream.RecvPacket(bytes.data(), bytes.size(), 0) == -1, "EOF must be an error");
    Require(stream.LastError() == "remote host closed the connection", "EOF error was lost");
    stream.Close();
    stream.Close();

    Require(stream.Open("127.0.0.1", port, 50), stream.LastError());
    Require(stream.ReadByte(bytes.data()) && bytes[0] == 42, "reopen/read byte failed");
    stream.Close();
    Require(!stream.Open("invalid-ip", port, 50), "invalid IPv4 address accepted");
    Require(stream.LastError().find("Invalid IPv4 address") != std::string::npos,
        "invalid address error was lost during cleanup");
    Require(!stream.Open("127.0.0.1", 0, 50), "connect to closed port succeeded");
    Require(stream.LastError().find("connect failed") != std::string::npos,
        "connect error was lost during cleanup");
    std::cout << "socket transport passed\n";
}
