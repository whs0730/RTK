#pragma once

#ifdef _WIN32
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <winsock2.h>
#include <ws2tcpip.h>
#ifdef _MSC_VER
#pragma comment(lib, "Ws2_32.lib")
#endif
#else
#include <arpa/inet.h>
#include <cerrno>
#include <fcntl.h>
#include <sys/select.h>
#include <sys/socket.h>
#include <unistd.h>
#endif

namespace socket_platform
{
#ifdef _WIN32
using Handle = SOCKET;
inline constexpr Handle invalid = INVALID_SOCKET;
inline int LastError() { return WSAGetLastError(); }
inline bool Interrupted(int error) { return error == WSAEINTR; }
inline bool WouldBlock(int error) { return error == WSAEWOULDBLOCK; }
inline void Close(Handle socket) { closesocket(socket); }
inline bool SetNonblocking(Handle socket)
{
    u_long enabled = 1;
    return ioctlsocket(socket, FIONBIO, &enabled) == 0;
}
#else
using Handle = int;
inline constexpr Handle invalid = -1;
inline int LastError() { return errno; }
inline bool Interrupted(int error) { return error == EINTR; }
inline bool WouldBlock(int error) { return error == EAGAIN || error == EWOULDBLOCK; }
inline void Close(Handle socket) { close(socket); }
inline bool SetNonblocking(Handle socket)
{
    const int flags = fcntl(socket, F_GETFL, 0);
    return flags != -1 && fcntl(socket, F_SETFL, flags | O_NONBLOCK) != -1;
}
#endif
}
