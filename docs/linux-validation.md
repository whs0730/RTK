# Linux 适配验证记录

日期：2026-11-02。工作区：`D:/Course Study/Satellite navigation algorithm2/RTK-linux`；分支 `linux-support`；起始 HEAD `39430206f13997e3d54c91b6dc177cfd5bbe17e6`。未提交、推送或合并分支。

## 修改范围

- `CMakeLists.txt`：显式列出 RTK 入口与现有计算模块，Windows 链接 Winsock，测试依赖只在 BUILD_TESTING 开启时启用。
- `RTK/socket_platform.h`、`stream_decode.h/.cpp`：隔离句柄、错误码、关闭与非阻塞设置；Linux 使用 POSIX 的 `select(fd + 1, ...)`，Windows 保留非重叠 Winsock 创建方式。
- `RTK/rtk_main.cpp`：Linux SIGINT / SIGTERM、标准 fopen、通用路径分隔符、退出时接收字节统计。模式、地址参数、坐标、Ratio 与估计器默认值不变。
- `RTK/obs.cpp`：七处卫星编号字符串格式化改用带缓冲区大小的 snprintf。
- `scripts/run_rooftop.sh`：独立运行目录、计时 SIGTERM、控制台日志、实际退出码。
- `portability_tests/`：Socket 及 Linux 退出/文件关闭测试；`.gitattributes` 保证 shell 脚本为 LF。

SPP、卫星位置、坐标、误差改正、准备层、双差函数模型、随机模型、LS/EKF、LAMBDA 的源文件未修改。`stream_decode.cpp` 中从 `TryDecodeOem4FromStreamBuff` 前的说明开始到文件末尾与起始版本逐字一致。

## 已通过的构建与测试

| 环境 | 工具 | 结果 |
| --- | --- | --- |
| WSL2 内 Ubuntu 24.04.3 x86_64 隔离用户空间 | GCC 13.3.0 / CMake 3.28.3 | Release 构建；CTest 2/2 通过 |
| WSL2 Ubuntu 22.04.5 x86_64 | GCC 11.4.0 / CMake 3.22.1 | Release 构建；CTest 2/2 通过 |
| Windows x64 | GCC 15.2.0 / CMake | Release 构建；Socket CTest 1/1 通过 |
| Windows x64 | MSVC 19.51.36260 / CMake / VS 2026 | Release 构建；Socket CTest 1/1 通过 |

Ubuntu 24.04 验证使用 [Ubuntu 官方 base rootfs](https://cdimage.ubuntu.com/ubuntu-base/releases/24.04/release/) 及已校验的软件包，运行于现有 WSL2 内的 chroot，不是 Azure 虚拟机实测。为避开 NTFS 大小写同名头文件冲突，头文件与 CMake 包解至 WSL 原生临时目录后参与构建；临时 rootfs 的 apt 安装没有完整配置完毕，不作为可部署系统。构建和 CTest 的日志为本地 `build/ubuntu2404/validation.log`。Azure 应按运行指南从源码正常安装工具并构建。

Socket 测试覆盖：静默超时、分片二进制接收、ReadExact、EOF、重复关闭、重开、非法 IPv4、连接失败及清理后错误保留。

Linux 退出测试分别向静默的双站实时进程发送 SIGINT / SIGTERM，验证退出码 0、最终统计、两份 CSV 表头及七份文本输出；另验证对端 EOF 返回非零并正常关闭输出。启动脚本用一秒本地双站连接验证了计时终止、日志、退出码 0，以及切换输出目录后的相对二进制路径解析。

已有本地回归测试四份全部通过：

- `rtk_model_tests.cpp`：默认单元测试分支。
- `rtk_estimation_tests.cpp`：默认单元测试分支。
- `rtk_ambiguity_tests.cpp`：默认单元测试分支。
- `rtk_shared_ephemeris_tests.cpp`：原零基线 OEM 输入的所需前缀，覆盖共同星历与 GEO 排除。

这些测试只读取原工程文件；副本、二进制和日志全部放在 Linux 工作区的 `build/regression/`。共享星历测试副本仅把测试入口中的 fopen_s 替换成 fopen。未修改原测试文件，未运行全量离线解算。本地回归测试和输入原本不随 Git 提交；新兼容性测试可由本分支 CTest 独立复现。

## WSL2 真实楼顶流短测

运行环境为上述 Ubuntu 22.04，连接 `8.148.22.229:2002/3002`，运行计时 90 秒后发送 SIGTERM。

| 项目 | 结果 |
| --- | ---: |
| 基站 / 流动站接收字节 | 833204 / 812712 |
| 基站 / 流动站观测历元 | 90 / 90 |
| 同步历元 | 90 |
| 两站 SPP 成功 / 准备 / 最终定位 | 88 / 88 / 88 |
| FIX / FLOAT | 86 / 2 |
| 基站 / 流动站丢弃历元 | 0 / 0 |
| 基站 / 流动站解码错误 | 0 / 0 |
| GPST 定位范围 | 2443 / 103482～103569 |

两份 CSV 均为 88 个数据行；定位时刻严格递增，坐标全部为有限值。日志有 `Stop signal: 15`、`Finished.` 及完整最终统计。原始短测使用未加 --preserve-status 的 timeout 包装器，未可靠采集该包装器退出码；不能补记为 0。RTK 的受控退出码 0 由独立信号测试和最终启动脚本测试验证。

证据位于 `build/live-smoke/run.log`、两份 CSV 与 `summary.json`，本地结果未纳入 Git。定位输出总量约 14.0 MB。该短测不构成四小时稳定性或绝对定位精度验收。

## 原工程未改动验证

适配前在 `build/original_snapshot.json` 记录原 `RTK` 工作区的 40 个已跟踪文件 SHA-256、分支、HEAD 和 Git 状态；适配后核对全部一致：

- 分支仍为 `main`。
- HEAD 仍为 `39430206f13997e3d54c91b6dc177cfd5bbe17e6`。
- 工作区状态与开始时相同。
- 40 个文件无哈希变化。

适配没有启动 Azure 四小时任务，也未创建或修改 Azure 资源。下一步按 [Linux 运行指南](linux.md) 在目标 Azure 主机上完成四小时测试并归档整个输出目录。
