# RTK 定位程序

基于 C++ 的 GPS/BDS 短基线 RTK 程序，支持 NovAtel OEM 二进制文件和双站 TCP 数据流。提供单历元最小二乘、逐历元 EKF 和 LAMBDA 整周模糊度固定，默认使用动态 EKF。

## 主要功能

- 两站独立解码与 SPP、历元同步、共同星历选择和共视卫星筛选。
- GPS、BDS 或双系统的单频/双频双差非组合模型。
- 等方差及高度角定权，保留共享参考星引入的双差相关性。
- EKF 参考星变换、卫星增减、LLI/GF/MW 周跳检测与模糊度重置。
- LAMBDA 候选搜索、Ratio 检验和固定后位置更新；未固定时保留有效浮点解。

## 编译

环境为 Windows、MSVC `v145`、C++20 和 Windows SDK，工程已配置 UTF-8 编译。

1. 用 Visual Studio 打开 [RTK.slnx](RTK.slnx)，或直接打开 [RTK.vcxproj](RTK/RTK.vcxproj)。
2. 选择 **Debug / x64**，生成解决方案。
3. 可执行文件位于 `x64\Debug\RTK.exe`。

Debug / x64 已启用 `rtk_main.cpp` 并排除旧的 `spp_main.cpp`。使用其他配置时，应确认只编译一个 `main()`。

## 运行

以下命令从仓库根目录开始。输入数据和输出文件使用相对路径，均以进程的工作目录为准。

### 离线文件

```powershell
Set-Location .\RTK
..\x64\Debug\RTK.exe
```

当前默认读取 `RTK/OEM7data/20240301/` 下的 `base.log`、`rove.log`。数据未随仓库提交，需要自行准备。

切换数据时，修改 `rtk_main.cpp` 中的 `DEFAULT_INPUT_BASE_FILE`、`DEFAULT_INPUT_ROVER_FILE`，并检查基站坐标与输出后缀，然后重新编译。当前入口不接受命令行文件路径参数。

### 双站实时流

```powershell
.\x64\Debug\RTK.exe --stream
.\x64\Debug\RTK.exe --stream-rooftop
```

`--stream` 默认使用 RTS1 两路数据流及其已知基站坐标；`--stream-rooftop` 使用楼顶两路数据流，并将首次有效基站 SPP 坐标固定为本次运行的概略基准。

两种模式均可追加四个参数覆盖连接地址：

```text
RTK.exe --stream base_ip base_port rover_ip rover_port
RTK.exe --stream-rooftop base_ip base_port rover_ip rover_port
```

覆盖地址不会自动更换基站坐标，接入其他数据源时需同时检查几何基准配置。运行中按一次 `Ctrl+C`，等待程序关闭文件并输出统计；当前没有自动重连或断点续跑。

查看命令帮助：`.\x64\Debug\RTK.exe --help`。

## 配置与代码入口

| 文件                                                                                                   | 职责与主要配置                                                                            |
| ------------------------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------- |
| [rtk_main.cpp](RTK/rtk_main.cpp)                                                                       | 输入路径、输出后缀、基站坐标及主流程；`ProcessPreparedEpoch()` 串联单历元构模、估计和固定 |
| [rtk_prepare.h](RTK/rtk_prepare.h) / [cpp](RTK/rtk_prepare.cpp)                                        | 两站同步与数据准备；默认同步差 `< 0.001 s`、高度角 `≥ 10°`、双频 SNR `> 30 dB-Hz`         |
| [rtk_stream.h](RTK/rtk_stream.h) / [cpp](RTK/rtk_stream.cpp)                                           | 双站 TCP 输入、独立解码缓存与配对；每站最多缓存 5 个历元，满时暂停接收                    |
| [rtk_mathematical_model.h](RTK/rtk_mathematical_model.h) / [cpp](RTK/rtk_mathematical_model.cpp)       | 系统、频点、参考星及双差模型；默认 GPS+BDS 双频，禁止 BDS GEO 作参考星                    |
| [rtk_stochastic_model.h](RTK/rtk_stochastic_model.h) / [cpp](RTK/rtk_stochastic_model.cpp)             | 等方差/高度角定权；主流程采用高度角定权，伪距/相位方差比为 `10000`                        |
| [rtk_estimation.h](RTK/rtk_estimation.h) / [cpp](RTK/rtk_estimation.cpp)                               | LS/EKF 估计、状态传递和周跳处理；默认 `ExtendedKalman`、`kinematic=true`                  |
| [rtk_ambiguity_resolution.h](RTK/rtk_ambiguity_resolution.h) / [cpp](RTK/rtk_ambiguity_resolution.cpp) | LAMBDA 搜索及固定；默认 `Ratio > 2.0` 才固定                                              |

估计方式、运动模式和系统/频点目前通过源码配置，没有命令行切换开关。使用单历元 LS 时，在主入口创建 `estimation_options` 后将其 `method` 设为 `RtkEstimationMethod::LeastSquares`；静态 EKF 将 `kinematic` 设为 `false`。

准备层目前按双频筛选共视集合，后续选择单频模型不会恢复因第二频点缺失而被剔除的卫星。两站 SPP 独立保留，RTK 卫星状态采用共同星历重算；已知基站坐标不会覆盖原 SPP 结果。

动态 EKF 每历元用 SPP 重置位置先验，跨历元传递可映射的模糊度及其协方差。整数固定只更新当前输出，不反馈滤波器，也不进行 fix-and-hold。异常解算被质量检查拒绝时，不提交新滤波状态。

## 输出

结果写入**运行时工作目录**。当前离线后缀为 `_3_1`，实时模式分别为 `_realtime`、`_realtime_rooftop`；离线后缀由源码设置，不随输入文件自动变化。

下表中的 `{suffix}` 表示上述后缀：

| 文件                                                              | 内容                                   |
| ----------------------------------------------------------------- | -------------------------------------- |
| `base_observations{suffix}.txt`、`rover_observations{suffix}.txt` | 解码观测值 P/L/D/SNR                   |
| `rtk_prepared{suffix}.txt`                                        | 共视星列表及两站 SPP 摘要              |
| `rtk_function_model{suffix}.txt`                                  | 参考星、线性化坐标、设计矩阵及双差残差 |
| `rtk_stochastic_model_Elevation{suffix}.txt`                      | 双差协方差阵 D 与权阵 P                |
| `rtk_estimation{suffix}.txt`、`rtk_float{suffix}.csv`             | 固定前浮点参数、协方差及估计统计       |
| `rtk_ambiguity{suffix}.txt`、`rtk_solution{suffix}.csv`           | 固定信息及最终选用的 FIX/FLOAT 结果    |

`rtk_float` 中的状态为 `FLOAT`，最终定位状态以 `rtk_solution` 为准。CSV 的 `dX_m/dY_m/dZ_m` 为流动站减基站的基线分量，ENU 以基站为原点。相位观测单位为周，模型残差和基线单位为米。

同一工作目录内重复运行会覆盖同名文件。保留不同算例或运行结果时，应使用独立工作目录或提前归档。离线观测文本在读取时导出，实时观测文本只保存准备成功的配对历元，两者历元数的含义不同。
