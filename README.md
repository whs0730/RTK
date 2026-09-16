# SPP 单点定位与测速程序说明

本工程用于解码 NovAtel OEM4 二进制数据，并完成 GPS/BDS 双系统双频单点定位与测速。程序支持两种数据来源：

- 离线文件：读取 `.log` 二进制观测文件。
- 实时流：通过 TCP 读取 NovAtel OEM4 实时数据流。

主流程为：

```text
读取 OEM4 数据
-> 解码观测值和广播星历
-> 按历元组织观测值
-> 计算卫星发射时刻位置、速度、钟差、钟速
-> 地球自转、对流层、TGD 等改正
-> 双频 IF 组合
-> GPS/BDS SPP 最小二乘定位
-> 多普勒测速
-> 输出结果文件
```

## 目录结构

```text
SPP/
  main.cpp                  主程序流程：输入、解码、卫星状态计算、SPP、测速、输出
  define.h                  常量定义：系统编号、频率、消息 ID、PRN 范围等

  obs.h / obs.cpp           观测值、星历、卫星状态、解算结果结构体和卫星编号转换
  decode.h / decode.cpp     NovAtel OEM4 二进制数据解码
  stream_decode.h / .cpp    实时 TCP 流读取
  timeTransform.h / .cpp    GPST、BDT、MJD、公历时间转换
  coordinate.h / .cpp       XYZ、BLH、ENU 坐标转换和高度角计算
  satpos.h / satpos.cpp     GPS/BDS 广播星历卫星位置、速度、钟差、钟速计算
  error_correction.h / .cpp 对流层、电离层、地球自转等误差改正
  spp.h / spp.cpp           IF 组合、质量控制、SPP 定位和测速
  matrix.h / matrix.cpp     简单矩阵运算和最小二乘求逆
```

## 编译方法

### 方法一：Visual Studio

1. 打开根目录下的 `SPP.slnx`。
2. 选择 `x64` 和 `Debug` 或 `Release`。
3. 生成解决方案。
4. 可执行文件通常位于项目目录下：

```text
SPP\x64\Debug\SPP.exe
```

### 方法二：命令行 g++

在工程根目录执行：

```powershell
E:\MSYS2\ucrt64\bin\g++.exe -std=c++17 -finput-charset=UTF-8 -fexec-charset=UTF-8 `
  .\SPP\main.cpp .\SPP\satpos.cpp .\SPP\decode.cpp .\SPP\spp.cpp `
  .\SPP\timeTransform.cpp .\SPP\matrix.cpp .\SPP\obs.cpp `
  .\SPP\coordinate.cpp .\SPP\error_correction.cpp .\SPP\stream_decode.cpp `
  -o .\build\SPP_current.exe -lws2_32
```

## 运行方法

### 离线文件模式

推荐进入源码项目目录后运行：

```powershell
cd "D:\Course Study\Satellite navigation algorithm\Algorithm\SPP\SPP"
.\x64\Debug\SPP.exe
```

如果不传文件名，程序默认读取：

```text
NovatelOEM20211114-01.log
```

注意：默认文件名是相对当前运行目录的。如果在上一级目录运行，需要显式写出文件路径：

```powershell
cd "D:\Course Study\Satellite navigation algorithm\Algorithm\SPP"
.\SPP\x64\Debug\SPP.exe .\SPP\NovatelOEM20211114-01.log
```

### 实时流模式

使用默认 IP 和端口：

```powershell
cd "D:\Course Study\Satellite navigation algorithm\Algorithm\SPP"
.\x64\Debug\SPP.exe --stream
```

指定 IP 和端口：

```powershell
cd "D:\Course Study\Satellite navigation algorithm\Algorithm\SPP"
.\x64\Debug\SPP.exe --stream 8.148.22.229 7003
```

实时模式下，终端会输出类似：

```text
Wk 2420 SOW 393399.000 XYZ ... ENU ... n=23 PDOP=1.109
```

其中：

- `Wk`：GPS 周
- `SOW`：周内秒
- `XYZ`：接收机 ECEF 坐标
- `ENU`：相对参考坐标的东、北、天误差
- `n`：当前可用卫星数
- `PDOP`：空间几何精度因子

## 输出文件

程序输出文件的位置与运行的那个 `SPP.exe` 对应的 `SPP` 目录有关。可以简单理解为：你运行哪一级 `SPP` 目录下的 `SPP.exe`，输出文件就放在哪一级 `SPP` 目录里。

这一点很重要：

- 如果在 Visual Studio 中直接调试离线程序，通常运行的是项目目录下的 `SPP\SPP\x64\Debug\SPP.exe`，输出文件会在 `D:\Course Study\Satellite navigation algorithm\Algorithm\SPP\SPP`。
- 如果在上一级目录终端运行实时流，例如在 `D:\Course Study\Satellite navigation algorithm\Algorithm\SPP>` 下执行 `.\x64\Debug\SPP.exe --stream`，输出文件会在上一级目录 `D:\Course Study\Satellite navigation algorithm\Algorithm\SPP`。
- 因此，同一个程序在 VS 调试和上一级终端运行时，输出位置可能不同。这是相对路径文件名的正常行为。
- 建议离线和实时流分别使用不同级别目录下的 `SPP.exe`：离线用 `SPP\SPP` 这一级，实时流用上一级 `SPP` 这一级。这样两种模式都会生成的 `observations1.txt`、`satpos1.txt`、`spp_solution.txt` 不会互相覆盖，便于对照分析。

当前代码在离线和实时模式都会输出：

```text
observations1.txt
satpos1.txt
spp_solution.txt
```

实时模式还会额外输出：

```text
stream.txt
realtime_raw_oem.log
```

### `observations1.txt`

输出每个历元解码出的观测值，主要用于检查解码是否正常。

字段大致为：

```text
week sow sat P L D SNR
```

其中：

- `P`：伪距，单位 m
- `L`：载波相位，单位周
- `D`：多普勒，单位 Hz
- `SNR`：信噪比，单位 dB-Hz

### `satpos1.txt`

低频抽样输出接收时刻的广播星历卫星位置，方便和参考卫星位置文件对照。

注意：这个文件输出的是“接收时刻广播星历结果”；SPP 解算内部使用的是经过发射时刻和地球自转改正后的卫星状态。

### `spp_solution.txt`

输出单点定位和测速结果，主要用于精度分析。

主要字段包括：

```text
Wk SOW
ECEF-X ECEF-Y ECEF-Z
REF-ECEF-X REF-ECEF-Y REF-ECEF-Z
EAST NORTH UP
B L H
VX VY VZ
PDOP SigmaP SigmaV
GS BS n
```

其中：

- `EAST/NORTH/UP`：相对参考坐标的 ENU 误差
- `B/L/H`：大地坐标
- `VX/VY/VZ`：接收机速度
- `GS`：参与解算的 GPS 卫星数
- `BS`：参与解算的 BDS 卫星数
- `n`：当前历元成功计算出卫星状态的 GPS/BDS 卫星数

### `realtime_raw_oem.log`

实时模式下保存接收到的有效 OEM4 观测/星历报文，后续可以作为离线文件回放。

### `stream.txt`

实时模式下保存控制台实时解算摘要，内容与控制台输出一致，格式类似：

```text
Wk 2421 SOW 42124.000 XYZ ... ENU ... n=23 PDOP=1.144
```

它比 `spp_solution.txt` 更简洁，适合实时观察。

## 离线对比参考文件

根目录下额外放了 3 个参考结果文件，用于离线运行后对照检查程序输出是否合理。

### `NovatelOEM20211114-01.obs`

参考观测值文件，用于和程序输出的 `observations1.txt` 对比。

主要检查内容：

- 历元时间是否一致
- 卫星号是否一致
- 伪距、载波、多普勒、信噪比是否在合理范围
- 是否出现明显缺频或错误卫星号

### `NovatelOEM20211114-01.pos`

参考定位结果文件，用于和程序输出的 `spp_solution.txt` 对比。

主要检查内容：

- `Wk`、`SOW` 历元是否对齐
- `ECEF-X/Y/Z` 是否接近
- `EAST/NORTH/UP` 误差变化趋势是否一致
- `PDOP` 和用星数是否在合理范围

程序中的离线参考坐标 `OFFLINE_REF_X/OFFLINE_REF_Y/OFFLINE_REF_Z` 与该 `.pos` 文件中的 `REF-ECEF` 保持一致，因此离线结果的 ENU 误差可以直接用于精度分析。实时流使用另一组测站参考坐标 `STREAM_REF_X/STREAM_REF_Y/STREAM_REF_Z`。

### `广播星历计算结果 - 241015.txt`

参考广播星历卫星位置文件，用于和程序输出的 `satpos1.txt` 对比。

主要检查内容：

- 历元是否按相同间隔抽样
- 卫星号是否一致
- 卫星位置、速度、钟差、钟速是否接近

注意：`satpos1.txt` 输出的是接收时刻的广播星历卫星状态，主要用于和该参考文件做格式和数值对照；SPP 解算内部使用的是发射时刻并经过地球自转改正后的卫星状态。

## 重要算法说明

### 1. GPS/BDS 时间系统

观测历元按 GPST 组织。GPS 星历直接使用 GPST 计算。BDS 星历原始 toe/toc 属于 BDT，解码时会保存原始 BDT 秒，同时转换为 GPST，方便与观测历元求时间差。

### 2. 卫星发射时刻

SPP 不直接使用接收时刻卫星位置，而是：

```text
PIF / c 得到近似传播时间
-> 得到近似发射时刻
-> 计算卫星钟差
-> 修正发射时刻
-> 重新计算卫星位置、速度、钟差、钟速
-> 做地球自转改正
```

### 3. TGD 改正

卫星钟差计算中不改 TGD。TGD 在 `GetPIF()` 中按伪距组合处理：

- GPS L1/L2 IF 组合不额外改 TGD
- BDS B1I/B3I IF 组合中处理 BDS TGD

### 4. 质量控制

SPP 使用统一筛选函数 `PassSppBasicCheck()`，主要检查：

- 只使用 GPS/BDS
- 双频伪距有效
- 星历有效
- 卫星位置有效
- SNR 不低于 30 dB-Hz
- 高度角不低于 10 度

## 常见问题

### PowerShell 提示找不到程序

错误示例：

```powershell
\x64\Debug\SPP.exe --stream
```

PowerShell 不会默认从当前目录查找程序，应写成：

```powershell
.\x64\Debug\SPP.exe --stream
```

### 实时模式出现 `time/stat skip`

类似输出：

```text
time/stat skip type=218 week=0 tow=0.000 stat=20
```

这表示收到的某些 OEM4 报文时间状态无效，程序跳过它们是正常的。只要后面持续输出 `Wk/SOW/XYZ/ENU/n/PDOP`，说明实时解算正在正常运行。

### 实时模式 ENU 误差较大

当前 ENU 是相对代码中的固定参考坐标计算(-2267810.173，5009324.109，3221016.632)的。如果实时流接收机位置不是该参考点，ENU 出现十几米或几十米是正常现象，不一定代表定位算法异常。
