# Ubuntu / Azure 实时 RTK 运行

目标环境：Ubuntu 24.04 x86_64、GCC、CMake。适配工作位于独立的 `RTK-linux` 工作区及 `linux-support` 分支。SPP、RTK、EKF、LAMBDA、双差模型及其默认参数保持不变。

## 编译

在 Ubuntu 24.04 上进入本分支的仓库根目录：

~~~bash
sudo apt-get update
sudo apt-get install -y build-essential cmake
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF
cmake --build build -j 4
./build/RTK --help
~~~

运行程序不需要 Python 或第三方数值库。CMake 只选择 `rtk_main.cpp`，旧单站入口 `spp_main.cpp` 不参与构建；两站 SPP 仍由 RTK 主流程调用。

兼容性测试需要 Python 标准库：

~~~bash
sudo apt-get install -y python3
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build -j 4
ctest --test-dir build --output-on-failure
~~~

Windows 可继续用本工作区的 Visual Studio 工程，或使用 CMake；Windows 分支链接 Winsock，Linux 分支使用 POSIX Socket。原始 `RTK` 工作区与 `main` 无需修改。

## 四小时楼顶测试

先完成编译，然后从仓库根目录运行：

~~~bash
bash scripts/run_rooftop.sh 4h
~~~

默认连接基站 `8.148.22.229:2002`、流动站 `8.148.22.229:3002`。覆盖 IP / 端口：

~~~bash
bash scripts/run_rooftop.sh 4h BASE_IP BASE_PORT ROVER_IP ROVER_PORT
~~~

每次自动创建唯一的 `runs/rooftop_<UTC时间>_<随机后缀>/`，包含原有 CSV / 文本、`run.log` 和 `exit_status.txt`。脚本先打印结果目录，可另开 SSH 终端执行 `tail -f <结果目录>/run.log`。

脚本使用 GNU `timeout --preserve-status --signal=TERM --kill-after=30s`：计时到达四小时后发送 SIGTERM；RTK 正常完成关闭时退出码为 0，接收或写入失败为非零。30 秒内不能退出才发送 SIGKILL。只有四小时内未提前结束且日志完整，才能记为完整四小时测试。脚本记录起止 UTC 和实际退出码。

断开 SSH 后继续运行：

~~~bash
nohup bash scripts/run_rooftop.sh 4h >rooftop_launcher.log 2>&1 &
~~~

也可直接运行程序，工作目录决定输出位置：

~~~bash
repo_root="$PWD"
mkdir -p runs/manual_rooftop
cd runs/manual_rooftop
"$repo_root/build/RTK" --stream-rooftop >run.log 2>&1
~~~

直接运行仍支持四个地址参数，`--stream` 模式也保留。`--stream-rooftop` 继续固定首次有效的基站 SPP 坐标；覆盖连接地址不会自动提供新的已知基站坐标。

Azure 虚拟机需要能向这两个地址发起 TCP 出站连接；无需为此接收器开放入站监听端口。运行期间不要覆盖输出目录或替换可执行文件。本次适配没有创建 Azure 资源，也没有启动四小时运行。

## 退出与结果检查

Linux 在连接前注册 SIGINT / SIGTERM。处理器只设置退出标志，主循环负责关闭两站 Socket、刷新并关闭所有输出，最后打印统计。运行中可按 Ctrl+C，或对 RTK 进程执行 `kill -TERM <RTK进程PID>`。SIGKILL、机器断电不经过该流程。

正常退出需检查：

- `exit_status.txt` 为 0，`run.log` 有 `Finished.` 和最终统计。
- `rtk_float_realtime_rooftop.csv` 与 `rtk_solution_realtime_rooftop.csv` 的数据行数符合最终统计。
- 最终定位状态以 `rtk_solution` 为准；`rtk_float` 保存固定前浮点结果。
- 同步数、SPP 失败数、解码错误、丢弃数、FIX 数分别记录，不用连接成功代替有效定位验证。

尚未建立成功的连接没有定位 CSV；连接错误保留在运行日志中。连接成功后真实 EOF / 网络错误仍停止并保存已有输出，不自动重连。两路独立缓存、每站五历元满队列暂停接收的策略不变。

结果文件较多，函数模型与协方差文本可能占用较大磁盘空间。先观察一次短测的目录大小，再为四小时输出预留空间。

## 离线输入

默认数据仍为 `OEM7data/20240303/base.log`、`rove.log`，后缀 `_3_3`，坐标配置不变；仅把路径分隔符换为 Windows / Linux 均支持的 `/`。数据未随仓库提交。

本次验证记录见 [linux-validation.md](linux-validation.md)。
