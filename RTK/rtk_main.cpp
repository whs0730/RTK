#include "decode.h"
#include "rtk_stream.h"
#include "satpos.h"
#include "coordinate.h"
#include "error_correction.h"
#include "spp.h"
#include "rtk_prepare.h"
#include "rtk_mathematical_model.h"
#include "rtk_stochastic_model.h"
#include "rtk_estimation.h"
#include "rtk_ambiguity_resolution.h"
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <atomic>
#include <memory>
#include <windows.h>
using namespace std;
// 无参数运行时只读取这两个默认文件；工作目录应能解析OEM7data相对路径。
static const char* DEFAULT_INPUT_BASE_FILE = R"(OEM7data\20240303\base.log)";
static const char* DEFAULT_INPUT_ROVER_FILE = R"(OEM7data\20240303\rove.log)";
// OEM7data/readme.txt：20240301和20240303两组算例的已知基准站ECEF坐标，单位m。
// 换其他基站数据时修改此处；没有已知坐标则将USE_KNOWN_BASE_XYZ设为false。
static const bool USE_KNOWN_BASE_XYZ = true;
static const double KNOWN_BASE_XYZ[3] = { -2267808.6227, 5009324.7699, 3221016.8928 };
// RTS1实时流的基站位于工学部，与离线算例的基站坐标不同。
// 覆盖IP/端口不会自动更换几何基准；改用其他基站时须一并核对本坐标。
static const double STREAM_BASE_XYZ[3] = { -2267335.856, 5008649.133, 3222374.912 };
static atomic<bool> stop_requested{false};

static BOOL WINAPI HandleConsoleControl(DWORD event)
{
	if (event == CTRL_C_EVENT || event == CTRL_BREAK_EVENT)
	{
		// 控制事件由系统线程回调，只置标志；输出关闭和连接释放交给主循环完成。
		stop_requested.store(true);
		return TRUE;
	}
	return FALSE;
}

static void PrintUsage(ostream& out)
{
	out << "Usage:\n"
		<< "  RTK.exe\n"
		<< "  RTK.exe --stream [base_ip base_port rover_ip rover_port]\n"
		<< "  RTK.exe --stream-rooftop [base_ip base_port rover_ip rover_port]\n"
		<< "Offline reads the default files configured in rtk_main.cpp.\n"
		<< "--stream uses the RTS1 base coordinates in STREAM_BASE_XYZ.\n"
		<< "Rooftop uses 8.148.22.229:2002/3002 and fixes the first valid base SPP coordinate.\n"
		<< "Press Ctrl+C to finish realtime processing and close output files.\n";
}

static bool ParsePort(const char* text, unsigned short& port)
{
	if (text[0] < '0' || text[0] > '9') return false;
	char* end = nullptr;
	const unsigned long value = strtoul(text, &end, 10);
	if (*end != '\0' || value == 0 || value > 65535) return false;
	port = static_cast<unsigned short>(value);
	return true;
}

// 仅切换默认文件/实时模式，实时可覆盖两路地址；不接受文件名位置参数。
static bool ParseCommandLine(
	int argc,
	char* argv[],
	bool& realtime,
	bool& rooftop,
	rtk_stream_options_t& stream_options)
{
	if (argc == 1) return true;
	if (strcmp(argv[1], "--stream") == 0)
	{
		realtime = true;
		if (argc == 2) return true;
		if (argc != 6) return false;
		stream_options.base_ip = argv[2];
		stream_options.rover_ip = argv[4];
		return !stream_options.base_ip.empty() && !stream_options.rover_ip.empty() &&
			ParsePort(argv[3], stream_options.base_port) && ParsePort(argv[5], stream_options.rover_port);
	}
	else if (strcmp(argv[1], "--stream-rooftop") == 0)
	{
		if (argc != 2 && argc != 6)
		{
			return false;
		}
		realtime = true;
		rooftop = true;
		stream_options.base_ip = "8.148.22.229";
		stream_options.base_port = 2002;
		stream_options.rover_ip = "8.148.22.229";
		stream_options.rover_port = 3002;
		if (argc == 2)
		{
			return true;
		}
		stream_options.base_ip = argv[2];
		stream_options.rover_ip = argv[4];
		return !stream_options.base_ip.empty() && !stream_options.rover_ip.empty() &&
			ParsePort(argv[3], stream_options.base_port) && ParsePort(argv[5], stream_options.rover_port);
	}
	return false;
}

static void InitializeRooftopBase(const rtk_epoch_t& epoch, rtk_model_options_t& options)
{
	if (options.use_fixed_spp_base_xyz)
	{
		return;
	}
	// 基站逐历元SPP会抖动；固定本次运行的概略基准，避免EKF因基站坐标变化反复重置。
	for (int axis = 0; axis < 3; ++axis)
	{
		options.fixed_spp_base_xyz[axis] = epoch.base_spp.XYZ[axis];
	}
	options.use_fixed_spp_base_xyz = true;
	const auto precision = cout.precision();
	cout << setprecision(15) << "Rooftop base reference: first prepared base SPP, fixed approximate ECEF(m)="
		<< options.fixed_spp_base_xyz[0] << ' ' << options.fixed_spp_base_xyz[1] << ' '
		<< options.fixed_spp_base_xyz[2] << endl;
	cout.precision(precision);
}

int main(int argc,char*argv[])
{
	if (argc == 2 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0))
	{
		PrintUsage(cout);
		return EXIT_SUCCESS;
	}
	bool realtime = false;
	bool rooftop = false;
	FILE* base_fp = nullptr;
	FILE* rover_fp = nullptr;
	const string input_base_file = DEFAULT_INPUT_BASE_FILE;
	const string input_rover_file = DEFAULT_INPUT_ROVER_FILE;
	rtk_stream_options_t stream_options;
	RtkStream stream;
	if (!ParseCommandLine(argc, argv, realtime, rooftop, stream_options))
	{
		PrintUsage(cerr);
		return EXIT_FAILURE;
	}
	if (realtime)
	{
		cout << "Opening base stream " << stream_options.base_ip << ':' << stream_options.base_port
			<< " and rover stream " << stream_options.rover_ip << ':' << stream_options.rover_port << endl;
		if (!stream.Open(stream_options))
		{
			cerr << stream.LastError() << endl;
			return EXIT_FAILURE;
		}
	}
	else if (fopen_s(&base_fp, input_base_file.c_str(), "rb") != 0)
	{
		cerr << "Can't open base file: " << input_base_file << endl;
		return -1;
	}
	if (!realtime && fopen_s(&rover_fp, input_rover_file.c_str(), "rb") != 0)
	{
		cerr << "Can't open rover file: " << input_rover_file << endl;
		fclose(base_fp);
		return -1;
	}
	// 实时输出使用独立后缀；同一模式再次运行会覆盖此前同名文件。
	const string output_suffix = rooftop ? "_realtime_rooftop" : realtime ? "_realtime" : "_3_3";
	ofstream base_obs_out("base_observations" + output_suffix + ".txt");
	ofstream rover_obs_out("rover_observations" + output_suffix + ".txt");
	ofstream prepared_out("rtk_prepared" + output_suffix + ".txt");
	//任务二：输出每个有效历元的参考星、行列映射、B矩阵和L向量
	ofstream model_out("rtk_function_model" + output_suffix + ".txt");
	//任务三：按任务二的观测行序保存方差阵D、权阵P及定权参数
	ofstream stochastic_out("rtk_stochastic_model_Elevation" + output_suffix + ".txt");
	// 浮点CSV始终保留估计器输出，最终CSV另存本历元选择的FIX/FLOAT，便于对照审查。
	ofstream estimation_out("rtk_estimation" + output_suffix + ".txt");
	ofstream solution_out("rtk_float" + output_suffix + ".csv");
	//任务五：输出整数候选、Ratio与FIX/FLOAT结果。
	ofstream ambiguity_out("rtk_ambiguity" + output_suffix + ".txt");
	ofstream resolved_out("rtk_solution" + output_suffix + ".csv");
	if (!base_obs_out.is_open() || !rover_obs_out.is_open() ||
		!prepared_out.is_open() || !model_out.is_open() || !stochastic_out.is_open() ||
		!estimation_out.is_open() || !solution_out.is_open() || !WriteRtkSolutionCsvHeader(solution_out) ||
		!ambiguity_out.is_open() || !resolved_out.is_open() || !WriteRtkAmbiguitySolutionCsvHeader(resolved_out))
	{
		cerr << "Cannot open one or more output files" << endl;
		if (base_fp != nullptr) fclose(base_fp);
		if (rover_fp != nullptr) fclose(rover_fp);
		return EXIT_FAILURE;
	}
	// 文件模式由主流程持有两站状态；实时模式由RtkStream持有，避免重复分配导航缓存。
	unique_ptr<station_state_t> base_station;
	unique_ptr<station_state_t> rover_station;
	if (!realtime)
	{
		base_station = make_unique<station_state_t>();
		rover_station = make_unique<station_state_t>();
	}
	rtk_config_t config = realtime ? stream_options.preparation : rtk_config_t{};
	//任务二的系统/频点选项；改成1可构建双系统单频模型
	rtk_model_options_t model_options;
	model_options.frequency_count = 2;
	model_options.use_known_base_xyz = realtime || USE_KNOWN_BASE_XYZ;
	for (int axis = 0; axis < 3; ++axis)
	{
		model_options.known_base_xyz[axis] = realtime ? STREAM_BASE_XYZ[axis] : KNOWN_BASE_XYZ[axis];
	}
	if (rooftop)
	{
		model_options.use_known_base_xyz = false;
	}
	rtk_function_model_t function_model;
	//主流程使用高度角定权，伪距/相位方差比10000
	rtk_stochastic_options_t stochastic_options;
	stochastic_options.method = RtkWeightingMethod::Elevation;
	rtk_stochastic_model_t stochastic_model;
	//任务四默认动态EKF；LeastSquares为历元间完全独立的浮点最小二乘
	rtk_estimation_options_t estimation_options;
	// filter必须跨成功历元持有；单历元构模/估计失败时，估计器保留上一成功状态。
	rtk_filter_state_t filter;
	rtk_solution_t solution;
	rtk_ambiguity_options_t ambiguity_options;
	rtk_ambiguity_solution_t resolved_solution;

	bool output_failed = false;
	bool input_failed = false;
	//读取首个历元并且写入每个站的观测数据
	bool has_base = !realtime && ReadAndExportNext(base_fp, *base_station, base_obs_out, output_failed);
	bool has_rover = !realtime && ReadAndExportNext(rover_fp, *rover_station, rover_obs_out, output_failed);
	// 计数按处理阶段区分：配对成功不等于SPP、准备、估计或固定成功。
	int synchronized_epochs = 0;
	//跳过的历元
	int skipped_base_epochs = 0;
	int skipped_rover_epochs = 0;
	//共同历元且都有SPP解
	int both_spp_ok_epochs = 0;
	//共同历元且都有SPP解且筛选后依旧有效
	int prepared_epochs = 0;
	//任务二：成功构模并写入文件的历元，以及未能构模的历元
	int model_epochs = 0;
	int model_failed_epochs = 0;
	int stochastic_epochs = 0;
	int stochastic_failed_epochs = 0;
	int estimation_epochs = 0;
	int estimation_failed_epochs = 0;
	int fixed_epochs = 0;
	int ambiguity_failed_epochs = 0;
	// 观测与导航结构较大，使用堆对象；本轮只有准备成功后才进入双差解算。
	auto rtk_epoch = make_unique<rtk_epoch_t>();
	if (realtime)
	{
		if (!SetConsoleCtrlHandler(HandleConsoleControl, TRUE))
		{
			cerr << "Cannot register realtime stop handler: " << GetLastError() << endl;
			return EXIT_FAILURE;
		}
		cout << "Realtime processing started. Press Ctrl+C to finish." << endl;
		cout << "Realtime output suffix: " << output_suffix << endl;
	}
	// 文件模式任一站耗尽即结束；实时模式持续轮询，直到Ctrl+C、接收错误或写入失败。
	while ((realtime || (has_base && has_rover)) && !output_failed && !stop_requested.load())
	{
		double dt = 0.0;
		int result = 0;
		bool prepared_ok = false;
		if (realtime)
		{
			const auto read_result = stream.ReadPreparedEpoch(*rtk_epoch);
			const auto& statistics = stream.Statistics();
			synchronized_epochs = static_cast<int>(statistics.synchronized_epochs);
			skipped_base_epochs = static_cast<int>(statistics.base.discarded_epochs);
			skipped_rover_epochs = static_cast<int>(statistics.rover.discarded_epochs);
			both_spp_ok_epochs = static_cast<int>(statistics.synchronized_epochs - statistics.spp_failed_epochs);
			if (read_result == RtkStreamReadResult::Waiting) continue;
			if (read_result == RtkStreamReadResult::Error)
			{
				input_failed = true;
				cerr << "Realtime processing stopped: " << stream.LastError() << endl;
				break;
			}
			// 实时观测文本仅保存成功准备的同步历元，并非全部接收报文的原始流存档。
			// 文件模式则在读入时导出观测，两种模式的观测文件历元数不能直接比较。
			if (!WriteStationObservations(base_obs_out, rtk_epoch->base_obs) ||
				!WriteStationObservations(rover_obs_out, rtk_epoch->rover_obs))
			{
				output_failed = true;
				cerr << "Failed to write realtime station observations" << endl;
				break;
			}
			dt = rtk_epoch->time_diff;
			prepared_ok = true;
		}
		else
		{
			result = CompareStationEpochTime(*base_station, *rover_station, config, dt);
		}
		if (result == 0)
		{
			const gtime_t& base_time = realtime ? rtk_epoch->base_obs.data[0].time : base_station->raw.obs.data[0].time;
			const gtime_t& rover_time = realtime ? rtk_epoch->rover_obs.data[0].time : rover_station->raw.obs.data[0].time;
			cout << "SYNC "
				<< "week=" << base_time.week
				<< " base=" << base_time.sec
				<< " rover=" << rover_time.sec
				<< " dt=" << dt
				<< " base_n=" << (realtime ? rtk_epoch->base_obs.n : base_station->raw.obs.n)
				<< " rover_n=" << (realtime ? rtk_epoch->rover_obs.n : rover_station->raw.obs.n)
				<< endl;
			bool base_ok = true;
			bool rover_ok = true;
			if (!realtime)
			{
				++synchronized_epochs;
				base_ok = ProcessStationEpoch(*base_station);
				rover_ok = ProcessStationEpoch(*rover_station);
				if (base_ok && rover_ok) ++both_spp_ok_epochs;
				prepared_ok = base_ok && rover_ok && BuildRtkEpoch(*base_station, *rover_station, config, *rtk_epoch);
			}
			if (prepared_ok)
			{
				if (rooftop)
				{
					InitializeRooftopBase(*rtk_epoch, model_options);
				}
				cout << "COMMON week=" << rtk_epoch->time.week
					<< " sow=" << rtk_epoch->time.sec
					<< " GPS=" << rtk_epoch->common_gps
					<< " BDS=" << rtk_epoch->common_bds
					<< " total=" << rtk_epoch->ncommon
					<< endl;
				//把共同有效历元数据写入文件
				if (!WritePreparedEpoch(prepared_out, *rtk_epoch))
				{
					output_failed = true;
					cerr << "Failed to write prepared epoch at "
						<< rtk_epoch->time.week << ' '
						<< rtk_epoch->time.sec << endl;
					break;
				}
				prepared_epochs++;
				//流动站用本程序SPP展开点，基站用已知坐标；epoch中的两站SPP结果不覆盖
				if (BuildRtkFunctionModel(*rtk_epoch, rtk_epoch->rover_spp.XYZ,
					config, model_options, function_model))
				{
					if (!WriteRtkFunctionModel(model_out, function_model))
					{
						output_failed = true;
						cerr << "Failed to write function model at epoch "
							<< rtk_epoch->time.week << ' ' << rtk_epoch->time.sec << endl;
						break;
					}
					model_epochs++;
					//只有任务二成功后才能构建任务三，保证历元、参考星和观测顺序一致
					if (BuildRtkStochasticModel(*rtk_epoch, function_model,
						stochastic_options, stochastic_model))
					{
						if (!WriteRtkStochasticModel(stochastic_out, stochastic_model))
						{
							output_failed = true;
							cerr << "Failed to write stochastic model at epoch "
								<< rtk_epoch->time.week << ' ' << rtk_epoch->time.sec << endl;
							break;
						}
						stochastic_epochs++;
						// 估计器在迭代/预报/更新位置处自行重建模型；上面导出的是SPP展开点的模型。
						if (EstimateRtkEpoch(*rtk_epoch, config, model_options, stochastic_options,
							estimation_options, filter, solution))
						{
							if (!WriteRtkSolution(estimation_out, solution) || !WriteRtkSolutionCsv(solution_out, solution))
							{
								output_failed = true;
								cerr << "Failed to write RTK float solution at " << rtk_epoch->time.sec << endl;
								break;
							}
							// 固定结果只用于当前历元输出，filter继续传递任务四的浮点状态。
							if (!ResolveRtkAmbiguities(solution, ambiguity_options, resolved_solution))
							{
								++ambiguity_failed_epochs;
								cerr << "RTK ambiguity resolution retained FLOAT at " << solution.time.sec
									<< ": " << resolved_solution.error_message << endl;
							}
							if (!WriteRtkAmbiguitySolution(ambiguity_out, resolved_solution) ||
								!WriteRtkAmbiguitySolutionCsv(resolved_out, resolved_solution))
							{
								output_failed = true;
								cerr << "Failed to write RTK ambiguity solution at " << solution.time.sec << endl;
								break;
							}
							if (resolved_solution.fixed) ++fixed_epochs;
							++estimation_epochs;
							if (realtime)
							{
								cout << "RTK " << (resolved_solution.fixed ? "FIX" : "FLOAT")
									<< " week=" << resolved_solution.time.week << " sow=" << resolved_solution.time.sec
									<< " ratio=" << resolved_solution.ratio << " baseline_m="
									<< resolved_solution.baseline_xyz[0] << ' ' << resolved_solution.baseline_xyz[1]
									<< ' ' << resolved_solution.baseline_xyz[2] << endl;
								// 长时间运行时逐历元刷新两份结果CSV，详细矩阵文本仍由文件流缓冲。
								solution_out.flush();
								resolved_out.flush();
								if (!solution_out || !resolved_out)
								{
									output_failed = true;
									cerr << "Failed to flush realtime solution output" << endl;
									break;
								}
							}
							else
							{
								cout << "RTK FLOAT week=" << solution.time.week << " sow=" << solution.time.sec
									<< " baseline_m=" << solution.baseline_xyz[0] << ' ' << solution.baseline_xyz[1]
									<< ' ' << solution.baseline_xyz[2] << " carried_N=" << solution.carried_ambiguities << endl;
							}
						}
						else
						{
							++estimation_failed_epochs;
							cerr << "RTK estimation skipped at " << rtk_epoch->time.sec << ": " << solution.error_message << endl;
						}
					}
					else
					{
						stochastic_failed_epochs++;
						cerr << "Stochastic model skipped at epoch "
							<< rtk_epoch->time.week << ' ' << rtk_epoch->time.sec
							<< ": " << stochastic_model.error_message << endl;
					}
				}
				else
				{
					//几何或观测不满足构模条件时跳过，避免输出无效/上一历元模型
					model_failed_epochs++;
					cerr << "Function model skipped at epoch "
						<< rtk_epoch->time.week << ' ' << rtk_epoch->time.sec
						<< ": " << function_model.error_message << endl;
				}
			}
			if (!realtime)
			{
				cout << "SPP "
					<< "week=" << base_time.week
					<< " sow=" << base_time.sec
					<< " base_pos=" << (base_ok ? "OK" : "FAIL")
					<< " base_speed=" << (base_station->speed_ok ? "OK" : "FAIL")
					<< " rover_pos=" << (rover_ok ? "OK" : "FAIL")
					<< " rover_speed=" << (rover_station->speed_ok ? "OK" : "FAIL")
					<< endl;
				//读取下一历元并且写入每个站的观测数据
				has_base = ReadAndExportNext(base_fp, *base_station, base_obs_out, output_failed);
				has_rover = ReadAndExportNext(rover_fp, *rover_station, rover_obs_out, output_failed);
			}

		}
		// 离线输入按时间有序，仅推进较早的一站，保留较晚一站等待下一次配对。
		else if (result < 0)
		{
			skipped_base_epochs++;//跳过一个基准站历元
			has_base = ReadAndExportNext(base_fp, *base_station, base_obs_out, output_failed);
		}
		//流动站历元提前
		else {
			skipped_rover_epochs++;
			has_rover = ReadAndExportNext(rover_fp, *rover_station, rover_obs_out, output_failed);
		}
	}
	if (realtime)
	{
		stream.Close();
		SetConsoleCtrlHandler(HandleConsoleControl, FALSE);
	}
	// close仍可能报告缓冲数据写入失败；正常跳出循环不代表所有输出已成功落盘。
	base_obs_out.close();
	rover_obs_out.close();
	if (!base_obs_out || !rover_obs_out)
	{
		output_failed = true;
		cerr << "Failed to close station observation files" << endl;
	}
	prepared_out.close();
	if (!prepared_out && !output_failed)
	{
		output_failed = true;
		cerr << "Failed to close prepared epoch output after writing" << endl;
	}
	//关闭任务二输出，检查缓冲区最后写入是否成功
	model_out.close();
	if (!model_out)
	{
		output_failed = true;
		cerr << "Failed to close function model output after writing" << endl;
	}
	//检查任务三输出缓冲区最终写入是否成功
	stochastic_out.close();
	if (!stochastic_out)
	{
		output_failed = true;
		cerr << "Failed to close stochastic model output after writing" << endl;
	}
	estimation_out.close();
	solution_out.close();
	if (!estimation_out || !solution_out)
	{
		output_failed = true;
		cerr << "Failed to close RTK estimation output after writing" << endl;
	}
	ambiguity_out.close();
	resolved_out.close();
	if (!ambiguity_out || !resolved_out)
	{
		output_failed = true;
		cerr << "Failed to close RTK ambiguity output after writing" << endl;
	}
	cout << (output_failed ? "Processing stopped: output write failed." :
		input_failed ? "Processing stopped: realtime input failed." : "Finished.") << endl;
	cout << "Synchronized epochs: "
		<< synchronized_epochs << endl;
	cout << "Skipped base epochs: "
		<< skipped_base_epochs << endl;
	cout << "Skipped rover epochs: "
		<< skipped_rover_epochs << endl;
	cout << "Both SPP successful epochs: " << both_spp_ok_epochs << endl;
	cout << "Prepared epochs: " << prepared_epochs << endl;
	cout << "Function model epochs: " << model_epochs << endl;
	cout << "Function model skipped epochs: " << model_failed_epochs << endl;
	cout << "Stochastic model epochs: " << stochastic_epochs << endl;
	cout << "Stochastic model skipped epochs: " << stochastic_failed_epochs << endl;
	cout << "RTK float solution epochs: " << estimation_epochs << endl;
	cout << "RTK estimation skipped epochs: " << estimation_failed_epochs << endl;
	cout << "RTK fixed solution epochs: " << fixed_epochs << endl;
	cout << "RTK ambiguity resolution failed epochs: " << ambiguity_failed_epochs << endl;
	cout << "RTK ratio threshold: " << ambiguity_options.ratio_threshold << endl;
	cout << "RTK estimator: " << (estimation_options.method == RtkEstimationMethod::ExtendedKalman ? "EKF" : "LS")
		<< " positioning=" << (estimation_options.kinematic ? "KINEMATIC" : "STATIC") << endl;
	cout << "Stochastic weighting: "
		<< (stochastic_options.method == RtkWeightingMethod::EqualVariance ? "EQUAL" : "ELEVATION")
		<< " code/phase variance ratio=" << stochastic_options.code_phase_variance_ratio << endl;
	cout << "Function model mode: GPS=" << model_options.use_gps
		<< " BDS=" << model_options.use_bds
		<< " frequencies=" << model_options.frequency_count << endl;
	cout << "Base coordinate source: "
		<< (rooftop ? (model_options.use_fixed_spp_base_xyz ? "SPP_FIXED_APPROX" : "SPP_NOT_INITIALIZED") :
			model_options.use_known_base_xyz ? "KNOWN" : "SPP") << endl;
	if (realtime)
	{
		const auto& statistics = stream.Statistics();
		cout << "Realtime received epochs: base=" << statistics.base.observation_epochs
			<< " rover=" << statistics.rover.observation_epochs << endl;
		cout << "Realtime duplicate epochs: base=" << statistics.base.duplicate_epochs
			<< " rover=" << statistics.rover.duplicate_epochs << endl;
		cout << "Realtime decode errors: base=" << statistics.base.decode_errors
			<< " rover=" << statistics.rover.decode_errors << endl;
		cout << "Realtime SPP/preparation failed epochs: " << statistics.spp_failed_epochs
			<< '/' << statistics.preparation_failed_epochs << endl;
	}
	cout << "Output files: base_observations" << output_suffix << ".txt, rover_observations" << output_suffix
		<< ".txt, rtk_prepared" << output_suffix << ".txt, rtk_function_model" << output_suffix
		<< ".txt, rtk_stochastic_model_Elevation" << output_suffix << ".txt, rtk_estimation" << output_suffix
		<< ".txt, rtk_float" << output_suffix << ".csv, rtk_ambiguity" << output_suffix
		<< ".txt, rtk_solution" << output_suffix << ".csv" << endl;
	if (base_fp != nullptr) fclose(base_fp);
	if (rover_fp != nullptr) fclose(rover_fp);
	return output_failed || input_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
