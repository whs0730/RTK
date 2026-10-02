#include "decode.h"
#include "stream_decode.h"
#include "satpos.h"
#include "coordinate.h"
#include "error_correction.h"
#include "spp.h"
#include "rtk_prepare.h"
#include "rtk_model.h"
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
using namespace std;
static const char* DEFAULT_INPUT_BASE_FILE = R"(OEM7data\20240301\base.log)";
static const char* DEFAULT_INPUT_ROVER_FILE = R"(OEM7data\20240301\rove.log)";
static const char* DEFAULT_STREAM_IP = "8.148.22.229";
static const unsigned short DEFAULT_STREAM_PORT = 7003;
static void ParseCommandLine(
	int argc,
	char* argv[],
	bool& realtime,
	string& input_base_file,
	string& input_rover_file)
{
	if (argc >= 3)
	{
		input_base_file = argv[1];
		input_rover_file = argv[2];
	}
}
int main(int argc,char*argv[])
{
	bool realtime = false;
	//创建文件地址指针
	FILE* base_fp = nullptr;
	FILE* rover_fp = nullptr;
	string input_base_file = DEFAULT_INPUT_BASE_FILE;
	string input_rover_file = DEFAULT_INPUT_ROVER_FILE;
	string stream_ip = DEFAULT_STREAM_IP;
	unsigned short stream_port = DEFAULT_STREAM_PORT;
	//解析命令参数
	ParseCommandLine(argc, argv, realtime, input_base_file, input_rover_file);
	//打开基准站和流动站输入文件
	if (fopen_s(&base_fp, input_base_file.c_str(), "rb") != 0 || base_fp == nullptr)
	{
		cerr << "Can't open base file: " << input_base_file << endl;
		return -1;
	}
	if (fopen_s(&rover_fp, input_rover_file.c_str(), "rb") != 0 || rover_fp == nullptr)
	{
		cerr << "Can't open rover file: " << input_rover_file << endl;
		fclose(base_fp);
		return -1;
	}
	//创建输出两个站观测文件和共视卫星写入流
	ofstream base_obs_out("base_observations_3_1.txt");
	ofstream rover_obs_out("rover_observations_3_1.txt");
	ofstream prepared_out("rtk_prepared_3_1.txt");
	//任务二：输出每个有效历元的参考星、行列映射、B矩阵和L向量
	ofstream model_out("rtk_function_model_3_1.txt");
	//输出文件创建失败返回错误信息
	if (!base_obs_out.is_open() || !rover_obs_out.is_open() ||
		!prepared_out.is_open() || !model_out.is_open())
	{
		cerr << "Cannot open one or more output files" << endl;
		fclose(base_fp);
		fclose(rover_fp);
		return EXIT_FAILURE;
	}
	//创建基准站和流动站状态结构体以及统一配置参数
	station_state_t *base_station=new station_state_t();
	station_state_t *rover_station=new station_state_t();
	rtk_config_t config;
	//任务二的系统/频点选项；改成1可构建双系统单频模型
	rtk_model_options_t model_options;
	model_options.frequency_count = 2;
	rtk_function_model_t function_model;

	bool output_failed = false;
	//读取首个历元并且写入每个站的观测数据
	bool has_base = ReadAndExportNext(base_fp, *base_station, base_obs_out, output_failed);
	bool has_rover = ReadAndExportNext(rover_fp, *rover_station, rover_obs_out, output_failed);
	//共同历元
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
	//为共同历元创建结构体
	rtk_epoch_t* rtk_epoch = new rtk_epoch_t();
	//遍历所有历元
	while (has_base && has_rover && !output_failed)
	{
		double dt = 0.0;
		int result = CompareStationEpochTime(*base_station, *rover_station, config, dt);
		if (result == 0)
		{
			const gtime_t& base_time =base_station->raw.obs.data[0].time;
			const gtime_t& rover_time =rover_station->raw.obs.data[0].time;
			cout << "SYNC "
				<< "week=" << base_time.week
				<< " base=" << base_time.sec
				<< " rover=" << rover_time.sec
				<< " dt=" << dt
				<< " base_n=" << base_station->raw.obs.n
				<< " rover_n=" << rover_station->raw.obs.n
				<< endl;
			synchronized_epochs++;
			bool base_ok = ProcessStationEpoch(*base_station);
			bool rover_ok = ProcessStationEpoch(*rover_station);
			if (base_ok && rover_ok)
			{
				both_spp_ok_epochs++;
			}
			const bool prepared_ok = base_ok && rover_ok &&BuildRtkEpoch(*base_station, *rover_station, config, *rtk_epoch);
			if (prepared_ok)
			{
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
				//流动站SPP坐标作为线性化展开点；基准站坐标由epoch提供
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
		//基准站历元提前
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
	//close观测文件写入流
	base_obs_out.close();
	rover_obs_out.close();
	//close失败返回错误信息
	if (!base_obs_out || !rover_obs_out)
	{
		output_failed = true;
		cerr << "Failed to close station observation files" << endl;
	}
	//close共视文件写入流
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
	//销毁rtk结构体
	delete rtk_epoch;
	//打印结果
	cout << (output_failed ? "Processing stopped: output write failed." : "Finished.") << endl;
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
	cout << "Function model mode: GPS=" << model_options.use_gps
		<< " BDS=" << model_options.use_bds
		<< " frequencies=" << model_options.frequency_count << endl;
	cout << "Output files: base_observations_3_1.txt, "
		<< "rover_observations_3_1.txt, rtk_prepared_3_1.txt, "
		<< "rtk_function_model_3_1.txt" << endl;
	//关闭文件
	fclose(base_fp);
	fclose(rover_fp);
	//返回结果
	return output_failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
