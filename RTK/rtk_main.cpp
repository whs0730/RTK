#include "decode.h"
#include "stream_decode.h"
#include "satpos.h"
#include "coordinate.h"
#include "error_correction.h"
#include "spp.h"
#include "rtk_prepare.h"
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
	FILE* base_fp = nullptr;
	FILE* rover_fp = nullptr;
	string input_base_file = DEFAULT_INPUT_BASE_FILE;
	string input_rover_file = DEFAULT_INPUT_ROVER_FILE;
	string stream_ip = DEFAULT_STREAM_IP;
	unsigned short stream_port = DEFAULT_STREAM_PORT;
	ParseCommandLine(argc, argv, realtime, input_base_file, input_rover_file);
	if (fopen_s(&base_fp, input_base_file.c_str(), "rb") != 0 || base_fp == nullptr)
	{
		cerr << "Can't open base file:" << input_base_file << endl;
		return -1;
	}
	if (fopen_s(&rover_fp, input_rover_file.c_str(), "rb") != 0 || rover_fp == nullptr)
	{
		cerr << "Can't open rover file:" << input_rover_file << endl;
		fclose(base_fp);
		return -1;
	}
	station_state_t *base_station=new station_state_t();
	station_state_t *rover_station=new station_state_t();
	rtk_config_t config;
	bool has_base = ReadNextObservationEpoch(base_fp, *base_station);
	bool has_rover = ReadNextObservationEpoch(rover_fp, *rover_station);
	int synchronized_epochs = 0;
	int skipped_base_epochs = 0;
	int skipped_rover_epochs = 0;
	while (has_base && has_rover)
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

			cout << "SPP "
				<< "week=" << base_time.week
				<< " sow=" << base_time.sec
				<< " base_pos=" << (base_ok ? "OK" : "FAIL")
				<< " base_speed=" << (base_station->speed_ok ? "OK" : "FAIL")
				<< " rover_pos=" << (rover_ok ? "OK" : "FAIL")
				<< " rover_speed=" << (rover_station->speed_ok ? "OK" : "FAIL")
				<< endl;
			has_base =ReadNextObservationEpoch(base_fp,*base_station);
			has_rover = ReadNextObservationEpoch(rover_fp, *rover_station);

		}
		else if (result < 0)
		{
			skipped_base_epochs++;//跳过一个基准站历元
			has_base= ReadNextObservationEpoch(base_fp, *base_station);
		}
		else {
			skipped_rover_epochs++;
			has_rover = ReadNextObservationEpoch(rover_fp, *rover_station);
		}
	}
	cout << "Finished." << endl;
	cout << "Synchronized epochs: "
		<< synchronized_epochs << endl;
	cout << "Skipped base epochs: "
		<< skipped_base_epochs << endl;
	cout << "Skipped rover epochs: "
		<< skipped_rover_epochs << endl;

	fclose(base_fp);
	fclose(rover_fp);

	return 0;
}