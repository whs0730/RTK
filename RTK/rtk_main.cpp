#include "decode.h"
#include "stream_decode.h"
#include "satpos.h"
#include "coordinate.h"
#include "error_correction.h"
#include "spp.h"
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
static void ParseCommandLine(int argc, char* argv[], bool& realtime, string& input_base_file_,string&input_rover_file_, string& stream_ip, unsigned short& stream_port);
static const char* DEFAULT_INPUT_BASE_FILE = ".\OEM7data\20240301\base.log";
static const char* DEFAULT_INPUT_ROVER_FILE = ".\OEM7data\20240301\rover.log";
static const char* DEFAULT_STREAM_IP = "8.148.22.229";
static const unsigned short DEFAULT_STREAM_PORT = 7003;
int main(int argc,char*argv[])
{
	bool realtime = false;
	FILE* BASE_fp = nullptr;
	FILE* ROVER_fp = nullptr;
	string input_base_file = DEFAULT_INPUT_BASE_FILE;
	string input_rover_file = DEFAULT_INPUT_ROVER_FILE;
	string stream_ip = DEFAULT_STREAM_IP;
	unsigned short stream_port = DEFAULT_STREAM_PORT;
	ParseCommandLine(argc, argv, realtime, input_base_file, input_rover_file, stream_ip, stream_port);
}