#include "rtk_prepare.h"
#include <cmath>
#include <cstring>
bool ReadNextObservationEpoch(FILE* fp, station_state_t& station)
{
	if (fp == nullptr || station.eof) 
	{
		return false;
	}
	while (true)
	{
		int ret = input_oem4f(&station.raw, fp);
		//文件结束
		if (ret == -2)
		{
			station.eof = true;
			return false;
		}
		//报文错误，例如CRC错误，跳过后继续读取
		if (ret < 0)
		{
			continue;
		}
		//ret==2表明星历已更新，ret==0无用报文
		if (ret != 1)
		{
			continue;
		}
		//观测历元卫星数不足，无效历元
		if (station.raw.obs.n <= 0)
		{
			continue;
		}
		//清空上一历元处理结果
		memset(station.sats, 0, sizeof(station.sats));
		station.spp_sol = sol_t();
		station.valid_sat = 0;
		return true;
	}
}
//比较基准站和流动站的历元时间
int CompareStationEpochTime(const station_state_t& base,
	const station_state_t& rover,
	const rtk_config_t& config,
	double& dt)
{
	dt = 0.0;
	if (base.raw.obs.n <= 0 || rover.raw.obs.n <= 0)
	{
		return 0;
	}
	const gtime_t& base_time = base.raw.obs.data[0].time;
	const gtime_t& rover_time = base.raw.obs.data[0].time;
	dt = timediff(base_time, rover_time);
	//时间差小于0.001
	if (fabs(dt) < config.sync_tolerance)
	{
		return 0;
	}
	if (dt < 0.0)
	{
		return -1;//基准站较早
	}
	return 1;//流动站较早
}
