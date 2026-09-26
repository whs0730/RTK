#include "rtk_prepare.h"
#include "coordinate.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
using namespace std;
bool ReadNextObservationEpoch(FILE* fp, station_state_t& station)
{
	if (fp == nullptr || station.eof) 
	{
		return false;
	}
	while (true)
	{
		int ret = input_oem4f(&station.raw, fp,true);
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
		station.speed_sol = solvel_t();
		station.valid_sat = 0;
		station.position_ok = false;
		station.speed_ok = false;
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
	const gtime_t& rover_time = rover.raw.obs.data[0].time;
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
bool ProcessStationEpoch(station_state_t& station) {
	// 当前测站、当前历元的观测数据
	obs_t& obs = station.raw.obs;
	//清空上一历元处理结果
	memset(station.sats, 0, sizeof(station.sats));
	station.spp_sol = sol_t();
	station.speed_sol = solvel_t();
	station.valid_sat = 0;
	station.position_ok = false;
	station.speed_ok = false;
	if (obs.n <= 0)
	{
		return false;
	}
	station.valid_sat = CaculateSatellitePositions(obs.data, obs.n, &station.raw.nav, station.sats);
	//spp至少四颗星
	if (station.valid_sat < 4)
	{
		return false;
	}
	//spp解算
	int position_nv = 0;
	const sol_t* initial_solution = nullptr;
	if (station.has_last_spp)
	{
		initial_solution = &station.last_spp_sol;
	}
	station.position_ok = SPP(obs.data, obs.n, &station.raw.nav, &station.spp_sol, station.sats, &position_nv, initial_solution);
	if (!station.position_ok || station.spp_sol.stat != 1)
	{
		station.spp_sol.stat = 0;
		return false;
	}
	//保存位置的成功解，作为下一历元的初值
	station.last_spp_sol = station.spp_sol;
	station.has_last_spp = true;
	//速度解算
	int speed_nv = 0;
	station.speed_ok = SPP_Speed(obs.data, obs.n, &station.spp_sol, station.sats, &station.speed_sol, &speed_nv);
	if (!station.speed_ok || station.speed_sol.stat != 1)
	{
		station.speed_ok = false;
		station.speed_sol.stat = 0;
	}
	return true;//计算出位置就返回成功
}

//计算高度角和方位角
static bool CalculateAzel(const sol_t& sol, const satpos_t& sat, double azel[2])
{
	double dx = sat.pos[0] - sol.XYZ[0];
	double dy = sat.pos[1] - sol.XYZ[1];
	double dz = sat.pos[2] - sol.XYZ[2];
	double rho = sqrt(dx * dx + dy * dy + dz * dz);

	if (!isfinite(rho) || rho < 1.0)
		return false;

	XYZ rec;
	rec.X = sol.XYZ[0];
	rec.Y = sol.XYZ[1];
	rec.Z = sol.XYZ[2];

	BLH* blh = XYZtoBLH(rec, 6378137.0, 1.0 / 298.257223563);
	if (blh == nullptr)
		return false;

	double los[3] = { dx / rho, dy / rho, dz / rho };
	satazel(blh, los, azel);
	delete blh;

	return isfinite(azel[0]) && isfinite(azel[1]);
}

bool PassRtkObservationCheck(
	const obsd_t& base_obs,
	const obsd_t& rover_obs,
	const satpos_t& base_sat,
	const satpos_t& rover_sat,
	const sol_t& base_sol,
	const sol_t& rover_sol,
	const rtk_config_t& config,
	double base_azel[2],
	double rover_azel[2]
)
{
	//不可用
	if (base_sol.stat != 1 || rover_sol.stat != 1)
	{
		return false;
	}

	if (base_obs.sat <= 0 || base_obs.sat != rover_obs.sat)
	{
		return false;
	}

	// sats[i] 必须对应 obs.data[i]
	if (base_sat.sat != base_obs.sat || rover_sat.sat != rover_obs.sat)
	{
		return false;

	}
	//历元时间对准
	if (fabs(timediff(base_obs.time, rover_obs.time)) >= config.sync_tolerance)
	{
		return false;
	}
	//获得卫星系统类型和prn
	int prn = 0;
	int sys = satsys(base_obs.sat, &prn);
	if (sys != SYS_GPS && sys != SYS_CMP)
		return false;

	
	for (int f = 0; f < 2; ++f)
	{
		//频点观测码有效且应相同
		if (base_obs.code[f] == CODE_NONE || base_obs.code[f] != rover_obs.code[f])
		{
			return false;
		}
		//伪距和载波应有效
		if (!std::isfinite(base_obs.P[f]) || base_obs.P[f] <= 0.0 ||
			!std::isfinite(rover_obs.P[f]) || rover_obs.P[f] <= 0.0 ||
			!std::isfinite(base_obs.L[f]) || base_obs.L[f] == 0.0 ||
			!std::isfinite(rover_obs.L[f]) || rover_obs.L[f] == 0.0)
		{
			return false;
		}

		// 信噪比大于 30 dB-Hz
		if (base_obs.SNR[f] <= config.snr_mask ||
			rover_obs.SNR[f] <= config.snr_mask)
			return false;
	}
	//BDS使用指定频点
	if (sys == SYS_CMP && (base_obs.code[0] != CODE_L2I || base_obs.code[1] != CODE_L6I))
	{
		return false;
	}
	//两站观测卫星的高度角和方位角应有效
	if (!CalculateAzel(base_sol, base_sat, base_azel) || !CalculateAzel(rover_sol, rover_sat, rover_azel))
	{
		return false;
	}
	//最后满足高度角限制
	const double mask_rad = config.elevation_mask * PI / 180.0;
	return base_azel[1] >= mask_rad &&rover_azel[1] >= mask_rad;
}

bool BuildRtkEpoch(const station_state_t& base,const station_state_t& rover,const rtk_config_t& config,rtk_epoch_t& epoch)
{
	epoch = rtk_epoch_t();

	if (!base.position_ok || !rover.position_ok || base.raw.obs.n <= 0 || rover.raw.obs.n <= 0)
	{
		return false;
	}

	const gtime_t base_time = base.raw.obs.data[0].time;
	const gtime_t rover_time = rover.raw.obs.data[0].time;
	const double dt = timediff(base_time, rover_time);
	//时间对准
	if (fabs(dt) >= config.sync_tolerance)
	{
		return false;

	}
	
	epoch.time = base_time;
	epoch.time_diff = dt;
	epoch.base_obs = base.raw.obs;
	epoch.rover_obs = rover.raw.obs;
	epoch.base_spp = base.spp_sol;
	epoch.rover_spp = rover.spp_sol;

	// 用卫星编号查流动站下标
	int rover_index[MAXSAT + 1];
	for (int i = 0; i < MAXSAT + 1; i++) {
		rover_index[i] = -1;
	}

	for (int j = 0; j < rover.raw.obs.n; ++j)
	{
		int sat = rover.raw.obs.data[j].sat;
		if (sat >= 1 && sat <= MAXSAT)
		{
			rover_index[sat] = j;
		}
	}

	for (int i = 0; i < base.raw.obs.n; ++i)
	{
		int sat = base.raw.obs.data[i].sat;
		if (sat < 1 || sat > MAXSAT)
		{
			continue;
		}

		int j = rover_index[sat];
		if (j < 0)
		{
			continue;

		}

		double base_azel[2] = {};
		double rover_azel[2] = {};
		//筛选是否可用
		if (!PassRtkObservationCheck(
			base.raw.obs.data[i], rover.raw.obs.data[j],
			base.sats[i], rover.sats[j],
			base.spp_sol, rover.spp_sol, config,
			base_azel, rover_azel))
		{
			continue;

		}
		//结构体赋值
		common_sat_t& item = epoch.common[epoch.ncommon++];
		item.sat = sat;
		item.base_index = i;
		item.rover_index = j;
		item.base_sat = base.sats[i];
		item.rover_sat = rover.sats[j];
		item.base_azel[0] = base_azel[0];
		item.base_azel[1] = base_azel[1];
		item.rover_azel[0] = rover_azel[0];
		item.rover_azel[1] = rover_azel[1];

		int prn = 0;
		int sys = satsys(sat, &prn);
		if (sys == SYS_GPS) ++epoch.common_gps;
		if (sys == SYS_CMP) ++epoch.common_bds;
	}

	return epoch.ncommon > 0;
}
// 保存一个已完成数据准备的同步历元
bool WritePreparedEpoch(ostream& out,const rtk_epoch_t& epoch)
{
	if (!out ||
		epoch.ncommon <= 0 ||
		epoch.ncommon > MAXSAT ||
		epoch.base_obs.n < 0 ||
		epoch.base_obs.n > MAXOBS ||
		epoch.rover_obs.n < 0 ||
		epoch.rover_obs.n > MAXOBS)
	{
		return false;
	}

	// 写入前检查下标
	for (int k = 0; k < epoch.ncommon; ++k)
	{
		const common_sat_t& item = epoch.common[k];

		if (item.base_index < 0 ||
			item.base_index >= epoch.base_obs.n ||
			item.rover_index < 0 ||
			item.rover_index >= epoch.rover_obs.n)
		{
			return false;
		}

		if (epoch.base_obs.data[item.base_index].sat != item.sat ||
			epoch.rover_obs.data[item.rover_index].sat != item.sat)
		{
			return false;
		}
	}
	out << setprecision(15);

	out << "EPOCH "
		<< epoch.time.week << ' '
		<< epoch.time.sec
		<< " dt_s=" << epoch.time_diff
		<< " base_n=" << epoch.base_obs.n
		<< " rover_n=" << epoch.rover_obs.n
		<< " common_gps=" << epoch.common_gps
		<< " common_bds=" << epoch.common_bds
		<< " common_total=" << epoch.ncommon
		<< '\n';

	out << "BASE_SPP "
		<< "XYZ_m="
		<< epoch.base_spp.XYZ[0] << ','
		<< epoch.base_spp.XYZ[1] << ','
		<< epoch.base_spp.XYZ[2]
		<< " clock_GPS_m=" << epoch.base_spp.dtr[0]
		<< " clock_BDS_m=" << epoch.base_spp.dtr[1]
		<< " ns=" << epoch.base_spp.ns
		<< " PDOP=" << epoch.base_spp.pdop
		<< '\n';

	out << "ROVER_SPP "
		<< "XYZ_m="
		<< epoch.rover_spp.XYZ[0] << ','
		<< epoch.rover_spp.XYZ[1] << ','
		<< epoch.rover_spp.XYZ[2]
		<< " clock_GPS_m=" << epoch.rover_spp.dtr[0]
		<< " clock_BDS_m=" << epoch.rover_spp.dtr[1]
		<< " ns=" << epoch.rover_spp.ns
		<< " PDOP=" << epoch.rover_spp.pdop
		<< '\n';

	for (int k = 0; k < epoch.ncommon; ++k)
	{
		const common_sat_t& item = epoch.common[k];

		const obsd_t& base_obs =epoch.base_obs.data[item.base_index];

		const obsd_t& rover_obs =epoch.rover_obs.data[item.rover_index];

		out << "SAT " << sat2id(item.sat)
			<< " base_az_deg=" << item.base_azel[0] * 180.0 / PI
			<< " base_el_deg=" << item.base_azel[1] * 180.0 / PI
			<< " rover_az_deg=" << item.rover_azel[0] * 180.0 / PI
			<< " rover_el_deg=" << item.rover_azel[1] * 180.0 / PI
			<< '\n';

		// [0]/[1]：GPS L1/L2，BDS B1I/B3I
		for (int f = 0; f < 2; ++f)
		{
			out << "BASE_OBS f=" << f
				<< " P_m=" << base_obs.P[f]
				<< " L_cycles=" << base_obs.L[f]
				<< " D_Hz=" << base_obs.D[f]
				<< " SNR_dBHz=" << base_obs.SNR[f]
				<< " code=" << static_cast<int>(base_obs.code[f])
				<< " LLI=" << static_cast<int>(base_obs.LLI[f])
				<< '\n';

			out << "ROVER_OBS f=" << f
				<< " P_m=" << rover_obs.P[f]
				<< " L_cycles=" << rover_obs.L[f]
				<< " D_Hz=" << rover_obs.D[f]
				<< " SNR_dBHz=" << rover_obs.SNR[f]
				<< " code=" << static_cast<int>(rover_obs.code[f])
				<< " LLI=" << static_cast<int>(rover_obs.LLI[f])
				<< '\n';
		}

		out << "BASE_SAT "
			<< "XYZ_m="
			<< item.base_sat.pos[0] << ','
			<< item.base_sat.pos[1] << ','
			<< item.base_sat.pos[2]
			<< " clock_s=" << item.base_sat.clk
			<< " tx_sow=" << item.base_sat.time.sec
			<< '\n';

		out << "ROVER_SAT "
			<< "XYZ_m="
			<< item.rover_sat.pos[0] << ','
			<< item.rover_sat.pos[1] << ','
			<< item.rover_sat.pos[2]
			<< " clock_s=" << item.rover_sat.clk
			<< " tx_sow=" << item.rover_sat.time.sec
			<< '\n';
	}

	out << "END_EPOCH\n";
	return out.good();
}

bool WriteStationObservations(ostream& out, const obs_t& obs)
{
	//判断有效性
	if (!out || obs.n <= 0 || obs.n > MAXOBS)
		return false;

	out << setprecision(15);

	const gtime_t& time = obs.data[0].time;
	out << "EPOCH week=" << time.week
		<< " sow=" << time.sec
		<< " satellites=" << obs.n << '\n';

	for (int i = 0; i < obs.n; ++i)
	{
		const obsd_t& item = obs.data[i];

		// 导出当前解码器保存的全部频点
		for (int f = 0; f < NFREQ + NEXOBS; ++f)
		{
			// 没有观测值，不写入
			if (item.code[f] == CODE_NONE &&
				item.P[f] == 0.0 &&
				item.L[f] == 0.0 &&
				item.D[f] == 0.0f &&
				item.SNR[f] == 0.0 &&
				item.LLI[f] == 0)
			{
				continue;
			}

			out << sat2id(item.sat)
				<< " slot=" << f
				<< " code=" << static_cast<int>(item.code[f])
				<< " P_m=" << item.P[f]
				<< " L_cycles=" << item.L[f]
				<< " D_Hz=" << item.D[f]
				<< " SNR_dBHz=" << item.SNR[f]
				<< " LLI=" << static_cast<int>(item.LLI[f])
				<< '\n';
		}
	}

	out << "END_EPOCH\n";
	return out.good();
}

bool ReadAndExportNext(FILE* fp,station_state_t& station,ostream& out,bool& output_failed)
{
	bool has_epoch = ReadNextObservationEpoch(fp, station);

	if (has_epoch &&!WriteStationObservations(out, station.raw.obs))
	{
		output_failed = true;
	}

	return has_epoch;
}