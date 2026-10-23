#include "rtk_prepare.h"
#include "coordinate.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <memory>
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
		// 只有文件结束才置eof；坏报文跳过后仍继续查找下一条观测。
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
		// 解码成功但没有观测时继续读；SPP所需的卫星数量在处理阶段检查。
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
	// 必须与共视筛选、实时缓存配对使用同一个严格不等式。
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
	// 前一成功历元只提供迭代初值，不把其位置当成本历元有效解。
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
	// RTK几何依赖位置，不依赖多普勒速度，因此速度失败不否定本历元位置。
	int speed_nv = 0;
	station.speed_ok = SPP_Speed(obs.data, obs.n, &station.spp_sol, station.sats, &station.speed_sol, &speed_nv);
	if (!station.speed_ok || station.speed_sol.stat != 1)
	{
		station.speed_ok = false;
		station.speed_sol.stat = 0;
	}
	return true;//计算出位置就返回成功
}

// 视线由接收机指向卫星；WGS-84大地坐标用于转到测站局部水平系，输出为rad。
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
	// common保留原始观测下标，位置与卫星编号都有效后才能用于双差几何。
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
	// 解码槽0/1中的北斗观测必须分别为B1I/B3I，不能用其他信号套用这两条波长。
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

// 用两站接收时刻距toe的较大值衡量星历年龄，选此值最小的一套；同龄时优先基站。
// 卫星状态计算函数继续检查7200 s有效期，不把共同星历写回任一站的SPP缓存。
static const eph_t* SelectRtkEphemeris(const nav_t& base_nav, const nav_t& rover_nav,
	int sat, const gtime_t& base_time, const gtime_t& rover_time)
{
	const eph_t* selected = nullptr;
	double best_age = 0.0;
	for (const nav_t* nav : { &base_nav, &rover_nav })
	{
		for (const auto& eph : nav->eph)
		{
			if (eph.sat != sat || !has_valid_eph(&eph) || !isfinite(eph.sqrtA)) continue;
			const double age = max(fabs(timediff(base_time, eph.toe)), fabs(timediff(rover_time, eph.toe)));
			if (!isfinite(age)) continue;
			if (selected == nullptr || age < best_age)
			{
				selected = &eph;
				best_age = age;
			}
		}
	}
	return selected;
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
	// 两站广播星历报文可能不在同一秒到达。RTK必须使用同一套星历，
	// 再按各站观测分别求发射时刻和自转改正，不能直接复制另一站的卫星坐标。
	nav_t shared_nav{};
	for (int i = 0; i < base.raw.obs.n; ++i)
	{
		const int sat = base.raw.obs.data[i].sat;
		if (sat < 1 || sat > MAXSAT || rover_index[sat] < 0) continue;
		const eph_t* eph = SelectRtkEphemeris(base.raw.nav, rover.raw.nav, sat, base_time, rover_time);
		if (eph != nullptr) shared_nav.eph[sat - 1] = *eph;
	}
	satpos_t base_sats[MAXOBS]{};
	satpos_t rover_sats[MAXOBS]{};
	CaculateSatellitePositions(epoch.base_obs.data, epoch.base_obs.n, &shared_nav, base_sats);
	CaculateSatellitePositions(epoch.rover_obs.data, epoch.rover_obs.n, &shared_nav, rover_sats);

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
			base_sats[i], rover_sats[j],
			base.spp_sol, rover.spp_sol, config,
			base_azel, rover_azel))
		{
			continue;

		}
		// 保存共同星历重算的状态；下标仍指向未压缩的两站obs数组。
		common_sat_t& item = epoch.common[epoch.ncommon++];
		item.sat = sat;
		item.base_index = i;
		item.rover_index = j;
		item.base_sat = base_sats[i];
		item.rover_sat = rover_sats[j];
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
    const auto flags = out.flags();
    const auto precision = out.precision();
    out << fixed << setprecision(3) << ">RTK: " << epoch.time.week << ' ' << epoch.time.sec << "\n\n";
    out << "共视星：GPS=" << epoch.common_gps << " BDS=" << epoch.common_bds << " CommonSats=";
    for (int k = 0; k < epoch.ncommon; ++k) out << sat2id(epoch.common[k].sat) << ' ';
    out << "\n\n";
    const sol_t* solutions[] = {&epoch.base_spp, &epoch.rover_spp};
    const char* names[] = {"基准站SPP：\n", "流动站SPP：\n"};
    for (int station = 0; station < 2; ++station)
    {
        const auto& sol = *solutions[station];
        XYZ xyz{sol.XYZ[0], sol.XYZ[1], sol.XYZ[2]};
        unique_ptr<BLH> blh(XYZtoBLH(xyz, 6378137.0, 1.0 / 298.257223563));
        out << names[station] << "SPP: " << epoch.time.week << ' ' << setprecision(3) << epoch.time.sec
            << " X:" << setw(14) << sol.XYZ[0] << " Y:" << setw(14) << sol.XYZ[1] << " Z:" << setw(14) << sol.XYZ[2]
            << " B:" << setprecision(8) << setw(14) << blh->B * 180.0 / PI
            << " L:" << setw(14) << blh->L * 180.0 / PI
            << " H:" << setprecision(3) << setw(10) << blh->H
            << " PDOP:" << setw(8) << sol.pdop << " Sigma:" << setw(8) << sol.sigma0 << "\n\n";
    }
    out.flags(flags);
    out.precision(precision);
    return out.good();
}

bool WriteStationObservations(ostream& out, const obs_t& obs)
{
    if (!out || obs.n <= 0 || obs.n > MAXOBS) return false;
    const auto flags = out.flags();
    const auto precision = out.precision();
    const auto& time = obs.data[0].time;
    out << fixed << setprecision(3) << ">OBS: " << time.week << ' ' << time.sec << "\n\n";
    out << "Sat  Freq          P(m)            L(cycle)       D(Hz)   SNR(dB-Hz)\n";
    for (int i = 0; i < obs.n; ++i)
    {
        const auto& item = obs.data[i];
        for (int f = 0; f < NFREQ + NEXOBS; ++f)
        {
            if (item.code[f] == CODE_NONE && item.P[f] == 0.0 && item.L[f] == 0.0 &&
                item.D[f] == 0.0f && item.SNR[f] == 0.0 && item.LLI[f] == 0) continue;
            out << sat2id(item.sat) << "    f" << f + 1 << setprecision(3) << setw(17) << item.P[f]
                << setprecision(6) << setw(20) << item.L[f] << setprecision(3) << setw(12) << item.D[f]
                << setprecision(1) << setw(13) << item.SNR[f] << '\n';
        }
    }
    out << '\n';
    out.flags(flags);
    out.precision(precision);
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
