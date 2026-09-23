#include "station_satpos.h"
#include "spp.h"
#include "satpos.h"
#include "decode.h"
#include "error_correction.h"

#include <cmath>

static const double MAX_EPH_AGE_SEC = 7200.0;
// 按卫星系统调用广播星历计算函数。
static bool CalculateBroadcastSatState(gtime_t calc_time, eph_t* eph, int sat_no, satpos_t& sat)
{
    if (!has_valid_eph(eph))
    {
        return false;
    }

    satpos_t tmp{};
    tmp.sat = sat_no;
    tmp.time = calc_time;

    int prn = 0;
    int sys = satsys(sat_no, &prn);

    if (sys == SYS_GPS)
    {
        CalculateGPS(calc_time, eph, &tmp);
    }
    else if (sys == SYS_CMP)
    {
        CalculateBDS(calc_time, eph, &tmp);
    }
    else
    {
        return false;
    }

    if (tmp.pos[0] == 0.0 && tmp.pos[1] == 0.0 && tmp.pos[2] == 0.0)
    {
        return false;
    }

    sat = tmp;
    return true;
}

// 计算信号发射时刻的卫星状态，并把卫星坐标旋转到接收时刻对应的地固系。
// 步骤：
// 1. 用双频IF伪距 PIF/c 得到近似传播时间；
// 2. 得到近似发射时刻，并计算一次卫星钟差；
// 3. 用卫星钟差修正发射时刻，再重新计算卫星位置、速度、钟差、钟速；
// 4. 按传播时间做地球自转改正。
static bool CaculateSatStateAtTransmitTime(const obsd_t& obs, eph_t* eph, satpos_t& sat)
{
    //双频无电离层组合伪距
    double pif = GetPIF((obsd_t*)&obs, eph);
    if (pif <= 0.0)
    {
        return false;
    }
    //近似传播时间
    double tau0 = pif / Clight;
    //近似发射时刻
    gtime_t tx0 = timeadd(obs.time, -tau0);

    satpos_t sat0{};
    if (!CalculateBroadcastSatState(tx0, eph, obs.sat, sat0))
    {
        return false;
    }

    // 卫星钟差为“卫星钟 - 系统时”，发射时刻需要再扣除该钟差。
    gtime_t tx = timeadd(tx0, -sat0.clk);

    if (!CalculateBroadcastSatState(tx, eph, obs.sat, sat))
    {
        return false;
    }
    //修正后的发射时刻重新计算卫星位置、速度、钟差、钟速
    double tau = timediff(obs.time, tx);
    ApplyEarthRotationCorrection(sat, tau);
    sat.time = tx;

    return true;
}
// 对当前历元每颗卫星计算位置、速度、钟差和钟速。validSat统计的是“有可用星历且成功算出发射时刻卫星状态”的GPS/BDS卫星数，
int CaculateSatellitePositions(obsd_t* obs, int n, nav_t* nav, satpos_t* sats)
{
    int validSat = 0;

    for (int i = 0; i < n; i++)
    {
        //取观测值和对应星历
        obsd_t* ob = &obs[i];
        eph_t* eph = find_eph(nav, ob->sat);

        if (!has_valid_eph(eph))
        {
            continue;
        }
        //测时刻和星历参考时刻 toe 相差超过 7200 秒跳过，防止过期
        if (fabs(timediff(ob->time, eph->toe)) > MAX_EPH_AGE_SEC)
        {
            continue;
        }

        int prn = 0;
        int sys = satsys(ob->sat, &prn);
        //只考虑GPS和BDS
        if (sys != SYS_GPS && sys != SYS_CMP)
        {
            continue;
        }
        satpos_t sat{};
        if (!CaculateSatStateAtTransmitTime(*ob, eph, sat))
        {
            continue;
        }
        sats[i] = sat;
        validSat++;
    }

    return validSat;
}