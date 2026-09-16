#include "timeTransform.h"
#include <cmath>

// 由 GPS 周和周内秒构造内部时间结构
gtime_t gpst2time(int week, double sec) {
	gtime_t t = { week, sec };
	return t;
}

// 计算两个 gtime_t 的时间差，单位：秒
double timediff(gtime_t t1, gtime_t t2) {
	return (t1.week - t2.week) * 604800.0 + (t1.sec - t2.sec);
}

// 从 gtime_t 中提取 GPS 周号，并返回周内秒
int time2gpst(gtime_t t, int* week) {
	if (week) *week = t.week;
	return t.sec;
}

// 调整过小的 GPS 周号，处理周号翻转后的旧数据。
int adjgpsweek(int week) {
	if (week < 300) week += 2048;
	return week;
}

// 公历时间转换为 MJD 时间
MJDTime* CommonTimeToMJDTime(CommonTime& common_time)
{
	long double JD = 0.0;
	MJDTime* MJD = new MJDTime();
	unsigned short M = common_time.Month;
	unsigned short Y = common_time.Year;
	unsigned short D = common_time.Day;
	unsigned short Hour = common_time.Hour;
	unsigned short Minute = common_time.Minute;
	double Second = common_time.Second;
	double UT = Hour + Minute / 60.0 + Second / 3600.0;
	unsigned short y;
	unsigned short m;
	if (M <= 2) {
		y = Y - 1;
		m = M + 12;
	}
	else {
		y = Y;
		m = M;
	}
	JD = (int)(365.25 * y) + (int)(30.6001 * (m + 1)) + D + UT / 24 + 1720981.5;
	MJD->Days = int(JD - 2400000.5);
	MJD->FracDay = JD - 2400000.5 - MJD->Days;
	return MJD;
}

// MJD 时间转换为公历时间
CommonTime* MJDToCommonTime(MJDTime& MJD)
{
	CommonTime* ct = new CommonTime();

	double jd = MJD.Days + MJD.FracDay + 2400000.5;
	int Z = (int)floor(jd + 0.5);
	double F = (jd + 0.5) - Z;

	int A = Z;
	int alpha = (int)((A - 1867216.25) / 36524.25);
	A = A + 1 + alpha - alpha / 4;

	int B = A + 1524;
	int C = (int)((B - 122.1) / 365.25);
	int D = (int)(365.25 * C);
	int E = (int)((B - D) / 30.6001);

	double day = B - D - (int)(30.6001 * E) + F;

	ct->Day = (int)day;

	if (E < 14) ct->Month = E - 1;
	else        ct->Month = E - 13;

	if (ct->Month > 2) ct->Year = C - 4716;
	else               ct->Year = C - 4715;

	double frac = day - ct->Day;
	double sod = frac * 86400.0;

	ct->Hour = (unsigned short)(sod / 3600.0);
	sod -= ct->Hour * 3600.0;
	ct->Minute = (unsigned short)(sod / 60.0);
	sod -= ct->Minute * 60.0;
	ct->Second = sod;

	return ct;
}

// MJD 时间转换为 GPS 周和周内秒
GPSTime* MJDToGPSTime(MJDTime& MJD) {
	GPSTime* GPS_Time = new GPSTime();
	GPS_Time->Week = (int)((MJD.Days + MJD.FracDay - 44244) / 7);
	GPS_Time->SecOfWeek = (MJD.Days + MJD.FracDay - 44244 - GPS_Time->Week * 7) * 86400.0;
	return GPS_Time;
}

// GPS 周和周内秒转换为 MJD 时间
MJDTime* GPSTimeToMJD(GPSTime& GPSTime) {
	MJDTime* MJD_Time = new MJDTime();
	double MJD = 44244 + GPSTime.Week * 7 + GPSTime.SecOfWeek / 86400.0;
	MJD_Time->Days = (int)(MJD);
	MJD_Time->FracDay = MJD - MJD_Time->Days;
	return MJD_Time;
}

// GPS 时间转换为公历时间
CommonTime* GPSTimeToCommonTime(GPSTime& GPSTime) {
	MJDTime* MJDTime = GPSTimeToMJD(GPSTime);
	return MJDToCommonTime(*MJDTime);
}

// 公历时间转换为 GPS 时间
GPSTime* CommonTimeToGPSTime(CommonTime& CommonTime) {
	MJDTime* MJDTime = CommonTimeToMJDTime(CommonTime);
	return MJDToGPSTime(*MJDTime);
}

// 由BDS 周和周内秒构造时间；内部周号先平移到 GPS 周号体系
gtime_t bdt2time(int week, double sec)
{
    gtime_t t;
    t.week = week + BDT_GPS_WEEK_OFFSET;
    t.sec = sec;
    return t;
}

// 给 gtime_t 加秒数，并处理跨周
gtime_t timeadd(gtime_t t, double dt)
{
    t.sec += dt;

    while (t.sec < 0.0) {
        t.sec += 604800.0;
        t.week--;
    }
    while (t.sec >= 604800.0) {
        t.sec -= 604800.0;
        t.week++;
    }
    return t;
}

// BDT 转 GPST：同一时刻 GPST 比 BDT 快 14 秒
gtime_t bdt2gpst(gtime_t t)
{
    return timeadd(t, 14.0);
}

// GPST 转 BDT：同一时刻 BDT 比 GPST 慢 14 秒
gtime_t gpst2bdt(gtime_t t)
{
    return timeadd(t, -14.0);
}
