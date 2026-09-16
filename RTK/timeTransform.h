#pragma once

// 公历时间。
struct CommonTime
{
    unsigned short Year = 0;
    unsigned short Month = 0;
    unsigned short Day = 0;
    unsigned short Hour = 0;
    unsigned short Minute = 0;
    double Second = 0.0;
};

// 简化 MJD 时间：整数天 + 日内小数。
struct MJDTime
{
    int Days;
    double FracDay;
    MJDTime()
    {
        Days = 0;
        FracDay = 0.0;
    }
};

// GPS 周和周内秒。
struct GPSTime
{
    unsigned short Week;
    double SecOfWeek;
    GPSTime(int week = 0, double secOfWeek = 0.0)
    {
        Week = static_cast<unsigned short>(week);
        SecOfWeek = secOfWeek;
    }
};

// 工程内部使用的轻量时间结构，week/sec 可表示 GPST 或 BDT
struct gtime_t
{
    int week;
    double sec;
};

// BDT 周号比 GPST 小 1356 周。
#define BDT_GPS_WEEK_OFFSET 1356

gtime_t gpst2time(int week, double sec);
double timediff(gtime_t t1, gtime_t t2);//时间差
int time2gpst(gtime_t t, int* week);
int adjgpsweek(int week);

// 公历、MJD、GPST 三类时间之间的转换函数。
MJDTime* CommonTimeToMJDTime(CommonTime& common_time);
CommonTime* MJDToCommonTime(MJDTime& MJD);
GPSTime* MJDToGPSTime(MJDTime& MJD);
MJDTime* GPSTimeToMJD(GPSTime& GPSTime);
CommonTime* GPSTimeToCommonTime(GPSTime& GPSTime);
GPSTime* CommonTimeToGPSTime(CommonTime& CommonTime);

gtime_t bdt2time(int week, double sec);
gtime_t timeadd(gtime_t t, double dt);//时间增加
gtime_t bdt2gpst(gtime_t t);
gtime_t gpst2bdt(gtime_t t);
