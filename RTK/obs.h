#pragma once

#include <cstdio>
#include <string>

#include "define.h"
#include "matrix.h"
#include "timeTransform.h"


// 单颗卫星在一个历元的观测值。
struct obsd_t {
    gtime_t time;                           // 观测时间
    int sat;                                // 内部卫星编号
    double L[NFREQ + NEXOBS];               // 载波相位，单位：周
    double P[NFREQ + NEXOBS];               // 伪距，单位：m
    float D[NFREQ + NEXOBS];                // 多普勒，单位：Hz
    double SNR[NFREQ + NEXOBS];             // 信噪比，单位：dB-Hz
    unsigned char LLI[NFREQ + NEXOBS];      // 周跳标记
    unsigned char code[NFREQ + NEXOBS];     // 观测码类型
};

// 一个历元的观测集合。
struct obs_t {
    int n;                  // 当前观测卫星数
    obsd_t data[MAXOBS];    // 观测数组
};

// 广播星历参数。
struct eph_t {
    int sat;        // 内部卫星编号
    int prn;        // 系统内 PRN
    int week;       // 周号
    double toes;    // 原始 toe 秒，BDS 为 BDT 周内秒
    double tocs;    // 原始 toc 秒
    gtime_t toe, toc;

    // 轨道参数
    double sqrtA;
    double A;       // 轨道长半轴
    double e;       // 偏心率
    double M0;      // 参考时刻平近点角
    double OMG0;    // 参考时刻升交点赤经
    double omg;     // 近地点角距
    double i0;      // 参考时刻轨道倾角

    double deln;    //平均角速度改正 Δn
    double OMGd;    // 升交点赤经变化率Ωdot
    double idot;     // 倾角变化率IDOT

    // 摄动改正参数
    double cuc, cus; //升交距角余弦和正弦改正项
    double crc, crs; //轨道半径余弦和正弦改正项
    double cic, cis; //轨道倾角余弦和正弦改正项

    // 卫星钟差和群延迟参数
    double af0, af1, af2; //卫星钟偏、钟速、钟漂
    double tgd[2]; //BDS群延迟改正
};

// 导航数据集合。
struct nav_t {
    eph_t eph[MAXSAT];
    double ion_gps[8]; //GPS Klobuchar 电离层模型参数
    double utc_gps[8]; 
    int glo_fcn[32];
};

// 卫星状态计算结果。
struct satpos_t {
    int sat;
    gtime_t time;      //发射时刻

    double pos[3];      // ECEF 坐标，单位：m
    double vel[3];      // ECEF 速度，单位：m/s
    double clk;         // 卫星钟差，单位：s
    double dclk;        // 卫星钟速，单位：s/s
};

// 单点定位结果。
struct sol_t {
    gtime_t time;       // 当前历元时间
    double XYZ[3];      // 接收机 ECEF 坐标，单位：m
    double dtr[2];      // 接收机钟差，单位：m；dtr[0] GPS，dtr[1] BDS
    Matrix Q;           // 协方差矩阵

    double pdop;        // PDOP
    double sigma0;      // 单位权中误差
    int ns;             // 参与解算卫星数
    int stat;           // 0 失败，1 成功
};

// 单点测速结果。
struct solvel_t {
    gtime_t time;       // 当前历元时间
    double V[3];        // 接收机 ECEF 速度，单位：m/s
    double dtrd;        // 接收机钟速，单位：m/s
    Matrix Q;           // 协方差矩阵
    double sigma0;      // 单位权中误差
    int ns;             // 参与测速卫星数
    int stat;           // 0 失败，1 成功
};
//根据卫星系统和prn计算卫星号
int satno(int sys, int prn);
// 根据内部 sat 编号反推出卫星系统和系统内 PRN。
int satsys(int sat, int* prn);
// 将内部卫星编号转换为 Gxx/Cxx 等输出名称。
string sat2id(int sat); 
// 判断卫星位置是否有效
bool IsValidSatpos(satpos_t* sat);
