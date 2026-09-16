#pragma once

#include "obs.h"

// 计算双频无电离层组合伪距。GPS 使用 L1/L2，BDS 使用 B1I/B3I，并在组合处处理 BDS TGD
double GetPIF(obsd_t* obs, const eph_t* eph = nullptr);
// 为 SPP 查找可用星历；无效星历返回 nullptr
const eph_t* GetSppEph(const nav_t* nav, int sat);
// SPP筛选：系统、星历、伪距、卫星位置、信噪比、高度角。
bool PassSppBasicCheck(const obsd_t& obs,const satpos_t& sat,const eph_t* eph,const double* rec_xyz,bool check_elev,double* pif_out = nullptr);

// 单点定位主入口，输出接收机坐标、钟差、PDOP 和单位权中误差。
bool SPP(obsd_t* obs,int n,const nav_t* nav,sol_t* sol,satpos_t* sat,int* nv,const sol_t* init_sol = nullptr);

// 单点测速主入口，使用多普勒观测估计接收机速度和钟速。
bool SPP_Speed(obsd_t* obs,int n,sol_t* sol,satpos_t* sat,solvel_t* vsol,int* nv);
