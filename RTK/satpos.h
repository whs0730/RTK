#pragma once
#include "obs.h"

// 卫星位置计算模块：根据广播星历计算卫星位置、速度、钟差和钟速

// GPS 广播星历计算，输入时间为 GPST
void CalculateGPS(gtime_t t, eph_t* eph, satpos_t* sat);

// BDS 广播星历计算，输入时间为 BDT
void CalculateBDS(gtime_t t, eph_t* eph, satpos_t* sat);
