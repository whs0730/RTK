#pragma once

#include "coordinate.h"
#include "obs.h"

// Hopfield 对流层延迟模型，输入高程 H(m) 和高度角 E(rad),返回延迟(m)
double Hopfield_delTrop(double H, double E);

// Klobuchar 电离层延迟模型，返回斜路径延迟(m)
double Klobuchar(gtime_t t, const double* pos, const double* azel, double* ionpara);

// 卫星坐标取发射时刻ECEF，接收机坐标取接收时刻ECEF；返回含一阶地球自转改正的距离(m)。
double GeometricRange(const double receiver[3], const double satellite[3]);
