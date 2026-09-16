#pragma once

#include "coordinate.h"
#include "obs.h"

// Hopfield 对流层延迟模型，输入高程 H(m) 和高度角 E(rad),返回延迟(m)
double Hopfield_delTrop(double H, double E);

// Klobuchar 电离层延迟模型，返回斜路径延迟(m)
double Klobuchar(gtime_t t, const double* pos, const double* azel, double* ionpara);

// 地球自转改正，把发射时刻卫星坐标旋转到接收时刻对应的地固系。
void ApplyEarthRotationCorrection(satpos_t& sat, double tau);
