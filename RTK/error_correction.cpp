#include "error_correction.h"
#include <cmath>
// Hopfield 对流层延迟模型，输入高程 H(m) 和高度角 E(rad),返回延迟(m)
double Hopfield_delTrop(double H, double E) {
    if (H < -1000.0 || H > 50000.0 || E <= 0.0) {
        return 0.0;
    }
    // 标准气象元素。
    double H0 = 0.0;                 // 海平面高，单位：m
    double T0 = 15.0 + 273.16;       // 温度，单位：K
    double P0 = 1013.25;             // 气压，单位：mbar
    double RH0 = 0.5;                // 相对湿度
   
	double hd = 40136.0 + 148.72 * (T0 - 273.16);//分量等效高度 m
	double hw = 11000.0;//湿分量等效高度m
	double T = T0 - 0.0065 * (H - H0);//测站温度
	double P = P0 * pow((1 - 0.0000226 * (H - H0)), 5.225);//根据测站高度估算气压
	double RH = RH0 * exp(-0.0006396 * (H - H0));//据高度估算相对湿度
	double e = RH * exp(-37.2465 + 0.213166 * T - 0.000256908 * T * T);//水汽压
	double Kd = 155.2e-7 * P * (hd - H) / T;//干延迟部分的系数
	double Kw = 155.2e-7 * 4810 * e * (hw - H)/(T*T);//湿延迟部分的系数
    // Hopfield 映射函数里的 E 使用角度。
    double elev_deg =E * 180.0 / PI;

    double md = 1.0 / sin(sqrt(elev_deg * elev_deg + 6.25) * PI / 180.0);//分量映射函数
    double mw = 1.0 / sin(sqrt(elev_deg * elev_deg + 2.25) * PI / 180.0);//湿分量映射函数

    return Kd * md + Kw * mw;
}
// Klobuchar 电离层延迟改正，返回单位：m
//t当前观测时间，通常是 GPST
//pos 测站大地坐标，pos[0] = 纬度B，pos[1] = 经度L，单位 rad
//azel     卫星方位角和高度角，azel[0] = az，azel[1] = el，单位 rad
//ionpara  GPS广播电离层参数，共8个：alpha0~alpha3, beta0~beta3
//返回值   电离层延迟改正，单位 m
double Klobuchar(gtime_t t, const double* pos,const double* azel, double* ionpara)
{
    // 参数全 0 时视为无效，不做电离层改正。
    if (!pos || !azel || !ionpara) return 0.0;
    if (azel[1] <= 0.0) return 0.0;
    int valid = 0;
    for (int i = 0; i < 8; i++) {
        if (fabs(ionpara[i]) > 1E-30) {
            valid = 1;
            break;
        }
    }
    if (!valid) return 0.0;

    double psi, phi, lam, phi_m;
    double tt, f, amp, per, x;
    double az = azel[0]; //方位角 rad
    double el = azel[1];//高度角 rad

    // 大地坐标转换为半周单位 semicircle。
    double phi_u = pos[0] / PI;
    double lam_u = pos[1] / PI;
    double E = el / PI;
    double A = az / PI;

    psi = 0.0137 / (E + 0.11) - 0.022; // 电离层穿刺点偏移角。
    // 电离层穿刺点 IPP。
    phi = phi_u + psi * cos(az);
    if (phi > 0.416) phi = 0.416;
    else if (phi < -0.416) phi = -0.416;
    lam = lam_u + psi * sin(az) / cos(phi * PI);
    // 地磁纬度。
    phi_m = phi + 0.064 * cos((lam - 1.617) * PI);
    // 地方时。
    int week = 0;
    double sow = time2gpst(t, &week);
    tt = 43200.0 * lam + sow;
    while (tt >= 86400.0) tt -= 86400.0;
    while (tt < 0.0)      tt += 86400.0;
    // 振幅。
    amp = ionpara[0] + phi_m * (ionpara[1] + phi_m * (ionpara[2] + phi_m * ionpara[3]));
    if (amp < 0.0) amp = 0.0;
    // 周期。
    per = ionpara[4] + phi_m * (ionpara[5] + phi_m * (ionpara[6] + phi_m * ionpara[7]));
    if (per < 72000.0) per = 72000.0;
    // 延迟计算。
    x = 2.0 * PI * (tt - 50400.0) / per;
    // 斜路径映射因子。
    f = 1.0 + 16.0 * pow(0.53 - E, 3.0);
    double delay_sec;
    if (fabs(x) < 1.57) {
        // 白天。
        delay_sec = f * (5E-9 + amp * (1.0 - x * x / 2.0 + x * x * x * x / 24.0));
    }
    else {
        // 夜间。
        delay_sec = f * 5E-9;
    }

    return Clight * delay_sec; // 单位：m
}
// 地球自转改正：把发射时刻卫星坐标旋转到接收时刻对应的地固系
// tau 为信号传播时间，单位：s
void ApplyEarthRotationCorrection(satpos_t& sat, double tau)
{
    double theta = OMGed_GPS * tau;

    double x = sat.pos[0];
    double y = sat.pos[1];

    sat.pos[0] = cos(theta) * x + sin(theta) * y;
    sat.pos[1] = -sin(theta) * x + cos(theta) * y;
}
