#pragma once
#define PI 3.1415926535898
struct XYZ
{
    double X = 0.0;
    double Y = 0.0;
    double Z = 0.0;
};
struct BLH
{
    double B = 0.0;
    double L = 0.0;
    double H = 0.0;
};
struct ENU
{
    double E = 0.0;
    double N = 0.0;
    double U = 0.0;
};
XYZ* BLHtoXYZ(BLH& blh, double a, double f);
double calculate_B(XYZ& xyz, double a, double f);
BLH* XYZtoBLH(XYZ& xyz, double a, double f);
ENU* XYZtoENU(BLH* blh, double* e);
// 由站点位置和视线单位向量计算方位角、高度角
double satazel(BLH* blh, double* e, double* azel);
