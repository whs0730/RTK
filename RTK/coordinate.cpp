#include "coordinate.h"
#include <cmath>

// 大地坐标转 ECEF 坐标，a 为椭球长半轴，f 为扁率。
XYZ* BLHtoXYZ(BLH& blh, double a, double f) {
	double e_2 = 2 * f - f * f;
	double sinb = sin(blh.B);
	double cosb = cos(blh.B);
	double sinl = sin(blh.L);
	double cosl = cos(blh.L);
	double N = a / sqrt(1 - e_2 * sinb * sinb);
	XYZ* xyz = new XYZ();
	xyz->X = (N + blh.H) * cosb * cosl;
	xyz->Y = (N + blh.H) * cosb * sinl;
	xyz->Z = (N * (1 - e_2) + blh.H) * sinb;
	return xyz;
}

// 迭代求解大地纬度 B。
double calculate_B(XYZ& xyz, double a, double f) {
	double e_2 = 2 * f - f * f;
	double p = sqrt(xyz.X * xyz.X + xyz.Y * xyz.Y);
	double tanb_0 = xyz.Z / p;
	while (true)
	{
		double tanb = (xyz.Z + (a * e_2 * tanb_0) / (sqrt(1 + (1 - e_2) * tanb_0 * tanb_0))) / p;
		if (fabs(tanb - tanb_0) < 5e-10) {
			return atan(tanb);
		}
		tanb_0 = tanb;
	}
}

// ECEF 坐标转大地坐标。
BLH* XYZtoBLH(XYZ& xyz, double a, double f) {
	double e_2 = 2 * f - f * f;
	BLH* blh = new BLH();
	blh->B = calculate_B(xyz, a, f);
	blh->L = atan2(xyz.Y, xyz.X);
	double N = a / sqrt(1 - e_2 * sin(blh->B) * sin(blh->B));
	blh->H = sqrt(xyz.X * xyz.X + xyz.Y * xyz.Y) / cos(blh->B) - N;
	return blh;
}

// 将 ECEF 坐标差向量旋转到当地 ENU 坐标系。
ENU* XYZtoENU(BLH* blh, double* e)
{
	ENU* enu = new ENU;

	double sinp = sin(blh->B), cosp = cos(blh->B);
	double sinl = sin(blh->L), cosl = cos(blh->L);
	double E[9];

	E[0] = -sinl;           E[1] = cosl;           E[2] = 0.0;
	E[3] = -sinp * cosl;    E[4] = -sinp * sinl;   E[5] = cosp;
	E[6] = cosp * cosl;     E[7] = cosp * sinl;    E[8] = sinp;

	enu->E = E[0] * e[0] + E[1] * e[1] + E[2] * e[2];
	enu->N = E[3] * e[0] + E[4] * e[1] + E[5] * e[2];
	enu->U = E[6] * e[0] + E[7] * e[1] + E[8] * e[2];

	return enu;
}

// 由站点位置和视线单位向量计算方位角、高度角。
double satazel(BLH* blh, double* e, double* azel)
{
	ENU* enu = XYZtoENU(blh, e);
	azel[0] = atan2(enu->E, enu->N);
	if (azel[0] < 0.0) azel[0] += 2.0 * PI;
	azel[1] = asin(enu->U);
	delete enu;
	return azel[1];
}
