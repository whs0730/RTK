#include "satpos.h"
#include <cmath>
#include "obs.h"
#include "matrix.h"
#include "define.h"

#define pi 3.1415926535898  
#define F -4.442807633e-10  //相对论改正常数
#define mu_gps 3.986005e14 //GPS 使用的地球引力常数
#define mu_bds 3.986004418e14 //BDS 使用的地球引力常数

// 根据 GPS 广播星历计算卫星 ECEF 位置、速度、钟差和钟速。
void CalculateGPS(gtime_t t, eph_t* eph, satpos_t* sat) {
	// GPS 星历和观测历元均按 GPST 组织，直接用 GPST 求 tk。
	double e = eph->e; //轨道偏心率
	double n0 = sqrt(mu_gps) / pow(eph->sqrtA, 3); //未改正的平均角速度
	double tk = timediff(t, eph->toe);//当前计算时刻 t 相对于星历参考时刻 toe 的时间差
	if (tk > 302400) {
		tk = tk - 604800;
	}
	else if (tk < -302400) {
		tk = tk + 604800;
	}
	else {
		tk = tk;
	}

	// 平近点角和偏近点角
	double n = n0 + eph->deln; //改正后的平均角速度
	double Mk = eph->M0 + n * tk; //平近点角
	//开普勒方程，求偏近点角 Ek
	double Ek = Mk;
	for (int iter = 0; iter < 30; iter++) {
		double E_old = Ek;
		Ek = Mk + e * sin(E_old);

		if (!isfinite(Ek)) {
			sat->pos[0] = sat->pos[1] = sat->pos[2] = 0.0;
			sat->vel[0] = sat->vel[1] = sat->vel[2] = 0.0;
			sat->clk = sat->dclk = 0.0;
			return;
		}

		if (fabs(Ek - E_old) < 1e-12) {
			break;
		}
		//30次未收敛清零返回
		if (iter == 29) {
			sat->pos[0] = sat->pos[1] = sat->pos[2] = 0.0;
			sat->vel[0] = sat->vel[1] = sat->vel[2] = 0.0;
			sat->clk = sat->dclk = 0.0;
			return;
		}
	}

	// 真近点角、升交角距和二阶调和改正。
	double fenZi = sqrt(1 - pow(e, 2)) * sin(Ek) / (1 - e * cos(Ek));
	double fenMu = (cos(Ek) - e) / (1 - e * cos(Ek));
	double vk = atan2(fenZi, fenMu); //真近点角
	double PHIk = vk + eph->omg;     //未改正的升交距角
	double phiuk = eph->cus * sin(2 * PHIk) + eph->cuc * cos(2 * PHIk); // 升交距角改正
	double phirk = eph->crs * sin(2 * PHIk) + eph->crc * cos(2 * PHIk);  // 轨道半径改正
	double phiik = eph->cis * sin(2 * PHIk) + eph->cic * cos(2 * PHIk);	  // 轨道倾角改正

	// 改正后的升交角距、向径和轨道倾角
	double uk = PHIk + phiuk;
	double rk = eph->A * (1 - e * cos(Ek)) + phirk;
	double ik = eph->i0 + phiik + eph->idot * tk;

	// 轨道平面坐标转 ECEF
	double xk_prime = rk * cos(uk);
	double yk_prime = rk * sin(uk);
	double OMGk = eph->OMG0 + (eph->OMGd - OMGed_GPS) * tk - OMGed_GPS * (eph->toe.sec);//升交点经度
	sat->pos[0] = xk_prime * cos(OMGk) - yk_prime * cos(ik) * sin(OMGk);
	sat->pos[1] = xk_prime * sin(OMGk) + yk_prime * cos(ik) * cos(OMGk);
	sat->pos[2] = yk_prime * sin(ik);

	// 对上述轨道参数求导，得到卫星速度
	double Ekd = n / (1 - e * cos(Ek));
	double PHIkd = sqrt((1 + e) / (1 - e)) * pow(cos(vk / 2) / cos(Ek / 2), 2) * Ekd;
	double ukd = 2 * (eph->cus * cos(2 * PHIk) - eph->cuc * sin(2 * PHIk)) * PHIkd + PHIkd;
	double rkd = eph->A * e * sin(Ek) * Ekd + 2 * (eph->crs * cos(2 * PHIk) - eph->crc * sin(2 * PHIk)) * PHIkd;
	double ikd = eph->idot + 2 * (eph->cis * cos(2 * PHIk) - eph->cic * sin(2 * PHIk)) * PHIkd;
	double OMGkd = eph->OMGd - OMGed_GPS;
	Matrix Rd = {
		{cos(OMGk), -sin(OMGk) * cos(ik), -(xk_prime * sin(OMGk) + yk_prime * cos(OMGk) * cos(ik)), yk_prime * sin(OMGk) * sin(ik)},
		{sin(OMGk), cos(OMGk) * cos(ik), xk_prime * cos(OMGk) - yk_prime * sin(OMGk) * cos(ik), yk_prime * cos(OMGk) * sin(ik)},
		{0, sin(ik), 0, yk_prime * cos(ik)}
	};
	double xk_primed = rkd * cos(uk) - rk * ukd * sin(uk);
	double yk_primed = rkd * sin(uk) + rk * ukd * cos(uk);
	Matrix M = {
		{xk_primed},
		{yk_primed},
		{OMGkd},
		{ikd}
	};
	Matrix speed = mul(Rd, M);
	sat->vel[0] = speed[0][0];
	sat->vel[1] = speed[1][0];
	sat->vel[2] = speed[2][0];

	// 卫星钟差和钟速。只改相对论项，不处理 TGD
	double delt = timediff(t, eph->toc);
	double deltr = F * e * eph->sqrtA * sin(Ek);//相对论改正
	sat->clk = eph->af0 + eph->af1 * delt + eph->af2 * delt * delt + deltr;//钟差
	double deltrd = F * e * eph->sqrtA * cos(Ek) * Ekd;
	sat->dclk = eph->af1 + 2.0 * eph->af2 * delt + deltrd;//钟速
}

// 根据 BDS 广播星历计算卫星 ECEF 位置、速度、钟差和钟速
void CalculateBDS(gtime_t t, eph_t* eph, satpos_t* sat) {
	// BDS 星历的 toe/toc 解码后已转换为 GPST，便于和观测历元求差
	// eph->toes 保留 BDT 周内秒，用于 BDS 升交点经度中的地球自转项
	double e = eph->e;
	double n0 = sqrt(mu_bds) / pow(eph->sqrtA, 3);
	double tk = timediff(t, eph->toe);
	if (tk > 302400) {
		tk = tk - 604800;
	}
	else if (tk < -302400) {
		tk = tk + 604800;
	}
	else {
		tk = tk;
	}

	double n = n0 + eph->deln;
	double Mk = eph->M0 + n * tk;
	double Ek = Mk;
	for (int iter = 0; iter < 30; iter++) {
		double E_old = Ek;
		Ek = Mk + e * sin(E_old);

		if (!isfinite(Ek)) {
			sat->pos[0] = sat->pos[1] = sat->pos[2] = 0.0;
			sat->vel[0] = sat->vel[1] = sat->vel[2] = 0.0;
			sat->clk = sat->dclk = 0.0;
			return;
		}

		if (fabs(Ek - E_old) < 1e-12) {
			break;
		}

		if (iter == 29) {
			sat->pos[0] = sat->pos[1] = sat->pos[2] = 0.0;
			sat->vel[0] = sat->vel[1] = sat->vel[2] = 0.0;
			sat->clk = sat->dclk = 0.0;
			return;
		}
	}

	// 真近点角、升交角距和二阶调和改正
	double fenZi = sqrt(1 - pow(e, 2)) * sin(Ek) / (1 - e * cos(Ek));
	double fenMu = (cos(Ek) - e) / (1 - e * cos(Ek));
	double vk = atan2(fenZi, fenMu);
	double PHIk = vk + eph->omg;
	double phiuk = eph->cus * sin(2 * PHIk) + eph->cuc * cos(2 * PHIk);
	double phirk = eph->crs * sin(2 * PHIk) + eph->crc * cos(2 * PHIk);
	double phiik = eph->cis * sin(2 * PHIk) + eph->cic * cos(2 * PHIk);

	double uk = PHIk + phiuk;
	double rk = eph->A * (1 - e * cos(Ek)) + phirk;
	double ik = eph->i0 + phiik + eph->idot * tk;
	double xk_prime = rk * cos(uk);
	double yk_prime = rk * sin(uk);

	// 速度计算所需的轨道参数导数
	double Ekd = n / (1 - e * cos(Ek));
	double vkd = sqrt(1 - e * e) * Ekd / (1 - e * cos(Ek));
	double PHIkd = vkd;
	double ukd = 2 * (eph->cus * cos(2 * PHIk) - eph->cuc * sin(2 * PHIk)) * PHIkd + PHIkd;
	double rkd = eph->A * e * sin(Ek) * Ekd + 2 * (eph->crs * cos(2 * PHIk) - eph->crc * sin(2 * PHIk)) * PHIkd;
	double ikd = eph->idot + 2 * (eph->cis * cos(2 * PHIk) - eph->cic * sin(2 * PHIk)) * PHIkd;

	// BDS MEO/IGSO 卫星
	if (eph->prn > 5 && eph->prn < 59) {
		double OMGk = eph->OMG0 + (eph->OMGd - OMGed_BDS) * tk - OMGed_BDS * (eph->toes);
		double Xk = xk_prime * cos(OMGk) - yk_prime * cos(ik) * sin(OMGk);
		double Yk = xk_prime * sin(OMGk) + yk_prime * cos(ik) * cos(OMGk);
		double Zk = yk_prime * sin(ik);
		sat->pos[0] = Xk;
		sat->pos[1] = Yk;
		sat->pos[2] = Zk;

		double OMGkd = eph->OMGd - OMGed_BDS;
		double xkd = rkd * cos(uk) - rk * ukd * sin(uk);
		double ykd = rkd * sin(uk) + rk * ukd * cos(uk);
		sat->vel[0] = -Yk * OMGkd - (ykd * cos(ik) - Zk * ikd) * sin(OMGk) + xkd * cos(OMGk);
		sat->vel[1] = Xk * OMGkd + (ykd * cos(ik) - Zk * ikd) * cos(OMGk) + xkd * sin(OMGk);
		sat->vel[2] = ykd * sin(ik) + yk_prime * ikd * cos(ik);
	}
	else {
		// BDS GEO 卫星需要从轨道坐标先转到倾斜地固系，再做 5 度倾角旋转。
		double OMGk = eph->OMG0 + eph->OMGd * tk - OMGed_BDS * (eph->toes);
		double XGK = xk_prime * cos(OMGk) - yk_prime * cos(ik) * sin(OMGk);
		double YGK = xk_prime * sin(OMGk) + yk_prime * cos(ik) * cos(OMGk);
		double ZGK = yk_prime * sin(ik);
		double phi1 = OMGed_BDS * tk;
		double phi2 = -5 * pi / 180.0;
		Matrix Rz = {
			{cos(phi1), sin(phi1), 0},
			{-sin(phi1), cos(phi1), 0},
			{0, 0, 1}
		};
		Matrix Rx = {
			{1, 0, 0},
			{0, cos(phi2), sin(phi2)},
			{0, -sin(phi2), cos(phi2)}
		};
		Matrix XYZGK = {
			{XGK},
			{YGK},
			{ZGK}
		};
		Matrix XYZ = mul(mul(Rz, Rx), XYZGK);
		sat->pos[0] = XYZ[0][0];
		sat->pos[1] = XYZ[1][0];
		sat->pos[2] = XYZ[2][0];

		double OMGkd = eph->OMGd;
		double xkd = rkd * cos(uk) - rk * ukd * sin(uk);
		double ykd = rkd * sin(uk) + rk * ukd * cos(uk);
		double XGKd = -YGK * OMGkd - (ykd * cos(ik) - ZGK * ikd) * sin(OMGk) + xkd * cos(OMGk);
		double YGKd = XGK * OMGkd + (ykd * cos(ik) - ZGK * ikd) * cos(OMGk) + xkd * sin(OMGk);
		double ZGKd = ykd * sin(ik) + yk_prime * ikd * cos(ik);
		Matrix Rzd = {
			{-sin(phi1) * OMGed_BDS, cos(phi1) * OMGed_BDS, 0},
			{-cos(phi1) * OMGed_BDS, -sin(phi1) * OMGed_BDS, 0},
			{0, 0, 0}
		};
		Matrix XYZGKd = {
			{XGKd},
			{YGKd},
			{ZGKd}
		};
		Matrix XYZkd = add(mul(mul(Rzd, Rx), XYZGK), mul(mul(Rz, Rx), XYZGKd));
		sat->vel[0] = XYZkd[0][0];
		sat->vel[1] = XYZkd[1][0];
		sat->vel[2] = XYZkd[2][0];
	}

	// 卫星钟差和钟速。TGD保留到伪距组合处处理，避免重复改正
	double delt = timediff(t, eph->toc);
	double deltr = F * e * eph->sqrtA * sin(Ek);
	sat->clk = eph->af0 + eph->af1 * delt + eph->af2 * delt * delt + deltr;
	double deltrd = F * e * eph->sqrtA * cos(Ek) * Ekd;
	sat->dclk = eph->af1 + 2.0 * eph->af2 * delt + deltrd;
}
