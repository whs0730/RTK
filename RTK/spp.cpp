#include "obs.h"
#include <cmath>
#include "matrix.h"
#include "define.h"
#include "error_correction.h"
// 计算双频无电离层组合伪距
double GetPIF(obsd_t* obs,const eph_t* eph=nullptr)
{
    if (obs == nullptr)
    {
        return 0.0;
    }
    int prn = 0;
    int sys = satsys(obs->sat, &prn);
    if (sys == SYS_GPS)
    {
        // GPS L1C + L2P(Y)
        if (obs->P[0] == 0.0 || obs->P[1] == 0.0)
        {
            return 0.0;
        }

        double f1 = FREQ_GPS_L1;
        double f2 = FREQ_GPS_L2;

        return (f1 * f1 * obs->P[0] - f2 * f2 * obs->P[1]) /
            (f1 * f1 - f2 * f2);
    }
    else if (sys == SYS_CMP)
    {
        // BDS B1I + B3I
        if (obs->P[0] == 0.0 || obs->P[1] == 0.0)
        {
            return 0.0;
        }

        double f1 = FREQ_BDS_B1;
        double f3 = FREQ_BDS_B3;
        double P1 = obs->P[0];  // B1I
        double P3 = obs->P[1];  // B3I
        // BDS B1I 要改 TGD，B3I 不改，TGD 单位是秒，乘光速变成米；这里处理后，后续最小二乘直接使用IF伪距。
        if (eph != nullptr)
        {
            P1 -= Clight * eph->tgd[0];
        }

        return (f1 * f1 * P1 - f3 * f3 * P3) /
            (f1 * f1 - f3 * f3);
    }

    return 0.0;
}
//判断一颗卫星的双频 SNR 是否合格,低于 30 dB-Hz 的观测不参与解算
static bool PassSnrCheck(const obsd_t& obs)
{
    int prn = 0;
    int sys = satsys(obs.sat, &prn);

    double snr1 = 0.0;
    double snr2 = 0.0;

    if (sys == SYS_GPS)
    {
        // GPS L1C + L2P(Y)
        snr1 =obs.SNR[0];
        snr2 =obs.SNR[1];
    }
    else if (sys == SYS_CMP)
    {
        // BDS B1I + B3I
        snr1 = obs.SNR[0];
        snr2 =obs.SNR[1];
    }
    else
    {
        return false;
    }
    if (snr1 > 0.0 && snr1 < 30.0)
    {
        return false;
    }

    if (snr2 > 0.0 && snr2 < 30.0)
    {
        return false;
    }

    return true;
}
//判断卫星位置是否有效
static bool IsValidSatposSPP(const satpos_t& sat)
{
    return !(sat.pos[0] == 0.0 &&
        sat.pos[1] == 0.0 &&
        sat.pos[2] == 0.0);
}
// 判断广播星历是否有效:星历是否存在、卫星号是否有效、轨道长半轴是否有效。
static bool HasValidEphSPP(const eph_t* eph)
{
    if (eph == nullptr)
    {
        return false;
    }

    if (eph->sat <= 0)
    {
        return false;
    }

    if (eph->sqrtA <= 0.0)
    {
        return false;
    }

    return true;
}
// 根据内部卫星号从导航数据中取星历
const eph_t* GetSppEph(const nav_t* nav, int sat)
{
    if (nav == nullptr || sat < 1 || sat > MAXSAT)
    {
        return nullptr;
    }

    const eph_t* eph = &nav->eph[sat - 1];
    return HasValidEphSPP(eph) ? eph : nullptr;
}

// 计算当前接收机位置下的卫星高度角，单位为弧度。
// 高度角筛选、对流层改正和最终用星统计复用
static bool PassElevationCheck(const satpos_t& sat,const double rec_xyz[3],double& elev_rad)
{
    elev_rad = 0.0;
    //因为首次迭代位置可能是0
    if (rec_xyz == nullptr)
    {
        return false;
    }
    //判断卫星是否有效
    if (!IsValidSatposSPP(sat))
    {
        return false;
    }
    //接收机数组形式的 XYZ 转成你自定义的 XYZ 结构体

    XYZ recXYZ;
    recXYZ.X = rec_xyz[0];
    recXYZ.Y = rec_xyz[1];
    recXYZ.Z = rec_xyz[2];
    //转为BLH

    BLH* recBLH = XYZtoBLH(recXYZ, 6378137.0, 1.0 / 298.257223563);
    //计算接收机到卫星的向量，以及几何距离 rho

    double dx = sat.pos[0] - rec_xyz[0];
    double dy = sat.pos[1] - rec_xyz[1];
    double dz = sat.pos[2] - rec_xyz[2];

    double rho = sqrt(dx * dx + dy * dy + dz * dz);

    if (rho < 1.0)
    {
        delete recBLH;
        return false;
    }
    //站星方向向量单位化，得到 LOS 单位向量

    double los[3];
    los[0] = dx / rho;
    los[1] = dy / rho;
    los[2] = dz / rho;
    //测站 BLH 和 LOS 计算方位角、高度角

    double azel[2] = { 0.0, 0.0 };
    elev_rad = satazel(recBLH, los, azel);

    delete recBLH;
    return true;
}


// =====================================================
// SPP统一质量控制：
// 1. GPS/BDS
// 2. 卫星位置有效
// 3. 双频PIF有效
// 4. SNR >= 30 dB-Hz
// 5.高度角 >= 10°
// =====================================================
bool PassSppBasicCheck(const obsd_t& obs,const satpos_t& sat,const eph_t* eph,const double* rec_xyz,bool check_elev,double* pif_out = nullptr)
{
    int prn = 0;
    int sys = satsys(obs.sat, &prn);

    // 当前SPP只解GPS和BDS；其他系统先不参与定位
    if (sys != SYS_GPS && sys != SYS_CMP)
    {
        return false;
    }

    // 没有卫星坐标时无法建立观测方程
    if (!IsValidSatposSPP(sat))
    {
        return false;
    }

    // 双频PIF有效
    double pif = GetPIF((obsd_t*)&obs, eph);
    if (pif == 0.0)
    {
        return false;
    }

    //SNR >= 30 dB-Hz
    if (!PassSnrCheck(obs))
    {
        return false;
    }
    //高度角 >= 10°
    // 初始迭代时还没有可靠接收机坐标，所以高度角检查可选,一旦有接收机坐标，就启用10度截止高度角
    if (check_elev && rec_xyz != nullptr)
    {
        double elev_rad = 0.0;
        if (!PassElevationCheck(sat, rec_xyz, elev_rad) ||elev_rad < 10.0 * PI / 180.0)
        {
            return false;
        }
    }
    //获取有效PIF
    if (pif_out != nullptr)
    {
        *pif_out = pif;
    }

    return true;
}
//利用枚举类型判断历元卫星种类数量：单系统和多系统
enum SppSolveMode
{
    SPP_MODE_GPS_ONLY,
    SPP_MODE_BDS_ONLY,
    SPP_MODE_GPS_BDS
};
//判断未知数数量，双系统返回 5 个未知数，单系统返回 4 个未知数
int SppUnknownCount(SppSolveMode mode)
{
    return mode == SPP_MODE_GPS_BDS ? 5 : 4;
}
//判断该卫星是GPS或者BDS，决定是否参与该模式下是否参与解算
bool IsSppSystemUsed(int sys, SppSolveMode mode)
{
    if (mode == SPP_MODE_GPS_BDS)
    {
        return sys == SYS_GPS || sys == SYS_CMP;
    }
    if (mode == SPP_MODE_GPS_ONLY)
    {
        return sys == SYS_GPS;
    }
    return sys == SYS_CMP;
}
//决定接收机钟差放在状态向量的第几列
int SppClockIndex(int sys, SppSolveMode mode)
{
    if (mode != SPP_MODE_GPS_BDS)
    {
        return 3;
    }
    if (sys == SYS_GPS)
    {
        return 3;
    }
    if (sys == SYS_CMP)
    {
        return 4;
    }
    return -1;
}

// 利用上一历元定位结果设置本历元迭代初值；首次定位或上一历元失败时仍从0开始
Matrix MakeInitialSppState(SppSolveMode mode, const sol_t* init_sol)
{
    Matrix X = zeros(SppUnknownCount(mode), 1);
    //没有上一历元解，或者上一历元定位失败，或者上一历元坐标全0，就返回全 0 初值
    if (init_sol == nullptr || init_sol->stat != 1)
    {
        return X;
    }
    if (init_sol->XYZ[0] == 0.0 && init_sol->XYZ[1] == 0.0 && init_sol->XYZ[2] == 0.0)
    {
        return X;
    }

    X[0][0] = init_sol->XYZ[0];
    X[1][0] = init_sol->XYZ[1];
    X[2][0] = init_sol->XYZ[2];
    //根据系统类型设置未知数
    if (mode == SPP_MODE_GPS_BDS)
    {
        X[3][0] = init_sol->dtr[0];
        X[4][0] = init_sol->dtr[1];
    }
    else if (mode == SPP_MODE_GPS_ONLY)
    {
        X[3][0] = init_sol->dtr[0];
    }
    else
    {
        X[3][0] = init_sol->dtr[1];
    }

    return X;
}
//#########################位置解算###################################
// 构建 SPP 线性化观测方程并求一次最小二乘改正
Matrix* LeastSquares(Matrix X, obsd_t* obs, satpos_t* sat, const nav_t* nav, int n, int nv, SppSolveMode mode)
{
    int nx = SppUnknownCount(mode);//nx 是未知数个数
    Matrix B = zeros(nv, nx);//B 是设计矩阵
    Matrix w = zeros(nv, 1);//w 是观测减计算向量
    int index = 0;//index 是当前已经写入第几颗可用卫星
    for (int i = 0; i < n; i++)
    {
        //获取第i颗卫星的位置
        satpos_t s = sat[i];
        // 获取第i颗卫星的星历，并交给PassSppBasicCheck()完成筛选
        const eph_t* eph = GetSppEph(nav, obs[i].sat);

        double pif_check = 0.0;
        if (!PassSppBasicCheck(obs[i], s, eph, nullptr, false, &pif_check))
        {
            continue;
        }

        int prn = 0;
        int sys = satsys(obs[i].sat, &prn);
        //当前解算模式，决定这颗卫星是否参与
        if (!IsSppSystemUsed(sys, mode))
        {
            continue;
        }
        //防止数组越界
        if (index >= nv)
        {
            return nullptr;
        }
        double Xs = s.pos[0], Ys = s.pos[1], Zs = s.pos[2];
        double delX = Xs - X[0][0], delY = Ys - X[1][0], delZ = Zs - X[2][0];
        //计算几何距离
        double p = sqrt(delX * delX + delY * delY + delZ * delZ);
        double Trop = 0.0;//对流层延迟初始为 0
        //如果当前接收机坐标不是全 0，就说明已经有近似位置，可以计算高度角和对流层改正
        if (!(X[0][0] == 0.0 && X[1][0] == 0.0 && X[2][0] == 0.0))
        {
            // 从第二轮迭代开始已有接收机近似坐标，可以计算高度角并进行10度截止
            double rec_xyz[3] = {
                X[0][0],
                X[1][0],
                X[2][0]
            };

            double elev_rad = 0.0;
            if (!PassElevationCheck(s, rec_xyz, elev_rad) ||
                elev_rad < 10.0 * PI / 180.0)
            {
                continue;
            }

            XYZ recXYZ;
            recXYZ.X = X[0][0];
            recXYZ.Y = X[1][0];
            recXYZ.Z = X[2][0];

            BLH* recBLH = XYZtoBLH(recXYZ, 6378137.0, 1.0 / 298.257223563);

            // 对流层改正
            Trop = Hopfield_delTrop(recBLH->H, elev_rad);

            delete recBLH;
        }
        //几何距离太小，说明计算异常
        if (p < 1.0)
        {
            return nullptr;
        }
        //设计矩阵
        double B1 = -delX / p, B2 = -delY / p, B3 = -delZ / p;
        int clk_index = SppClockIndex(sys, mode);
        if (clk_index < 0)
        {
            continue;
        }
        //前三列是坐标偏导数，钟差列填 1
        B[index][0] = B1;
        B[index][1] = B2;
        B[index][2] = B3;
        B[index][clk_index] = 1.0;
        double rec_clk = X[clk_index][0];//取当前近似的接收机钟差 m
        w[index][0] = pif_check - (p + rec_clk - Clight * s.clk + Trop);
        index++;
    }
    //有效观测数必须不少于未知数个数。双系统至少 5 颗，单系统至少 4 颗
    if (index < nx)
    {
        return nullptr;
    }
    //把矩阵裁剪到真实可用卫星数量。因为前面可能剔除了一些卫星，实际 index 可能小于预分配的 nv
    B.resize(index);
    w.resize(index);
    Matrix N, dx;
    N = mul(transpose(B), B);
    N = inverse(N);
    dx = mul(N, transpose(B));
    dx = mul(dx, w);
    Matrix *m = new Matrix[4];
    m[0] = add(dx, X);//更新后的状态 X+dx
    m[1] = B;//设计矩阵 B
    m[2] = w;//w 向量
    m[3] = dx;//本次改正数 dx
    return m;
}
// 迭代执行最小二乘，直到接收机位置改正量足够小或达到最大迭代次数。
Matrix* IterativeSolution(obsd_t* obs, satpos_t* sat, const nav_t* nav, int n, int nv, SppSolveMode mode, const sol_t* init_sol, int times = 1, double Threshold = 1e-6)
{
    //先用上一历元结果设置初值。如果没有上一历元，就用全 0
    Matrix X = MakeInitialSppState(mode, init_sol);
    Matrix *m = nullptr; // X,B,w
    while (times <= 10)
    {
        if (m != nullptr)
        {
            delete[] m;
            m = nullptr;
        }
        m = LeastSquares(X, obs, sat, nav,n, nv, mode);
        if (m == nullptr)
        {
            return nullptr;
        }
        Matrix dX = sub(m[0], X);
        //计算本轮坐标改正量的模长
        double dpos = sqrt(dX[0][0] * dX[0][0] + dX[1][0] * dX[1][0] + dX[2][0] * dX[2][0]);
        if (dpos < Threshold)
        {
            return m;
        }

        X = m[0];
        times++;
    }
    return m;
}
// 计算存储可用卫星的 IF 组合伪距，并统计 GPS/BDS 数量。
void PIF(obsd_t* obs,int n,satpos_t* sat,const nav_t* nav,int* nv,double* pif,int* a = nullptr,int* b = nullptr)
{
    int index = 0;
    //初始化总可用卫星数、GPS 数量、BDS 数量
    if (nv) *nv = 0;
    if (a)  *a = 0;
    if (b)  *b = 0;

    for (int i = 0; i < n; i++)
    {
        obsd_t o = obs[i];
        satpos_t s = sat[i];
        // PIF预筛也传入星历，使BDS TGD改正和后续最小二乘保持一致。
        const eph_t* eph = GetSppEph(nav, obs[i].sat);
        double pif_value = 0.0;

        // 高度角筛选会在最小二乘迭代中有近似坐标后再执行
        if (!PassSppBasicCheck(o, s, eph, nullptr, false, &pif_value))
        {
            continue;
        }

        int prn = 0;
        int sys = satsys(o.sat, &prn);

        if (sys == SYS_GPS)
        {
            if (a) (*a)++;
        }
        else if (sys == SYS_CMP)
        {
            if (b) (*b)++;
        }

        pif[index++] = pif_value;
    }

    if (nv) *nv = index;
}
// SPP 主函数：双系统估计 x/y/z/dtG/dtC，单系统估计 x/y/z/dt。
bool SPP(obsd_t *obs, int n, const nav_t *nav, sol_t *sol, satpos_t *sat, int *nv, const sol_t* init_sol)
{
    if (obs == nullptr || sol == nullptr || sat == nullptr || nv == nullptr)
    {
        return false;
    }
    double *pif = new double[n];
    int a = 0;
    int b = 0;
    // 计算双频组合观测值获取可用卫星数
    PIF(obs, n, sat, nav, nv, pif, &a, &b);
    //判断历元解算模式,判断用星数 nv
    SppSolveMode mode = SPP_MODE_GPS_BDS;
    if (a > 0 && b > 0 && *nv >= 5)
    {
        mode = SPP_MODE_GPS_BDS;
    }
    else if (a >= 4)
    {
        mode = SPP_MODE_GPS_ONLY;
        *nv = a;
    }
    else if (b >= 4)
    {
        mode = SPP_MODE_BDS_ONLY;
        *nv = b;
    }
    else
    {
        sol->stat = 0;
        delete[] pif;
        return false;
    }
    // 迭代解算
    Matrix *m = IterativeSolution(obs, sat, nav, n, *nv, mode, init_sol);
    if (m == nullptr)
    {
        sol->stat = 0;
        delete[] pif;
        return false;
    }
    Matrix X = m[0]; // 坐标和接收机钟差
    Matrix B = m[1];
    Matrix w = m[2];
    Matrix dx = m[3];  
    *nv = (int)B.size();// 最后一次迭代真实值与参考值的差
    Matrix v = sub(mul(B, dx), w); // 残差
    double vtv = mul(transpose(v), v)[0][0];
    double sigma0 = 0.0;
    int nx = SppUnknownCount(mode);
    if (*nv > nx)
    {
        sigma0 = sqrt(vtv / double(*nv - nx));//单位权中误差
    }
    Matrix Qxx = inverse(mul(transpose(B), B));//协因数阵
    Matrix Dxx = scalar_mul(sigma0*sigma0, Qxx);//方差阵
    double PDOP = sqrt(Qxx[0][0] + Qxx[1][1] + Qxx[2][2]);//几何精度因子
    sol->time = obs[0].time;
    sol->XYZ[0] = X[0][0];
    sol->XYZ[1] = X[1][0];
    sol->XYZ[2] = X[2][0];
    //GPS/BDS 接收机钟差清零
    sol->dtr[0] = 0.0;
    sol->dtr[1] = 0.0;
    if (mode == SPP_MODE_GPS_BDS)
    {
        sol->dtr[0] = X[3][0]; // GPS
        sol->dtr[1] = X[4][0]; // BDS
    }
    else if (mode == SPP_MODE_GPS_ONLY)
    {
        sol->dtr[0] = X[3][0]; // GPS
    }
    else
    {
        sol->dtr[1] = X[3][0]; // BDS
    }
    sol->Q = Qxx;
    sol->pdop = PDOP;
    sol->sigma0 = sigma0;
    sol->ns = *nv;
    sol->stat = 1;//标记定位成功
    delete[] pif;
    delete[] m;
    return true;
}
//#########################速度解算###################################
// 多普勒转距离率函数
double DopplerToRangeRate(obsd_t o)
{
    int prn = 0;
    int sys = satsys(o.sat, &prn);

    double f = 0.0;

    if (sys == SYS_GPS)
    {
        f = FREQ_GPS_L1;
    }
    else if (sys == SYS_CMP)
    {
        f = FREQ_BDS_B1;
    }
    else
    {
        return 0.0;
    }

    double lambda = Clight / f;
    // 卫星接近时 Doppler 为正，而距离率为负
    double D_mps = -lambda * o.D[0];

    return D_mps;
}
// 统计可参与测速的观测数量。
int CountSpeedObs(obsd_t *obs, satpos_t *sat, int n, sol_t *sol)
{
    if (obs == nullptr || sat == nullptr || sol == nullptr)
    {
        return 0;
    }

    if (sol->stat != 1)
    {
        return 0;
    }

    int nv = 0;

    for (int i = 0; i < n; i++)
    {
        obsd_t o = obs[i];
        satpos_t s = sat[i];

        int prn = 0;
        int sys = satsys(o.sat, &prn);

        if (sys != SYS_GPS && sys != SYS_CMP)
        {
            continue;
        }

        // 多普勒不能为0
        if (fabs(o.D[0]) < 1e-9)
        {
            continue;
        }

        // 卫星位置不能为0
        if (s.pos[0] == 0.0 && s.pos[1] == 0.0 && s.pos[2] == 0.0)
        {
            continue;
        }

        // 卫星速度不能为0
        if (s.vel[0] == 0.0 && s.vel[1] == 0.0 && s.vel[2] == 0.0)
        {
            continue;
        }

        nv++;
    }

    return nv;
}
// 速度最小二乘求解，未知数为 Vx、Vy、Vz 和接收机钟速
Matrix *SpeedLeastSquares(obsd_t *obs, satpos_t *sat, int n, sol_t *sol, int nv)
{
    Matrix B = zeros(nv, 4);
    Matrix w = zeros(nv, 1);

    int index = 0;

    for (int i = 0; i < n; i++)
    {
        obsd_t o = obs[i];
        satpos_t s = sat[i];

        int prn = 0;
        int sys = satsys(o.sat, &prn);

        if (sys != SYS_GPS && sys != SYS_CMP)
        {
            continue;
        }

        if (fabs(o.D[0]) < 1e-9)
        {
            continue;
        }

        if (s.pos[0] == 0.0 && s.pos[1] == 0.0 && s.pos[2] == 0.0)
        {
            continue;
        }

        if (s.vel[0] == 0.0 && s.vel[1] == 0.0 && s.vel[2] == 0.0)
        {
            continue;
        }

        double Xs = s.pos[0];
        double Ys = s.pos[1];
        double Zs = s.pos[2];

        double Xr = sol->XYZ[0];
        double Yr = sol->XYZ[1];
        double Zr = sol->XYZ[2];

        double delX = Xs - Xr;
        double delY = Ys - Yr;
        double delZ = Zs - Zr;

        double rho = sqrt(delX * delX + delY * delY + delZ * delZ);

        if (rho < 1.0)
        {
            continue;
        }

        // 方向余弦
        double l = delX / rho;
        double m = delY / rho;
        double nn = delZ / rho;

        // 多普勒转距离率，单位 m/s
        double D_mps = DopplerToRangeRate(o);

        // 卫星速度沿视线方向投影
        double sat_rate = l * s.vel[0] + m * s.vel[1] + nn * s.vel[2];
        if (index >= nv)
        {
            return nullptr;
        }
        B[index][0] = -l;
        B[index][1] = -m;
        B[index][2] = -nn;
        B[index][3] = 1.0;

        // w矩阵
        w[index][0] = D_mps - sat_rate + Clight * s.dclk;

        index++;
    }
    if (index != nv)
    {
        return nullptr;
    }
    Matrix N, x;

    N = mul(transpose(B), B);
    N = inverse(N);

    x = mul(N, transpose(B));
    x = mul(x, w);

    Matrix v = sub(mul(B, x), w);

    Matrix *result = new Matrix[4];

    result[0] = x; // Vx Vy Vz dtrd
    result[1] = B; 
    result[2] = w;
    result[3] = v; // 残差

    return result;
}
// 单点测速入口：定位成功后，用多普勒观测估计速度
bool SPP_Speed(obsd_t *obs, int n, sol_t *sol, satpos_t *sat, solvel_t *vsol, int *nv)
{
    if (obs == nullptr || sat == nullptr || sol == nullptr || vsol == nullptr || nv == nullptr)
    {
        return false;
    }

    if (sol->stat != 1)
    {
        vsol->stat = 0;
        *nv = 0;
        return false;
    }

    *nv = CountSpeedObs(obs, sat, n, sol);

    // 测速至少4颗卫星
    if (*nv < 4)
    {
        vsol->stat = 0;
        return false;
    }

    Matrix *m = nullptr;

    try
    {
        m = SpeedLeastSquares(obs, sat, n, sol, *nv);
    }
    catch (...)
    {
        vsol->stat = 0;
        if (m != nullptr)
        {
            delete[] m;
        }
        return false;
    }

    if (m == nullptr)
    {
        vsol->stat = 0;
        return false;
    }

    Matrix x = m[0];
    Matrix B = m[1];
    Matrix w = m[2];
    Matrix v = m[3];

    double vtv = mul(transpose(v), v)[0][0];

    double sigma0 = 0.0;//测速单位权中误差
    if (*nv > 4)
    {
        sigma0 = sqrt(vtv / double(*nv - 4));
    }

    Matrix Qxx = inverse(mul(transpose(B), B));//测速参数协因数阵

    vsol->time = sol->time;

    vsol->V[0] = x[0][0];
    vsol->V[1] = x[1][0];
    vsol->V[2] = x[2][0];

    // 接收机钟速，单位 m/s
    vsol->dtrd = x[3][0];

    vsol->Q = Qxx;
    vsol->sigma0 = sigma0;
    vsol->ns = *nv;
    vsol->stat = 1;

    delete[] m;

    return true;
}
