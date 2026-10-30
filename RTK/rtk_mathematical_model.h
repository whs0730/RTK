#pragma once
#include <string>
#include <vector>

#include "rtk_prepare.h"
#include "matrix.h"

enum class RtkObservationType
{
    Phase,
    Code
};

// 短基线非组合模型不估计电离层、对流层或系统间偏差，假定差分后残差可忽略。
// GPS和BDS分别选参考星，频点共享三维位置，每个星对/频点各估计一个双差模糊度，单位周。
struct rtk_model_options_t
{
    bool use_gps = true;
    bool use_bds = true;
    int frequency_count = 2;       // 1：L1/B1I，2：L1/L2和B1I/B3I
    bool exclude_bds_geo_reference = true; // 仅禁止GEO作参考星，不剔除GEO非参考星
    // 已知坐标优先于固定SPP基准；两者均关闭时使用本历元基站SPP。
    bool use_known_base_xyz = false;
    double known_base_xyz[3]{};    // ECEF XYZ，单位m；不覆盖epoch.base_spp
    bool use_fixed_spp_base_xyz = false;
    double fixed_spp_base_xyz[3]{}; // ECEF XYZ，单位m；固定的SPP概略值，不作为精确真值标记
};

// B/L的行映射；共视下标指向epoch.common，不是卫星编号。
struct dd_row_info_t
{
    int sys = SYS_NONE;
    int frequency = 0;             // 0：L1/B1I，1：L2/B3I
    int sat = 0;                   // 非参考星
    int reference_sat = 0;
    RtkObservationType type = RtkObservationType::Phase;

    int ambiguity_column = -1;     // 伪距行没有模糊度，保持-1
    int common_index = -1;         // epoch.common中的非参考星下标
    int reference_common_index = -1;
    double wavelength = 0.0;       // m/cycle
    unsigned char lli = 0;         // 两站、两颗卫星的LLI按位或；供后续半周固定检查
};

// B矩阵第3列起的模糊度列映射，列号从0开始。
struct dd_ambiguity_info_t
{
    int column = -1;
    int sys = SYS_NONE;
    int frequency = 0;
    int sat = 0;
    int reference_sat = 0;
    double wavelength = 0.0;
};

struct rtk_function_model_t
{
    gtime_t time{};

    int gps_reference_sat = 0;
    int bds_reference_sat = 0;
    int gps_satellite_count = 0;    // 实际参与模型的共视星数，包含参考星
    int bds_satellite_count = 0;
    int frequency_count = 0;
    int double_difference_count = 0; // 非参考星对数，不乘频点/观测类型
    double time_diff = 0.0;        // 基准站接收时刻减流动站接收时刻，s
    double base_xyz[3]{};          // 基准站坐标作为已知量，m
    bool base_xyz_is_known = false; // 区分已知基站坐标与SPP概略坐标
    double rover_xyz[3]{};         // 流动站线性化展开点，m

    Matrix B;                       // 设计矩阵
    Matrix L;                       // 观测减计算值，单位：m

    int observation_count = 0;
    int parameter_count = 0;

    vector<dd_row_info_t> rows;
    vector<dd_ambiguity_info_t> ambiguities;
    string error_message;      // 构模失败时的原因；失败模型的B/L为空
};

// 短基线单历元非组合模型：L = B * [dX,dY,dZ,N_DD...]^T + noise。
// dX/dY/dZ为坐标改正(m)，N_DD为本历元的绝对双差模糊度(cycle)。
// L为原始双差观测减双差几何距离，未扣除任何模糊度初值，单位统一为m。
// 单差取流动站减基准站，双差再取目标星减参考星；相位方程的模糊度系数为+lambda。
// 各站使用epoch.common中各自的卫星位置；输入epoch不被修改。
// 基站几何坐标由options选择；流动站展开点始终使用调用者传入的rover_xyz。
// 行序：GPS相位f0/f1、GPS伪距f0/f1，再BDS相位f0/f1、BDS伪距f0/f1。
// 每块按sat排序；默认GPS/BDS双系统双频，单频时省去f1。
bool BuildRtkFunctionModel(
    const rtk_epoch_t& epoch,
    const double rover_xyz[3],
    const rtk_config_t& config,
    rtk_function_model_t& model
);

// 可选择单/双系统和单/双频；在输入epoch已有共视集合中重新检查所选频点。
// 准备阶段未保存到common的卫星不会由此接口恢复。
bool BuildRtkFunctionModel(
    const rtk_epoch_t& epoch,
    const double rover_xyz[3],
    const rtk_config_t& config,
    const rtk_model_options_t& options,
    rtk_function_model_t& model
);

// 按课程参考格式输出参考星、共视星、坐标、设计矩阵和双差残差。
bool WriteRtkFunctionModel(ostream& out,const rtk_function_model_t& model);
