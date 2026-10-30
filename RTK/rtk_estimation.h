#pragma once

#include "rtk_stochastic_model.h"

enum class RtkEstimationMethod
{
    LeastSquares,
    ExtendedKalman
};

struct rtk_estimation_options_t
{
    RtkEstimationMethod method = RtkEstimationMethod::ExtendedKalman;
    bool kinematic = true;                 // 动态：位置每历元用SPP重置；静态：传递位置及协方差
    int max_iterations = 10;                // 单历元最小二乘的重新线性化次数
    double position_tolerance_m = 1e-4;      // LS三维坐标改正范数的收敛限，m
    double initial_position_sigma_m = 30.0; // EKF首次或动态位置重置时的先验标准差，m
    double initial_ambiguity_sigma_m = 30.0; // 除以波长后才是模糊度初始标准差(周)
    double position_process_sigma = 1e-4;   // 静态位置随机游走，m/sqrt(s)
    double ambiguity_process_sigma = 1e-4;  // 双差模糊度随机游走，cycle/sqrt(s)
    double reset_gap_s = 300.0;            // 距上一成功历元超过此间隔时重置整个滤波器，s
    double geometry_free_threshold_m = 0.05; // 站间单差GF在两次成功历元间的跳变限，m
    double melbourne_wubbena_threshold_cycles = 5.0; // 站间单差MW跳变限，周
};

// 周跳检查保留站间单差组合量，不依赖星间参考星，换参考星不会引入组合量跳变。
struct rtk_phase_history_t
{
    int sat = 0;
    int frequency_count = 0;
    unsigned char base_code[2]{};
    unsigned char rover_code[2]{};
    unsigned char base_half_cycle[2]{};    // 两站分别记录，避免按位或掩盖其中一站的状态变化
    unsigned char rover_half_cycle[2]{};
    double geometry_free_m = 0.0;
    double melbourne_wubbena_cycles = 0.0;
};

// 仅保存上一成功历元；模糊度列映射随状态一起传递，不能把前后历元的同列当成同一颗星。
struct rtk_filter_state_t
{
    bool initialized = false;
    gtime_t time{};
    double base_xyz[3]{};
    bool base_xyz_is_known = false;
    bool kinematic = true;
    Matrix x;                            // [绝对ECEF坐标(m), 绝对双差模糊度(周)]，列向量
    Matrix covariance;                   // 与x同序，包含位置和模糊度的交叉协方差
    vector<dd_ambiguity_info_t> ambiguities;
    vector<rtk_phase_history_t> phase_history;
};

// 当前历元的浮点结果。与filter.x不同，solution.x前三项按线性模型表示坐标改正。
struct rtk_solution_t
{
    bool valid = false;
    gtime_t time{};
    RtkEstimationMethod method = RtkEstimationMethod::LeastSquares;
    int observation_count = 0;
    int parameter_count = 0;
    int degrees_of_freedom = 0;           // nv-nx，仅用于有多余观测的LS后验中误差
    int iterations = 0;
    int initialized_ambiguities = 0;      // 含首历元、新星、周跳及无法变换的双差状态
    int carried_ambiguities = 0;
    bool filter_reinitialized = false;
    double base_xyz[3]{};
    double rover_xyz[3]{};
    double baseline_xyz[3]{};
    double baseline_enu[3]{};               // E,N,U，基站处的局部坐标系
    Matrix x;                             // [dX,dY,dZ,N...]；EKF前三项为0，更新后坐标见rover_xyz
    Matrix covariance;                    // 与x同序：m²、m*cycle、cycle²
    Matrix ambiguity_covariance;          // 直接截取完整协方差的NN块，供LAMBDA使用
    Matrix position_ambiguity_covariance;  // XN块，供固定后位置更新使用
    Matrix residuals;                     // B*x-L，m；EKF为更新后非线性模型残差
    vector<dd_row_info_t> rows;
    vector<dd_ambiguity_info_t> ambiguities;
    double phase_rms_m = 0.0;
    double code_rms_m = 0.0;
    double weighted_residual_sum = 0.0;   // v'Pv，使用双差观测完整权阵，含非对角相关项
    double sigma0 = 0.0;                  // 仅LS且nv>nx时有效：sqrt(v'Pv/(nv-nx))，无量纲
    double innovation_nis = 0.0;          // 仅EKF：innovation' * S^-1 * innovation
    string error_message;
};

// 输入B/L的行序必须与D/P完全一致，前三列为坐标改正，其余列由ambiguities定义。
// 只解当前线性模型；模糊度是绝对双差值，不累加初值。协方差沿用D的绝对噪声尺度。
bool SolveRtkLeastSquares(const rtk_function_model_t& function_model,
    const rtk_stochastic_model_t& stochastic_model, rtk_solution_t& solution);

// LS每次更新展开点后重建B/L；EKF在预报位置处线性化，并传递双差状态及完整协方差。
// EKF要求时刻严格递增；filter由调用方跨历元持有，LS分支不读写其中的历史状态。
// 输入历元、两站SPP均不修改；失败时solution清空，filter保留上一成功历元。
// EKF拒绝NIS/nv或更新后v'Pv/nv超过25的历元，不提交滤波状态。
// 动态EKF不含速度状态：每历元用SPP重置位置先验，仅连续传递可映射的模糊度。
// 两种估计方式均输出浮点解；模糊度固定在独立模块完成，不反馈到此处的filter。
bool EstimateRtkEpoch(const rtk_epoch_t& epoch, const rtk_config_t& config,
    const rtk_model_options_t& model_options, const rtk_stochastic_options_t& stochastic_options,
    const rtk_estimation_options_t& options, rtk_filter_state_t& filter, rtk_solution_t& solution);

// TXT按课程参考格式显示浮点参数、完整协方差和基线总结；固定模块直接使用内存中的QNN/QXN。
// CSV中的dX/dY/dZ是流动站减基站的基线分量，不是solution.x中的线性化改正。
bool WriteRtkSolution(ostream& out, const rtk_solution_t& solution);
bool WriteRtkSolutionCsvHeader(ostream& out);
bool WriteRtkSolutionCsv(ostream& out, const rtk_solution_t& solution);
