#pragma once

#include "rtk_mathematical_model.h"
#include <ostream>

// 两站、不同卫星、频点和观测类型的非差误差假定独立。
// 双差之间仍然相关：同一块的双差共用了参考星的站间单差。
enum class RtkWeightingMethod
{
    EqualVariance,  // 各站各颗卫星的非差相位方差相同
    Elevation       // 非差相位方差 = a² + b² / sin²(E)
};

struct rtk_stochastic_options_t
{
    RtkWeightingMethod method = RtkWeightingMethod::EqualVariance;
    double phase_sigma_m = 0.005;             // 等方差模型：非差相位标准差，m
    double code_phase_variance_ratio = 10000.0; // 伪距方差/相位方差，不是标准差比
    double elevation_a_m = 0.004;             // 高度角模型系数，m
    double elevation_b_m = 0.003;             // 高度角模型系数，m；E使用各站高度角(弧度)
};

struct rtk_stochastic_model_t
{
    gtime_t time{};
    int observation_count = 0;
    int frequency_count = 0;
    int gps_reference_sat = 0;
    int bds_reference_sat = 0;
    rtk_stochastic_options_t options;

    Matrix D;  // 双差观测方差协方差阵，m²；供滤波量测更新使用
    Matrix P;  // D的逆，m⁻²；供加权最小二乘使用，不能逐元素取倒数

    // 与函数模型B/L同序；两组单差方差用于构建各相关块。
    vector<dd_row_info_t> rows;
    vector<double> target_sd_variance;
    vector<double> reference_sd_variance;
    string error_message;  // 失败时D、P及行映射清空，只保留错误原因
};

// 沿用函数模型的筛星结果、参考星和行序；相位按米定权，不再乘波长平方。
// 默认等方差，伪距/相位方差比为10000；矩阵维数为观测数。
// 调用前应完成函数模型构建；B/L及参数个数由函数模型负责校验。
// 每块D=diag(v_target)+v_reference*1*1'，v是两站非差方差之和；跨块置零。
// 等方差相位块的双差对角为4*sigma²、非对角为2*sigma²；伪距块再乘方差比。
bool BuildRtkStochasticModel(
    const rtk_epoch_t& epoch,
    const rtk_function_model_t& function_model,
    rtk_stochastic_model_t& model
);

bool BuildRtkStochasticModel(
    const rtk_epoch_t& epoch,
    const rtk_function_model_t& function_model,
    const rtk_stochastic_options_t& options,
    rtk_stochastic_model_t& model
);

// 按课程参考格式导出历元、双差协方差阵和权阵；返回写入是否成功。
bool WriteRtkStochasticModel(ostream& out, const rtk_stochastic_model_t& model);
