#pragma once

#include "rtk_estimation.h"

struct rtk_ambiguity_options_t
{
    double ratio_threshold = 2.0;   // 任务五课件：次优/最优残差二次型严格大于阈值才固定
    int max_search_iterations = 5000;
};

struct rtk_integer_candidates_t
{
    bool valid = false;
    Matrix integers;                // n*2，按残差二次型从小到大排列；单位cycle
    double squared_norms[2]{};      // (N_float-N_integer)' QNN^-1 (N_float-N_integer)
    int search_iterations = 0;
    string error_message;
};

struct rtk_ambiguity_solution_t
{
    bool valid = false;
    bool fixed = false;
    gtime_t time{};
    double ratio = 0.0;             // 最优二次型为0且次优大于0时取正无穷
    double ratio_threshold = 2.0;
    rtk_integer_candidates_t candidates;
    Matrix float_ambiguities;       // 与输入solution.ambiguities同序，cycle
    vector<dd_ambiguity_info_t> ambiguities;
    double rover_xyz[3]{};
    double baseline_xyz[3]{};
    double baseline_enu[3]{};
    Matrix position_covariance;     // 3*3，m²；FLOAT时保留原QXX
    string error_message;
};

// 联合搜索全部模糊度的两组整数候选，保留QNN中的跨系统、跨频率相关项。
// 输入为n*1的周数及n*n的cycle²协方差；失败时不返回未完成搜索的候选。
bool SearchRtkIntegerAmbiguities(const Matrix& float_ambiguities, const Matrix& covariance,
    int max_search_iterations, rtk_integer_candidates_t& candidates);

// 仅更新当前历元输出，不修改浮点solution或跨历元filter，也不进行fix-and-hold。
// Ratio未通过属于正常FLOAT结果。搜索/更新失败返回false；若浮点位置本身有效，
// result仍保留可输出的FLOAT位置，原因见error_message。无效浮点输入使result.valid=false。
bool ResolveRtkAmbiguities(const rtk_solution_t& float_solution,
    const rtk_ambiguity_options_t& options, rtk_ambiguity_solution_t& result);

bool WriteRtkAmbiguitySolution(ostream& out, const rtk_ambiguity_solution_t& solution);
bool WriteRtkAmbiguitySolutionCsvHeader(ostream& out);
bool WriteRtkAmbiguitySolutionCsv(ostream& out, const rtk_ambiguity_solution_t& solution);
