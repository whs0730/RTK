#pragma once

#include "rtk_estimation.h"

struct rtk_ambiguity_options_t
{
    double ratio_threshold = 2.0;   // 任务五课件：次优/最优残差二次型严格大于阈值才固定
    int max_search_iterations = 5000; // 搜索树循环次数上限；耗尽后失败，不接受暂存候选
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
    bool valid = false;             // 至少已有通过输入检查、可输出的浮点位置
    bool fixed = false;             // 整数搜索、Ratio及固定位置/协方差检查均通过
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
// 通过整数可逆变换降相关，在变换空间搜索，再回原空间计算两组候选的二次型。
bool SearchRtkIntegerAmbiguities(const Matrix& float_ambiguities, const Matrix& covariance,
    int max_search_iterations, rtk_integer_candidates_t& candidates);

// 仅更新当前历元输出，不修改浮点solution或跨历元filter，也不进行fix-and-hold。
// Ratio未通过属于正常FLOAT结果。搜索/更新失败返回false；若浮点位置本身有效，
// result在已复制有效浮点位置后仍保留FLOAT，原因见error_message；输入基础检查失败时valid=false。
// 存在未解决半周标志时不搜索，返回true且fixed=false；Ratio等于阈值也不固定。
bool ResolveRtkAmbiguities(const rtk_solution_t& float_solution,
    const rtk_ambiguity_options_t& options, rtk_ambiguity_solution_t& result);

// 详细文本含候选与列映射；CSV的dX/dY/dZ为基线，标准差为所选FIX或FLOAT位置的标准差。
// 未进行或未完成搜索时，CSV的Ratio及候选二次型留空，不能把空值当作Ratio=0。
bool WriteRtkAmbiguitySolution(ostream& out, const rtk_ambiguity_solution_t& solution);
bool WriteRtkAmbiguitySolutionCsvHeader(ostream& out);
bool WriteRtkAmbiguitySolutionCsv(ostream& out, const rtk_ambiguity_solution_t& solution);
