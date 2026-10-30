#include "rtk_stochastic_model.h"
#include <cmath>
#include <stdexcept>
#include <iomanip>
using namespace std;

// 同一系统、同一频点、同一种观测共用参考星单差，因此构成一个相关块。
static bool SameObservationBlock(const dd_row_info_t &first, const dd_row_info_t &second)
{
    return first.sys == second.sys && first.frequency == second.frequency && first.type == second.type;
}

static bool Fail(rtk_stochastic_model_t &model, const string &message)
{
    // 清空失败模型，避免调用方沿用上一历元的矩阵和行映射。
    model = rtk_stochastic_model_t{};
    model.error_message = message;
    return false;
}

static bool Positive(double value)
{
    return isfinite(value) && value > 0.0;
}

static bool ValidOptions(const rtk_stochastic_options_t &options)
{
    if (!Positive(options.code_phase_variance_ratio))
    {
        return false;
    }
    if (options.method == RtkWeightingMethod::EqualVariance)
    {
        return Positive(options.phase_sigma_m);
    }
    if (options.method != RtkWeightingMethod::Elevation)
    {
        return false;
    }

    if (!isfinite(options.elevation_a_m) || options.elevation_a_m < 0.0)
    {
        return false;
    }
    if (!isfinite(options.elevation_b_m) || options.elevation_b_m < 0.0)
    {
        return false;
    }
    return options.elevation_a_m > 0.0 || options.elevation_b_m > 0.0;
}

static bool ValidObservationInfo(const dd_row_info_t &row, int frequency_count)
{
    if (row.sys != SYS_GPS && row.sys != SYS_CMP)
    {
        return false;
    }
    if (row.frequency < 0 || row.frequency >= frequency_count)
    {
        return false;
    }
    if (row.type != RtkObservationType::Phase && row.type != RtkObservationType::Code)
    {
        return false;
    }
    return true;
}

// 两站非差方差相加得到站间单差方差；伪距方差按配置比例放大。
static bool SingleDifferenceVariance(const common_sat_t &sat, RtkObservationType type,
                                     const rtk_stochastic_options_t &options, double &variance)
{
    if (options.method == RtkWeightingMethod::EqualVariance)
    {
        variance = 2.0 * options.phase_sigma_m * options.phase_sigma_m;
    }
    else
    {
        variance = 0.0;
        double elevations[2] = {sat.base_azel[1], sat.rover_azel[1]};
        for (int station = 0; station < 2; ++station)
        {
            double elevation = elevations[station];
            // 高度角使用弧度，零高度角会使定权公式除零。
            if (!isfinite(elevation) || elevation <= 0.0 || elevation > PI / 2.0)
            {
                return false;
            }
            const double sine = sin(elevation);
            const double a = options.elevation_a_m;
            const double b = options.elevation_b_m;
            variance += a * a + (b / sine) * (b / sine);
        }
    }
    if (type == RtkObservationType::Code)
    {
        variance *= options.code_phase_variance_ratio;
    }
    return Positive(variance);
}

static bool ValidOutput(const rtk_stochastic_model_t &model)
{
    const int n = model.observation_count;
    if (!model.error_message.empty())
    {
        return false;
    }
    if (n <= 0 || n > 4 * MAXSAT)
    {
        return false;
    }
    if (!ValidOptions(model.options) || !isfinite(model.time.sec))
    {
        return false;
    }
    if (model.frequency_count != 1 && model.frequency_count != 2)
    {
        return false;
    }

    if (static_cast<int>(model.D.size()) != n || static_cast<int>(model.P.size()) != n)
    {
        return false;
    }
    if (static_cast<int>(model.rows.size()) != n)
    {
        return false;
    }
    if (static_cast<int>(model.target_sd_variance.size()) != n ||
        static_cast<int>(model.reference_sd_variance.size()) != n)
    {
        return false;
    }
    for (int i = 0; i < n; ++i)
    {
        const dd_row_info_t &row = model.rows[i];
        if (static_cast<int>(model.D[i].size()) != n || static_cast<int>(model.P[i].size()) != n)
        {
            return false;
        }
        if (!Positive(model.target_sd_variance[i]) || !Positive(model.reference_sd_variance[i]))
        {
            return false;
        }
        if (!Positive(model.D[i][i]) || !Positive(model.P[i][i]))
        {
            return false;
        }
        if (!ValidObservationInfo(row, model.frequency_count))
        {
            return false;
        }
        for (int j = 0; j < n; ++j)
        {
            if (!isfinite(model.D[i][j]) || !isfinite(model.P[i][j]))
            {
                return false;
            }
        }
    }
    return true;
}

bool BuildRtkStochasticModel(const rtk_epoch_t &epoch, const rtk_function_model_t &function_model,
                             rtk_stochastic_model_t &model)
{
    return BuildRtkStochasticModel(epoch, function_model, rtk_stochastic_options_t{}, model);
}

bool BuildRtkStochasticModel(const rtk_epoch_t &epoch, const rtk_function_model_t &function_model,
                             const rtk_stochastic_options_t &options, rtk_stochastic_model_t &model)
{
    // B/L由函数模型负责；定权只使用其历元和行映射。
    model = rtk_stochastic_model_t{};
    if (!ValidOptions(options))
    {
        return Fail(model, "随机模型定权参数无效");
    }

    const int n = function_model.observation_count;
    if (!function_model.error_message.empty())
    {
        return Fail(model, "任务二构模失败，无法建立随机模型");
    }
    if (epoch.ncommon <= 0 || epoch.ncommon > MAXSAT)
    {
        return Fail(model, "共视卫星数量无效");
    }
    if (n <= 0 || n > 4 * MAXSAT || static_cast<int>(function_model.rows.size()) != n)
    {
        return Fail(model, "观测行数量与行映射不一致");
    }
    if (function_model.frequency_count != 1 && function_model.frequency_count != 2)
    {
        return Fail(model, "频点个数必须为1或2");
    }

    if (!isfinite(epoch.time.sec) || !isfinite(function_model.time.sec))
    {
        return Fail(model, "历元时间无效");
    }
    if (epoch.time.week != function_model.time.week || fabs(epoch.time.sec - function_model.time.sec) > 1e-9)
    {
        return Fail(model, "任务二模型与输入共视数据不属于同一历元");
    }

    // 方差数组沿用B/L的行序，不重新筛星或选择参考星。
    model.target_sd_variance.resize(n);
    model.reference_sd_variance.resize(n);
    for (int i = 0; i < n; ++i)
    {
        const dd_row_info_t &row = function_model.rows[i];
        if (!ValidObservationInfo(row, function_model.frequency_count))
        {
            return Fail(model, "观测行的系统、频点或类型无效");
        }

        if (row.common_index < 0 || row.common_index >= epoch.ncommon)
        {
            return Fail(model, "目标星共视下标越界");
        }
        if (row.reference_common_index < 0 || row.reference_common_index >= epoch.ncommon)
        {
            return Fail(model, "参考星共视下标越界");
        }
        if (row.sat == row.reference_sat)
        {
            return Fail(model, "目标星与参考星不能相同");
        }

        const common_sat_t &target = epoch.common[row.common_index];
        const common_sat_t &reference = epoch.common[row.reference_common_index];
        const int reference_sat =
            row.sys == SYS_GPS ? function_model.gps_reference_sat : function_model.bds_reference_sat;
        if (target.sat != row.sat || reference.sat != row.reference_sat)
        {
            return Fail(model, "观测行与共视卫星编号不一致");
        }
        if (row.reference_sat != reference_sat)
        {
            return Fail(model, "观测行与任务二选定的参考星不一致");
        }
        if (satsys(row.sat, nullptr) != row.sys || satsys(row.reference_sat, nullptr) != row.sys)
        {
            return Fail(model, "目标星与参考星必须属于当前系统");
        }
        for (int j = 0; j < i; ++j)
        {
            const dd_row_info_t &previous = function_model.rows[j];
            if (!SameObservationBlock(row, previous))
            {
                continue;
            }
            if (row.sat == previous.sat)
            {
                return Fail(model, "任务二模型包含重复双差观测行");
            }
            if (row.reference_common_index != previous.reference_common_index)
            {
                return Fail(model, "同一双差块的参考星共视下标不一致");
            }
        }
        if (!SingleDifferenceVariance(target, row.type, options, model.target_sd_variance[i]) ||
            !SingleDifferenceVariance(reference, row.type, options, model.reference_sd_variance[i]))
        {
            return Fail(model, "单差方差无效，请检查高度角、标准差及方差比例");
        }
    }

    // 保留原行号填入相关块，跨块元素保持零。
    model.D = zeros(n, n);
    vector<int> processed(n, 0);
    for (int first = 0; first < n; ++first)
    {
        if (processed[first] != 0)
        {
            continue;
        }
        vector<int> indices;
        for (int i = first; i < n; ++i)
        {
            if (SameObservationBlock(function_model.rows[first], function_model.rows[i]))
            {
                indices.push_back(i);
                processed[i] = 1;
            }
        }
        int block_size = static_cast<int>(indices.size());
        const double reference_variance = model.reference_sd_variance[indices[0]];
        // D块=diag(v_target)+v_reference*1*1'；共同减去参考星产生正协方差。
        for (int a = 0; a < block_size; ++a)
        {
            int i = indices[a];
            for (int b = 0; b < block_size; ++b)
            {
                int j = indices[b];
                model.D[i][j] = reference_variance + (i == j ? model.target_sd_variance[i] : 0.0);
            }
            if (!Positive(model.D[i][i]))
            {
                return Fail(model, "随机模型数值溢出或矩阵无效");
            }
        }
    }

    // 求逆必须保留非对角项；奇异矩阵按构模失败处理。
    try
    {
        model.P = inverse(model.D);
    }
    catch (const runtime_error &error)
    {
        return Fail(model, string("权阵求逆失败：") + error.what());
    }

    // 通用求逆可能带来极小的非对称舍入误差，两侧取平均保持权阵对称。
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            double value = 0.5 * model.P[i][j] + 0.5 * model.P[j][i];
            model.P[i][j] = value;
            model.P[j][i] = value;
        }
    }
    // 维数与行映射由本函数构造；求逆后的有限性和正方差仍需检查。
    for (int i = 0; i < n; ++i)
    {
        if (!Positive(model.P[i][i]))
        {
            return Fail(model, "随机模型数值溢出或矩阵无效");
        }
        for (double value : model.P[i])
        {
            if (!isfinite(value))
            {
                return Fail(model, "随机模型数值溢出或矩阵无效");
            }
        }
    }

    model.time = function_model.time;
    model.observation_count = n;
    model.frequency_count = function_model.frequency_count;
    model.gps_reference_sat = function_model.gps_reference_sat;
    model.bds_reference_sat = function_model.bds_reference_sat;
    model.options = options;
    model.rows = function_model.rows;
    return true;
}

bool WriteRtkStochasticModel(ostream& out, const rtk_stochastic_model_t& model)
{
    if (!out.good() || !ValidOutput(model)) return false;
    const auto flags = out.flags();
    const auto precision = out.precision();
    out << fixed << setprecision(3) << ">RTK: " << model.time.week << ' ' << model.time.sec << "\n\n";
    writeMatrix(out, model.D, "双差协方差阵：");
    writeMatrix(out, model.P, "权阵：");
    out.flags(flags);
    out.precision(precision);
    return out.good();
}
