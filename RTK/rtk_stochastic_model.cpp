#include "rtk_stochastic_model.h"
#include <cmath>
#include <exception>
#include <iomanip>
using namespace std;

// 同一系统、同一频点、同一种观测属于一个块，不使用tuple或类型别名。
static bool SameObservationBlock(const dd_row_info_t &first, const dd_row_info_t &second)
{
    return first.sys == second.sys && first.frequency == second.frequency && first.type == second.type;
}

static bool Fail(rtk_stochastic_model_t &model, const string &message)
{
    // 失败时清空矩阵和映射，仅保留诊断信息，防止继续使用上一历元的结果。
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

// 一条观测只需确认系统、频点和类型；构模与导出共用这三个判断。
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
            // 防止角度单位错误、零高度角除零或无效高度角进入模型。
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

// 导出前检查维数和数值，拒绝失败结果、残缺矩阵及不匹配的行信息。
static bool ValidOutput(const rtk_stochastic_model_t &model)
{
    const int n = model.observation_count;
    // 先检查基本状态，再检查后面需要访问的数组长度。
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

// 按行输出一个矩阵，D和P分别调用此函数，避免使用pair及初始化列表遍历。
static void WriteStochasticMatrix(ostream &out, const Matrix &matrix, const char *name)
{
    int n = static_cast<int>(matrix.size());
    out << name << ' ' << n << ' ' << n << '\n';
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            out << (j == 0 ? "" : " ") << matrix[i][j];
        }
        out << '\n';
    }
}

bool BuildRtkStochasticModel(const rtk_epoch_t &epoch, const rtk_function_model_t &function_model,
                             rtk_stochastic_model_t &model)
{
    return BuildRtkStochasticModel(epoch, function_model, rtk_stochastic_options_t{}, model);
}

bool BuildRtkStochasticModel(const rtk_epoch_t &epoch, const rtk_function_model_t &function_model,
                             const rtk_stochastic_options_t &options, rtk_stochastic_model_t &model)
{
    // 1. 只检查任务三使用的定权配置、行映射和历元。
    // B、L的维数与数值、参数个数等由任务二负责，这里不再重复检查。
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

    // 2. 读取任务二行映射，计算每行的目标星和参考星单差方差，并按观测类型分块。
    model.target_sd_variance.resize(n);
    model.reference_sd_variance.resize(n);
    for (int i = 0; i < n; ++i)
    {
        const dd_row_info_t &row = function_model.rows[i];
        if (!ValidObservationInfo(row, function_model.frequency_count))
        {
            return Fail(model, "观测行的系统、频点或类型无效");
        }

        // 必须先确认下标有效，才能读取共视卫星数组。
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
        // 与前面的行逐一比较，检查重复观测以及同一块是否使用了同一参考星下标。
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

    // 3. 将每个块的协方差填入D的原行号，跨块元素保持零。
    model.D = zeros(n, n);
    // processed记录哪些行已完成；用普通循环收集同一块的原行号。
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
        // 一个块：D = diag(v_目标星) + v_参考星 * 1*1ᵀ。
        // 对角元素为目标星与参考星单差方差之和，非对角元素为参考星单差方差。
        for (int a = 0; a < block_size; ++a)
        {
            int i = indices[a];
            for (int b = 0; b < block_size; ++b)
            {
                int j = indices[b];
                model.D[i][j] = reference_variance + (i == j ? model.target_sd_variance[i] : 0.0);
            }
        }
    }

    // 4. 完整的D构建好后，直接调用matrix.cpp中的矩阵求逆函数，得到权阵P。
    // inverse()在矩阵奇异等情况下会抛出异常，这里转换为构模失败信息。
    try
    {
        model.P = inverse(model.D);
    }
    catch (const exception &error)
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

    // 全局矩阵始终按任务二原行号填入；即使行序被交错排列，也无需重新排B和L。
    model.time = function_model.time;
    model.observation_count = n;
    model.frequency_count = function_model.frequency_count;
    model.gps_reference_sat = function_model.gps_reference_sat;
    model.bds_reference_sat = function_model.bds_reference_sat;
    model.options = options;
    model.rows = function_model.rows;
    if (!ValidOutput(model))
    {
        return Fail(model, "随机模型数值溢出或矩阵无效");
    }
    return true;
}

bool WriteRtkStochasticModel(ostream &out, const rtk_stochastic_model_t &model)
{
    if (!out || !ValidOutput(model))
    {
        return false;
    }
    const ios::fmtflags flags = out.flags();
    const streamsize precision = out.precision();
    out << scientific << setprecision(15);
    out << "STOCHASTIC_MODEL week=" << model.time.week << " sow=" << model.time.sec << " nv=" << model.observation_count
        << " frequencies=" << model.frequency_count << '\n';
    out << "WEIGHTING method=" << (model.options.method == RtkWeightingMethod::EqualVariance ? "EQUAL" : "ELEVATION")
        << " phase_sigma_m=" << model.options.phase_sigma_m
        << " code_phase_variance_ratio=" << model.options.code_phase_variance_ratio
        << " elevation_a_m=" << model.options.elevation_a_m << " elevation_b_m=" << model.options.elevation_b_m << '\n';
    out << "UNITS D=m^2 P=m^-2 ASSUMPTION independent_undifferenced_errors\n";
    out << "REFERENCE GPS=" << (model.gps_reference_sat ? sat2id(model.gps_reference_sat) : "NONE")
        << " BDS=" << (model.bds_reference_sat ? sat2id(model.bds_reference_sat) : "NONE") << '\n';
    for (int i = 0; i < model.observation_count; ++i)
    {
        const dd_row_info_t &row = model.rows[i];
        out << "ROW " << i << ' ' << (row.sys == SYS_GPS ? "GPS" : "BDS") << " f=" << row.frequency << ' '
            << (row.type == RtkObservationType::Phase ? "PHASE" : "CODE") << " sat=" << sat2id(row.sat)
            << " ref=" << sat2id(row.reference_sat) << " target_sd_variance_m2=" << model.target_sd_variance[i]
            << " reference_sd_variance_m2=" << model.reference_sd_variance[i] << '\n';
    }
    WriteStochasticMatrix(out, model.D, "D");
    WriteStochasticMatrix(out, model.P, "P");
    out << "END_STOCHASTIC_MODEL\n";
    out.flags(flags);
    out.precision(precision);
    return out.good();
}
