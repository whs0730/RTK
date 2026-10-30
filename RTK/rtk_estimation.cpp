#include "rtk_estimation.h"
#include "coordinate.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <stdexcept>
using namespace std;

// 失败历元不能留下可被误用的上一份结果；跨历元滤波状态由调用方另外持有。
static bool Fail(rtk_solution_t& solution, const string& message)
{
    solution = rtk_solution_t{};
    solution.error_message = message;
    return false;
}

static bool ValidMatrix(const Matrix& matrix, int rows, int columns)
{
    if (static_cast<int>(matrix.size()) != rows || rows <= 0 || columns <= 0)
    {
        return false;
    }
    for (const auto& row : matrix)
    {
        if (static_cast<int>(row.size()) != columns)
        {
            return false;
        }
        for (double value : row)
        {
            if (!isfinite(value))
            {
                return false;
            }
        }
    }
    return true;
}

// 调用方须保证方阵。只消除求逆、矩阵乘法带来的非对称舍入差，不修复负方差或秩亏。
static void Symmetrize(Matrix& matrix)
{
    for (int i = 0; i < static_cast<int>(matrix.size()); ++i)
    {
        for (int j = 0; j < i; ++j)
        {
            const double value = 0.5 * matrix[i][j] + 0.5 * matrix[j][i];
            matrix[i][j] = matrix[j][i] = value;
        }
    }
}

// 用对角尺度归一化后检查对称正定性，避免米、周混合的矩阵直接使用绝对主元阈值。
// 此接口也用于法矩阵和新息协方差；Cholesky只负责检查，实际求逆复用matrix.cpp。
static Matrix CovarianceInverse(const Matrix& matrix)
{
    const int n = static_cast<int>(matrix.size());
    if (!ValidMatrix(matrix, n, n))
    {
        throw runtime_error("Invalid covariance dimensions or numeric values");
    }
    vector<double> scale(n);
    Matrix normalized = zeros(n, n);
    Matrix lower = zeros(n, n);
    for (int i = 0; i < n; ++i)
    {
        if (matrix[i][i] <= 0.0)
        {
            throw runtime_error("Covariance diagonal must be positive");
        }
        scale[i] = sqrt(matrix[i][i]);
    }
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            normalized[i][j] = matrix[i][j] / scale[i] / scale[j];
            if (fabs(matrix[i][j] - matrix[j][i]) / scale[i] / scale[j] > 1e-9)
            {
                throw runtime_error("Covariance is not symmetric");
            }
        }
    }
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            double value = normalized[i][j];
            for (int k = 0; k < j; ++k)
            {
                value -= lower[i][k] * lower[j][k];
            }
            if (i == j)
            {
                if (!isfinite(value) || value <= 1e-12)
                {
                    throw runtime_error("Covariance is singular or not positive definite");
                }
                lower[i][j] = sqrt(value);
            }
            else
            {
                lower[i][j] = value / lower[j][j];
            }
        }
    }
    // A=diag(scale)*normalized*diag(scale)，求逆后需在两侧分别除回尺度。
    Matrix result = inverse(normalized);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            result[i][j] = result[i][j] / scale[i] / scale[j];
        }
    }
    Symmetrize(result);
    return result;
}

// 相同维数不足以证明模型配套：参考星、观测类型或模糊度列不同都会使定权错位。
static bool SameRow(const dd_row_info_t& a, const dd_row_info_t& b)
{
    return a.sys == b.sys && a.frequency == b.frequency && a.sat == b.sat &&
        a.reference_sat == b.reference_sat && a.type == b.type &&
        a.ambiguity_column == b.ambiguity_column && a.wavelength == b.wavelength;
}

static bool ValidModels(const rtk_function_model_t& model, const rtk_stochastic_model_t& stochastic)
{
    const int nv = model.observation_count;
    const int nx = model.parameter_count;
    if (!model.error_message.empty() || !stochastic.error_message.empty() || nx < 4 || nv < nx ||
        !ValidMatrix(model.B, nv, nx) || !ValidMatrix(model.L, nv, 1) ||
        static_cast<int>(stochastic.D.size()) != nv || static_cast<int>(stochastic.P.size()) != nv ||
        static_cast<int>(model.rows.size()) != nv || static_cast<int>(stochastic.rows.size()) != nv ||
        static_cast<int>(model.ambiguities.size()) != nx - 3 || stochastic.observation_count != nv ||
        model.time.week != stochastic.time.week || model.time.sec != stochastic.time.sec ||
        model.gps_reference_sat != stochastic.gps_reference_sat ||
        model.bds_reference_sat != stochastic.bds_reference_sat ||
        model.frequency_count != stochastic.frequency_count)
    {
        return false;
    }
    for (int i = 0; i < nv; ++i)
    {
        if (!SameRow(model.rows[i], stochastic.rows[i]))
        {
            return false;
        }
    }
    for (int i = 0; i < nx - 3; ++i)
    {
        if (model.ambiguities[i].column != i + 3 || model.ambiguities[i].wavelength <= 0.0 ||
            !isfinite(model.ambiguities[i].wavelength) || model.ambiguities[i].frequency < 0 ||
            model.ambiguities[i].frequency >= 2 || model.ambiguities[i].sat < 1 ||
            model.ambiguities[i].sat > MAXSAT || model.ambiguities[i].reference_sat < 1 ||
            model.ambiguities[i].reference_sat > MAXSAT)
        {
            return false;
        }
    }
    for (int k = 0; k < 3; ++k)
    {
        if (!isfinite(model.base_xyz[k]) || !isfinite(model.rover_xyz[k]))
        {
            return false;
        }
    }
    if (hypot(model.base_xyz[0], model.base_xyz[1], model.base_xyz[2]) <= 1.0)
    {
        return false;
    }
    return true;
}

// EKF调用时展开点已经更新，前三项改正由调用方置零，模糊度仍是绝对双差周数。
static void FillSolution(const rtk_function_model_t& model, const rtk_stochastic_model_t& stochastic,
    const Matrix& x, const Matrix& covariance, rtk_solution_t& solution)
{
    solution.time = model.time;
    solution.observation_count = model.observation_count;
    solution.parameter_count = model.parameter_count;
    solution.degrees_of_freedom = model.observation_count - model.parameter_count;
    solution.x = x;
    solution.covariance = covariance;
    solution.rows = model.rows;
    solution.ambiguities = model.ambiguities;
    solution.residuals = sub(mul(model.B, x), model.L);
    const int na = model.parameter_count - 3;
    // 保留NN块内的相关性和XN交叉项，供后续整数固定及固定后位置更新使用。
    solution.ambiguity_covariance = zeros(na, na);
    solution.position_ambiguity_covariance = zeros(3, na);
    for (int i = 0; i < na; ++i)
    {
        for (int j = 0; j < na; ++j)
        {
            solution.ambiguity_covariance[i][j] = covariance[i + 3][j + 3];
        }
        for (int j = 0; j < 3; ++j)
        {
            solution.position_ambiguity_covariance[j][i] = covariance[j][i + 3];
        }
    }
    for (int k = 0; k < 3; ++k)
    {
        solution.base_xyz[k] = model.base_xyz[k];
        solution.rover_xyz[k] = model.rover_xyz[k] + x[k][0];
        solution.baseline_xyz[k] = solution.rover_xyz[k] - solution.base_xyz[k];
    }
    // 基线方向为流动站减基站；ENU的旋转基准也取基站，而非流动站。
    XYZ base{model.base_xyz[0], model.base_xyz[1], model.base_xyz[2]};
    BLH* blh = XYZtoBLH(base, 6378137.0, 1.0 / 298.257223563);
    ENU* enu = XYZtoENU(blh, solution.baseline_xyz);
    solution.baseline_enu[0] = enu->E;
    solution.baseline_enu[1] = enu->N;
    solution.baseline_enu[2] = enu->U;
    delete enu;
    delete blh;
    double phase_sum = 0.0, code_sum = 0.0;
    int phase_count = 0, code_count = 0;
    for (int i = 0; i < model.observation_count; ++i)
    {
        const double v = solution.residuals[i][0];
        if (model.rows[i].type == RtkObservationType::Phase)
        {
            phase_sum += v * v;
            ++phase_count;
        }
        else
        {
            code_sum += v * v;
            ++code_count;
        }
    }
    // 分类型RMS只描述观测拟合，不是坐标真误差；v'Pv另用完整权阵计入行间相关。
    solution.phase_rms_m = phase_count ? sqrt(phase_sum / phase_count) : 0.0;
    solution.code_rms_m = code_count ? sqrt(code_sum / code_count) : 0.0;
    solution.weighted_residual_sum = mul(mul(transpose(solution.residuals), stochastic.P), solution.residuals)[0][0];
    if (!isfinite(solution.weighted_residual_sum) || solution.weighted_residual_sum < -1e-8 ||
        !ValidMatrix(x, model.parameter_count, 1))
    {
        throw runtime_error("Nonfinite parameter estimate or residual statistic");
    }
    // D已给出绝对观测方差；sigma0只作为拟合统计量，不再用于缩放输出协方差。
    if (solution.method == RtkEstimationMethod::LeastSquares && solution.degrees_of_freedom > 0)
    {
        solution.sigma0 = sqrt(max(0.0, solution.weighted_residual_sum) / solution.degrees_of_freedom);
    }
    solution.valid = true;
}

bool SolveRtkLeastSquares(const rtk_function_model_t& model,
    const rtk_stochastic_model_t& stochastic, rtk_solution_t& solution)
{
    solution = rtk_solution_t{};
    if (!ValidModels(model, stochastic))
    {
        return Fail(solution, "Function/stochastic model dimensions, epochs or row mappings disagree");
    }
    try
    {
        // CovarianceInverse统一检查方阵结构、有限性和正定性；再核对P=D^-1，避免混用定权参数。
        const Matrix inverse_d = CovarianceInverse(stochastic.D);
        CovarianceInverse(stochastic.P);
        for (int i = 0; i < model.observation_count; ++i)
        {
            for (int j = 0; j < model.observation_count; ++j)
            {
                const double scale = sqrt(stochastic.P[i][i]) * sqrt(stochastic.P[j][j]);
                if (fabs(inverse_d[i][j] - stochastic.P[i][j]) > 1e-8 * scale)
                {
                    return Fail(solution, "Stochastic model P is not the inverse of D");
                }
            }
        }
        // L中未扣除模糊度初值，因此法方程直接解出[dX,dY,dZ,N...]，N无需再次累加。
        // x=(B'PB)^-1 B'PL，输出协方差为(B'PB)^-1，不乘后验sigma0²。
        const Matrix bt_p = mul(transpose(model.B), stochastic.P);
        const Matrix normal = mul(bt_p, model.B);
        const Matrix covariance = CovarianceInverse(normal);
        const Matrix x = mul(covariance, mul(bt_p, model.L));
        FillSolution(model, stochastic, x, covariance, solution);
        solution.iterations = 1;
        return true;
    }
    catch (const runtime_error& error)
    {
        return Fail(solution, string("RTK least squares failed: ") + error.what());
    }
}

static bool ValidOptions(const rtk_estimation_options_t& options)
{
    return (options.method == RtkEstimationMethod::LeastSquares || options.method == RtkEstimationMethod::ExtendedKalman) &&
        options.max_iterations > 0 && options.max_iterations <= 100 &&
        isfinite(options.position_tolerance_m) && options.position_tolerance_m > 0.0 &&
        isfinite(options.initial_position_sigma_m) && options.initial_position_sigma_m > 0.0 &&
        isfinite(options.initial_ambiguity_sigma_m) && options.initial_ambiguity_sigma_m > 0.0 &&
        isfinite(options.position_process_sigma) && options.position_process_sigma >= 0.0 &&
        isfinite(options.ambiguity_process_sigma) && options.ambiguity_process_sigma >= 0.0 &&
        isfinite(options.reset_gap_s) && options.reset_gap_s > 0.0 &&
        isfinite(options.geometry_free_threshold_m) && options.geometry_free_threshold_m > 0.0 &&
        isfinite(options.melbourne_wubbena_threshold_cycles) && options.melbourne_wubbena_threshold_cycles > 0.0;
}

static bool BuildModels(const rtk_epoch_t& epoch, const double xyz[3], const rtk_config_t& config,
    const rtk_model_options_t& model_options, const rtk_stochastic_options_t& stochastic_options,
    rtk_function_model_t& model, rtk_stochastic_model_t& stochastic, rtk_solution_t& solution)
{
    // 每次线性化均让随机模型沿用同一次函数模型的行映射，不能另行排序或换参考星。
    if (!BuildRtkFunctionModel(epoch, xyz, config, model_options, model))
    {
        return Fail(solution, model.error_message);
    }
    if (!BuildRtkStochasticModel(epoch, model, stochastic_options, stochastic))
    {
        return Fail(solution, stochastic.error_message);
    }
    return true;
}

static int FindHistory(const vector<rtk_phase_history_t>& history, int sat)
{
    for (int i = 0; i < static_cast<int>(history.size()); ++i)
    {
        if (history[i].sat == sat)
        {
            return i;
        }
    }
    return -1;
}

static vector<rtk_phase_history_t> PhaseHistory(const rtk_epoch_t& epoch, int frequencies)
{
    vector<rtk_phase_history_t> history;
    for (int i = 0; i < epoch.ncommon; ++i)
    {
        const common_sat_t& item = epoch.common[i];
        const obsd_t& base = epoch.base_obs.data[item.base_index];
        const obsd_t& rover = epoch.rover_obs.data[item.rover_index];
        const int sys = satsys(item.sat, nullptr);
        const double f0 = sys == SYS_GPS ? FREQ_GPS_L1 : FREQ_BDS_B1;
        const double f1 = sys == SYS_GPS ? FREQ_GPS_L2 : FREQ_BDS_B3;
        rtk_phase_history_t track;
        track.sat = item.sat;
        track.frequency_count = frequencies;
        for (int f = 0; f < frequencies; ++f)
        {
            track.base_code[f] = base.code[f];
            track.rover_code[f] = rover.code[f];
            track.base_half_cycle[f] = base.LLI[f] & (LLI_HALFC | LLI_HALFA);
            track.rover_half_cycle[f] = rover.LLI[f] & (LLI_HALFC | LLI_HALFA);
        }
        if (frequencies == 2)
        {
            // 先按站作差，再构成几何无关相位组合(GF)和MW组合；两者都不依赖星间参考星。
            const double l0 = Clight / f0 * (rover.L[0] - base.L[0]);
            const double l1 = Clight / f1 * (rover.L[1] - base.L[1]);
            const double p0 = rover.P[0] - base.P[0];
            const double p1 = rover.P[1] - base.P[1];
            track.geometry_free_m = l0 - l1;
            // 宽巷相位减窄巷伪距，再除以宽巷波长，MW最终以周为单位。
            track.melbourne_wubbena_cycles = ((f0 * l0 - f1 * l1) / (f0 - f1) -
                (f0 * p0 + f1 * p1) / (f0 + f1)) / (Clight / (f0 - f1));
        }
        history.push_back(track);
    }
    return history;
}

// slip也包括信号码和半周状态变化，标记后不能继续沿用旧模糊度。
// GF/MW仅在双频时可用，比较对象为上一成功滤波历元；单频依赖LLI及跟踪状态。
static vector<vector<bool>> DetectSlips(const rtk_epoch_t& epoch,
    const vector<rtk_phase_history_t>& current, const rtk_filter_state_t& filter,
    const rtk_estimation_options_t& options)
{
    vector<vector<bool>> slip(MAXSAT + 1, vector<bool>(2, false));
    // 按卫星/频点标记，目标星与参考星分别查询；参考星失锁会影响涉及它的全部双差。
    for (int i = 0; i < epoch.ncommon; ++i)
    {
        const common_sat_t& item = epoch.common[i];
        const obsd_t& b = epoch.base_obs.data[item.base_index];
        const obsd_t& r = epoch.rover_obs.data[item.rover_index];
        const auto& now = current[i];
        const int previous = FindHistory(filter.phase_history, item.sat);
        for (int f = 0; f < now.frequency_count; ++f)
        {
            slip[item.sat][f] = ((b.LLI[f] | r.LLI[f]) & (LLI_SLIP | LLI_HALFC)) != 0;
            if (previous >= 0)
            {
                const auto& old = filter.phase_history[previous];
                slip[item.sat][f] = slip[item.sat][f] || f >= old.frequency_count ||
                    now.base_code[f] != old.base_code[f] || now.rover_code[f] != old.rover_code[f] ||
                    now.base_half_cycle[f] != old.base_half_cycle[f] ||
                    now.rover_half_cycle[f] != old.rover_half_cycle[f];
            }
        }
        if (previous >= 0 && now.frequency_count == 2 && filter.phase_history[previous].frequency_count == 2)
        {
            const auto& old = filter.phase_history[previous];
            if (fabs(now.geometry_free_m - old.geometry_free_m) > options.geometry_free_threshold_m ||
                fabs(now.melbourne_wubbena_cycles - old.melbourne_wubbena_cycles) > options.melbourne_wubbena_threshold_cycles)
            {
                // 组合跳变无法确定是哪一频点周跳，两个频点都重新初始化，不修补原始载波。
                slip[item.sat][0] = slip[item.sat][1] = true;
            }
        }
    }
    return slip;
}

// 旧状态只存N(s,r_old)，其中N(r_old,r_old)=0。
// 新双差N(s,r_new)=N(s,r_old)-N(r_new,r_old)，两项分别以+1、-1写入变换行。
// sat等于旧参考星时该项为0，仍算可映射；旧状态找不到所需星时返回false并改走初始化。
static bool AddOldAmbiguity(const rtk_filter_state_t& filter, const dd_ambiguity_info_t& current,
    int sat, double sign, vector<double>& transform)
{
    for (const auto& old : filter.ambiguities)
    {
        if (old.sys == current.sys && old.frequency == current.frequency && old.wavelength == current.wavelength)
        {
            if (sat == old.reference_sat)
            {
                return true;
            }
            if (sat == old.sat)
            {
                transform[old.column] += sign;
                return true;
            }
        }
    }
    return false;
}

static double InitialAmbiguity(const rtk_function_model_t& model, const dd_ambiguity_info_t& ambiguity)
{
    int phase_row = -1, code_row = -1;
    for (int i = 0; i < model.observation_count; ++i)
    {
        const auto& row = model.rows[i];
        if (row.sys == ambiguity.sys && row.frequency == ambiguity.frequency &&
            row.sat == ambiguity.sat && row.reference_sat == ambiguity.reference_sat)
        {
            if (row.type == RtkObservationType::Phase)
            {
                phase_row = i;
            }
            else
            {
                code_row = i;
            }
        }
    }
    if (phase_row < 0 || code_row < 0)
    {
        throw logic_error("Missing code/phase pair for ambiguity initialization");
    }
    // 同一双差的相位和伪距含相同几何项，相减后除波长得到绝对模糊度初值(周)。
    return (model.L[phase_row][0] - model.L[code_row][0]) / ambiguity.wavelength;
}

// 状态顺序随当前模型变化，T把旧状态映射到当前卫星/频点/参考星组合。
// 不可传递的行全零：其初值、方差随后独立设置，与其他状态的先验交叉项保持零。
static void Predict(const rtk_epoch_t& epoch, const rtk_function_model_t& model,
    const rtk_estimation_options_t& options, const rtk_filter_state_t& filter,
    const vector<vector<bool>>& slips, bool restart, double dt, Matrix& prior, Matrix& covariance,
    rtk_solution_t& solution)
{
    const int nx = model.parameter_count;
    prior = zeros(nx, 1);
    covariance = zeros(nx, nx);
    vector<bool> carried(nx, false);
    if (!restart)
    {
        Matrix transform = zeros(nx, static_cast<int>(filter.x.size()));
        if (!options.kinematic)
        {
            for (int k = 0; k < 3; ++k)
            {
                transform[k][k] = 1.0;
                carried[k] = true;
            }
        }
        for (const auto& ambiguity : model.ambiguities)
        {
            const int c = ambiguity.column;
            if (slips[ambiguity.sat][ambiguity.frequency] || slips[ambiguity.reference_sat][ambiguity.frequency])
            {
                continue;
            }
            vector<double> row(filter.x.size(), 0.0);
            if (AddOldAmbiguity(filter, ambiguity, ambiguity.sat, 1.0, row) &&
                AddOldAmbiguity(filter, ambiguity, ambiguity.reference_sat, -1.0, row))
            {
                transform[c] = row;
                carried[c] = true;
            }
        }
        prior = mul(transform, filter.x);
        // T同时完成参考星变换和消失星裁剪，协方差用T*C*T'一起变换NN、XN、XX。
        covariance = mul(mul(transform, filter.covariance), transpose(transform));
    }
    for (int k = 0; k < 3; ++k)
    {
        if (!carried[k])
        {
            // 动态位置每历元独立初始化；对应变换行全零，旧的XN相关项也随之清除。
            prior[k][0] = epoch.rover_spp.XYZ[k];
            covariance[k][k] = options.initial_position_sigma_m * options.initial_position_sigma_m;
        }
        else
        {
            covariance[k][k] += options.position_process_sigma * options.position_process_sigma * dt;
        }
    }
    for (const auto& ambiguity : model.ambiguities)
    {
        const int c = ambiguity.column;
        if (carried[c])
        {
            // 配置值是随机游走标准差率，方差增长为sigma^2*dt，不是sigma*dt。
            covariance[c][c] += options.ambiguity_process_sigma * options.ambiguity_process_sigma * dt;
            ++solution.carried_ambiguities;
        }
        else
        {
            // 新星、周跳或无法用旧参考组成的双差均走初始化，不沿用同列的旧值。
            prior[c][0] = InitialAmbiguity(model, ambiguity);
            const double sigma = options.initial_ambiguity_sigma_m / ambiguity.wavelength;
            covariance[c][c] = sigma * sigma;
            ++solution.initialized_ambiguities;
        }
    }
    Symmetrize(covariance);
}

bool EstimateRtkEpoch(const rtk_epoch_t& epoch, const rtk_config_t& config,
    const rtk_model_options_t& model_options, const rtk_stochastic_options_t& stochastic_options,
    const rtk_estimation_options_t& options, rtk_filter_state_t& filter, rtk_solution_t& solution)
{
    solution = rtk_solution_t{};
    if (!ValidOptions(options))
    {
        return Fail(solution, "Invalid RTK estimation options");
    }
    try
    {
        double xyz[3];
        copy_n(epoch.rover_spp.XYZ, 3, xyz);
        rtk_function_model_t model;
        rtk_stochastic_model_t stochastic;
        if (!BuildModels(epoch, xyz, config, model_options, stochastic_options, model, stochastic, solution))
        {
            return false;
        }
        if (options.method == RtkEstimationMethod::LeastSquares)
        {
            // 只将位置估值带入下一次线性化；每次求出的N都是绝对双差值。
            for (int iteration = 1; iteration <= options.max_iterations; ++iteration)
            {
                if (!SolveRtkLeastSquares(model, stochastic, solution))
                {
                    return false;
                }
                solution.iterations = iteration;
                if (hypot(solution.x[0][0], solution.x[1][0], solution.x[2][0]) <= options.position_tolerance_m)
                {
                    return true;
                }
                copy_n(solution.rover_xyz, 3, xyz);
                if (!BuildModels(epoch, xyz, config, model_options, stochastic_options, model, stochastic, solution))
                {
                    return false;
                }
            }
            return Fail(solution, "RTK least squares did not converge");
        }

        const double dt = filter.initialized ? timediff(epoch.time, filter.time) : 0.0;
        if (filter.initialized && (!isfinite(dt) || dt <= 0.0))
        {
            return Fail(solution, "EKF requires strictly increasing epoch times");
        }
        // 间隔从上一成功更新算起；初始化模式或基站基准变化时，整组状态重新建立。
        bool restart = !filter.initialized || dt > options.reset_gap_s || filter.kinematic != options.kinematic ||
            filter.base_xyz_is_known != model.base_xyz_is_known;
        // 基站几何基准变化后旧模糊度不能直接沿用；使用逐历元基站SPP时也按此规则重置。
        for (int k = 0; k < 3 && filter.initialized; ++k)
        {
            restart = restart || fabs(filter.base_xyz[k] - model.base_xyz[k]) > 1e-4;
        }
        if (!restart && (!ValidMatrix(filter.x, static_cast<int>(filter.ambiguities.size()) + 3, 1) ||
            !ValidMatrix(filter.covariance, static_cast<int>(filter.x.size()), static_cast<int>(filter.x.size()))))
        {
            return Fail(solution, "Invalid saved RTK filter state");
        }
        if (!restart)
        {
            for (int i = 0; i < static_cast<int>(filter.ambiguities.size()); ++i)
            {
                if (filter.ambiguities[i].column != i + 3)
                {
                    return Fail(solution, "Invalid saved ambiguity column mapping");
                }
            }
        }
        const auto history = PhaseHistory(epoch, model.frequency_count);
        const auto slips = DetectSlips(epoch, history, filter, options);
        Matrix prior, prior_covariance;
        solution.method = RtkEstimationMethod::ExtendedKalman;
        solution.filter_reinitialized = restart;
        Predict(epoch, model, options, filter, slips, restart, dt, prior, prior_covariance, solution);
        for (int k = 0; k < 3; ++k)
        {
            xyz[k] = prior[k][0];
        }
        if (!BuildModels(epoch, xyz, config, model_options, stochastic_options, model, stochastic, solution))
        {
            return false;
        }
        // 参考星与筛选沿用epoch保存的高度角，改变展开点不会改变状态列映射。
        // L已在prior的位置处扣除几何量，B前三列作用于位置改正，不能再乘绝对ECEF坐标。
        Matrix linear_prior = prior;
        for (int k = 0; k < 3; ++k)
        {
            linear_prior[k][0] = 0.0;
        }
        // linear_prior只在模糊度分量保留预报值；新息还需扣掉lambda*N的预报贡献。
        const Matrix innovation = sub(model.L, mul(model.B, linear_prior));
        const Matrix c_ht = mul(prior_covariance, transpose(model.B));
        // S=B*C_prior*B'+D；D保留共享参考星产生的观测相关性。
        const Matrix s = add(mul(model.B, c_ht), stochastic.D);
        const Matrix s_inverse = CovarianceInverse(s);
        const Matrix gain = mul(c_ht, s_inverse);
        const Matrix correction = mul(gain, innovation);
        Matrix posterior = add(prior, correction);
        // Joseph形式保留两项协方差之和，减小直接计算(I-KB)*C_prior时的舍入影响。
        const Matrix a = sub(eye(model.parameter_count), mul(gain, model.B));
        Matrix covariance = add(mul(mul(a, prior_covariance), transpose(a)),
            mul(mul(gain, stochastic.D), transpose(gain)));
        Symmetrize(covariance);
        CovarianceInverse(covariance);
        // 用完整S保留双差相关性。25倍观测数对应整历元约5倍模型标准差，作为宽松的异常拒绝门限。
        solution.innovation_nis = mul(mul(transpose(innovation), s_inverse), innovation)[0][0];
        if (!isfinite(solution.innovation_nis) || solution.innovation_nis < -1e-8)
        {
            return Fail(solution, "Invalid EKF innovation statistic");
        }
        if (solution.innovation_nis > 25.0 * model.observation_count)
        {
            return Fail(solution, "RTK innovation rejected: NIS/observations = " +
                to_string(solution.innovation_nis / model.observation_count) + " exceeds 25");
        }
        // 在更新后坐标处重算几何量，仅检查后验残差，不再用同一批观测做第二次量测更新。
        for (int k = 0; k < 3; ++k)
        {
            xyz[k] = posterior[k][0];
        }
        if (!BuildModels(epoch, xyz, config, model_options, stochastic_options, model, stochastic, solution))
        {
            return false;
        }
        // FillSolution按“展开点+改正”输出坐标；展开点已是posterior，改正必须置零。
        Matrix final_x = posterior;
        for (int k = 0; k < 3; ++k)
        {
            final_x[k][0] = 0.0;
        }
        FillSolution(model, stochastic, final_x, covariance, solution);
        // 一次线性化可能在远离真值的SPP处通过先验检验，仍需检查更新后的非线性残差。
        if (solution.weighted_residual_sum > 25.0 * model.observation_count)
        {
            return Fail(solution, "RTK posterior residual rejected: vPv/observations = " +
                to_string(solution.weighted_residual_sum / model.observation_count) + " exceeds 25");
        }
        solution.iterations = 1;
        rtk_filter_state_t next;
        next.initialized = true;
        next.time = epoch.time;
        next.x = posterior;
        next.covariance = covariance;
        next.ambiguities = model.ambiguities;
        next.phase_history = history;
        next.kinematic = options.kinematic;
        next.base_xyz_is_known = model.base_xyz_is_known;
        copy_n(model.base_xyz, 3, next.base_xyz);
        // 通过全部检查后才提交，失败历元保留上一成功滤波状态。
        filter = next;
        return true;
    }
    catch (const runtime_error& error)
    {
        return Fail(solution, string("RTK estimation failed: ") + error.what());
    }
}

static const char* MethodName(RtkEstimationMethod method)
{
    return method == RtkEstimationMethod::LeastSquares ? "LS" : "EKF";
}

bool WriteRtkSolution(ostream& out, const rtk_solution_t& solution)
{
    if (!out.good() || !solution.valid || !ValidMatrix(solution.x, solution.parameter_count, 1) ||
        !ValidMatrix(solution.covariance, solution.parameter_count, solution.parameter_count)) return false;
    const auto flags = out.flags();
    const auto precision = out.precision();
    out << fixed << setprecision(3) << ">RTK: " << solution.time.week << ' ' << solution.time.sec << "\n\n";
    out << "浮点解信息：\n流动站坐标（XYZ）：";
    for (double value : solution.rover_xyz) out << setw(15) << value;
    out << "\n\n基线+模糊度参数：浮点\ndxyz:";
    for (double value : solution.baseline_xyz) out << setw(11) << value;
    out << '\n';
    for (const auto& item : solution.ambiguities)
        out << sat2id(item.sat) << " f" << item.frequency + 1 << " DDAmb:"
            << setw(11) << solution.x[item.column][0] << '\n';
    out << '\n';
    writeMatrix(out, solution.covariance, "浮点解方差：");
    out << "浮点解总结：\n"
        << "Week  Sow(sec)      dX(m)      dY(m)      dZ(m)      dN(m)      dE(m)      dU(m)"
        << "  dX_std(m) dY_std(m) dZ_std(m) obsL_RMS(m) obsP_RMS(m)\n";
    out << setw(4) << solution.time.week << setw(12) << solution.time.sec << setprecision(4);
    for (double value : solution.baseline_xyz) out << setw(11) << value;
    out << setw(11) << solution.baseline_enu[1] << setw(11) << solution.baseline_enu[0] << setw(11) << solution.baseline_enu[2];
    for (int i = 0; i < 3; ++i) out << setw(11) << sqrt(solution.covariance[i][i]);
    out << setw(12) << solution.phase_rms_m << setw(12) << solution.code_rms_m << "\n\n";
    out.flags(flags);
    out.precision(precision);
    return out.good();
}

bool WriteRtkSolutionCsvHeader(ostream& out)
{
    out << "week,sow,status,method,nv,nx,iterations,filter_reset,initialized_N,carried_N,"
        "X_m,Y_m,Z_m,dX_m,dY_m,dZ_m,E_m,N_m,U_m,std_X_m,std_Y_m,std_Z_m,"
        "phase_rms_m,code_rms_m,vPv,ls_sigma0,innovation_nis\n";
    return out.good();
}

bool WriteRtkSolutionCsv(ostream& out, const rtk_solution_t& solution)
{
    if (!out.good() || !solution.valid)
    {
        return false;
    }
    const auto flags = out.flags();
    const auto precision = out.precision();
    out << defaultfloat << setprecision(15);
    out << solution.time.week << ',' << solution.time.sec << ",FLOAT," << MethodName(solution.method)
        << ',' << solution.observation_count << ',' << solution.parameter_count << ',' << solution.iterations
        << ',' << solution.filter_reinitialized << ',' << solution.initialized_ambiguities << ',' << solution.carried_ambiguities;
    for (double value : solution.rover_xyz) { out << ',' << value; }
    for (double value : solution.baseline_xyz) { out << ',' << value; }
    for (double value : solution.baseline_enu) { out << ',' << value; }
    for (int k = 0; k < 3; ++k) { out << ',' << sqrt(solution.covariance[k][k]); }
    out << ',' << solution.phase_rms_m << ',' << solution.code_rms_m << ',' << solution.weighted_residual_sum << ',';
    // LS中误差和EKF新息统计各填自己的列，不适用的列留空，避免把0误当成有效统计量。
    if (solution.method == RtkEstimationMethod::LeastSquares && solution.degrees_of_freedom > 0) { out << solution.sigma0; }
    out << ',';
    if (solution.method == RtkEstimationMethod::ExtendedKalman) { out << solution.innovation_nis; }
    out << '\n';
    out.flags(flags);
    out.precision(precision);
    return out.good();
}
