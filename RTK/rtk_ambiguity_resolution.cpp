#include "rtk_ambiguity_resolution.h"
#include "coordinate.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <set>
#include <stdexcept>
#include <tuple>


bool ValidMatrix(const Matrix& a, int rows, int columns)
{
    if (rows <= 0 || columns <= 0 || static_cast<int>(a.size()) != rows) return false;
    for (const auto& row : a)
    {
        if (static_cast<int>(row.size()) != columns) return false;
        for (double v : row) if (!isfinite(v)) return false;
    }
    return true;
}

// 在相关系数阵上分解；米和周的量级不作为秩判断依据。
Matrix FactorCovariance(const Matrix& q, vector<double>& scale)
{
    const int n = static_cast<int>(q.size());
    if (!ValidMatrix(q, n, n)) throw runtime_error("Invalid covariance dimensions or values");
    scale.resize(n);
    for (int i = 0; i < n; ++i)
    {
        if (q[i][i] <= 0.0) throw runtime_error("Nonpositive covariance diagonal");
        scale[i] = sqrt(q[i][i]);
    }
    Matrix lower = zeros(n, n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            if (fabs(q[i][j] - q[j][i]) / scale[i] / scale[j] > 1e-9)
                throw runtime_error("Asymmetric covariance");
            double v = (0.5 * q[i][j] + 0.5 * q[j][i]) / scale[i] / scale[j];
            for (int k = 0; k < j; ++k) v -= lower[i][k] * lower[j][k];
            if (i == j)
            {
                if (!isfinite(v) || v <= 1e-12) throw runtime_error("Singular or indefinite covariance");
                lower[i][j] = sqrt(v);
            }
            else lower[i][j] = v / lower[j][j];
        }
    }
    return lower;
}

Matrix SolveCovariance(const Matrix& lower, const vector<double>& scale, const Matrix& rhs)
{
    const int n = static_cast<int>(lower.size());
    const int m = static_cast<int>(rhs[0].size());
    Matrix x = zeros(n, m);
    for (int c = 0; c < m; ++c)
    {
        for (int i = 0; i < n; ++i)
        {
            double v = rhs[i][c] / scale[i];
            for (int j = 0; j < i; ++j) v -= lower[i][j] * x[j][c];
            x[i][c] = v / lower[i][i];
        }
        for (int i = n - 1; i >= 0; --i)
        {
            double v = x[i][c];
            for (int j = i + 1; j < n; ++j) v -= lower[j][i] * x[j][c];
            x[i][c] = v / lower[i][i];
        }
        for (int i = 0; i < n; ++i) x[i][c] /= scale[i];
    }
    if (!ValidMatrix(x, n, m)) throw runtime_error("Nonfinite covariance solution");
    return x;
}

double Sign(double x) { return x <= 0.0 ? -1.0 : 1.0; }
double RoundInteger(double x) { return floor(x + 0.5); }

// 降相关和MLAMBDA搜索沿用课程配套LAMBDA/LAMBDA/lambda.c的算法。
// 内部仍按列存储，a[i+j*n]是第i行第j列；不与项目Matrix的行索引混用。
void LdFactor(int n, const vector<double>& q, vector<double>& lower, vector<double>& d)
{
    auto a = q;
    for (int i = n - 1; i >= 0; --i)
    {
        d[i] = a[i + i * n];
        if (!isfinite(d[i]) || d[i] <= 0.0) throw runtime_error("LAMBDA LD factorization failed");
        const double root = sqrt(d[i]);
        for (int j = 0; j <= i; ++j) lower[i + j * n] = a[i + j * n] / root;
        for (int j = 0; j < i; ++j)
            for (int k = 0; k <= j; ++k) a[j + k * n] -= lower[i + k * n] * lower[i + j * n];
        for (int j = 0; j <= i; ++j) lower[i + j * n] /= lower[i + i * n];
    }
}

void IntegerGauss(int n, vector<double>& lower, vector<double>& z, int i, int j)
{
    const double mu = RoundInteger(lower[i + j * n]);
    if (!isfinite(mu) || fabs(mu) > 1e12) throw runtime_error("LAMBDA integer transform overflow");
    if (mu == 0.0) return;
    for (int k = i; k < n; ++k) lower[k + j * n] -= mu * lower[k + i * n];
    for (int k = 0; k < n; ++k)
    {
        z[k + j * n] -= mu * z[k + i * n];
        if (!isfinite(z[k + j * n]) || fabs(z[k + j * n]) > 1e12)
            throw runtime_error("LAMBDA transform lost integer precision");
    }
}

void Permute(int n, vector<double>& lower, vector<double>& d, int j, double delta, vector<double>& z)
{
    const double eta = d[j] / delta;
    const double lambda = d[j + 1] * lower[j + 1 + j * n] / delta;
    d[j] = eta * d[j + 1];
    d[j + 1] = delta;
    for (int k = 0; k < j; ++k)
    {
        const double a = lower[j + k * n], b = lower[j + 1 + k * n];
        lower[j + k * n] = -lower[j + 1 + j * n] * a + b;
        lower[j + 1 + k * n] = eta * a + lambda * b;
    }
    lower[j + 1 + j * n] = lambda;
    for (int k = j + 2; k < n; ++k) swap(lower[k + j * n], lower[k + (j + 1) * n]);
    for (int k = 0; k < n; ++k) swap(z[k + j * n], z[k + (j + 1) * n]);
}

void Reduce(int n, vector<double>& lower, vector<double>& d, vector<double>& z)
{
    int j = n - 2, last = j, iterations = 0;
    while (j >= 0)
    {
        if (++iterations > 100000) throw runtime_error("LAMBDA reduction limit exceeded");
        if (j <= last) for (int i = j + 1; i < n; ++i) IntegerGauss(n, lower, z, i, j);
        const double v = lower[j + 1 + j * n];
        const double delta = d[j] + v * v * d[j + 1];
        if (!isfinite(delta) || delta <= 0.0) throw runtime_error("Invalid reduced variance");
        if (delta + 1e-6 < d[j + 1])
        {
            Permute(n, lower, d, j, delta, z);
            last = j;
            j = n - 2;
        }
        else --j;
    }
}

Matrix Search(int n, const vector<double>& lower, const vector<double>& d, const vector<double>& mean,
    int limit, int& iterations)
{
    vector<double> sums(n * n), distance(n), conditional(n), integer(n), step(n);
    Matrix candidates = zeros(n, 2);
    double norms[2]{};
    double radius = numeric_limits<double>::infinity();
    int k = n - 1, count = 0, worst = 0;
    conditional[k] = mean[k];
    integer[k] = RoundInteger(conditional[k]);
    double y = conditional[k] - integer[k];
    step[k] = Sign(y);
    for (iterations = 1; iterations <= limit; ++iterations)
    {
        const double next_distance = distance[k] + y * y / d[k];
        if (!isfinite(next_distance)) throw runtime_error("Nonfinite integer search distance");
        if (next_distance < radius)
        {
            if (k > 0)
            {
                distance[--k] = next_distance;
                for (int i = 0; i <= k; ++i)
                    sums[k + i * n] = sums[k + 1 + i * n] + (integer[k + 1] - conditional[k + 1]) * lower[k + 1 + i * n];
                conditional[k] = mean[k] + sums[k + k * n];
                integer[k] = RoundInteger(conditional[k]);
                y = conditional[k] - integer[k];
                step[k] = Sign(y);
            }
            else
            {
                if (count < 2)
                {
                    if (count == 0 || next_distance > norms[worst]) worst = count;
                    for (int i = 0; i < n; ++i) candidates[i][count] = integer[i];
                    norms[count++] = next_distance;
                }
                else if (next_distance < norms[worst])
                {
                    for (int i = 0; i < n; ++i) candidates[i][worst] = integer[i];
                    norms[worst] = next_distance;
                    worst = norms[0] > norms[1] ? 0 : 1;
                }
                if (count == 2) radius = norms[worst];
                integer[0] += step[0];
                y = conditional[0] - integer[0];
                step[0] = -step[0] - Sign(step[0]);
            }
        }
        else
        {
            if (k == n - 1)
            {
                if (count != 2) throw runtime_error("Fewer than two integer candidates");
                return candidates;
            }
            ++k;
            integer[k] += step[k];
            y = conditional[k] - integer[k];
            step[k] = -step[k] - Sign(step[k]);
        }
    }
    throw runtime_error("LAMBDA search iteration limit exceeded");
}

void CopyFloatPosition(const rtk_solution_t& input, rtk_ambiguity_solution_t& result)
{
    result.time = input.time;
    result.position_covariance = zeros(3, 3);
    for (int i = 0; i < 3; ++i)
    {
        result.rover_xyz[i] = input.rover_xyz[i];
        result.baseline_xyz[i] = input.baseline_xyz[i];
        result.baseline_enu[i] = input.baseline_enu[i];
        for (int j = 0; j < 3; ++j) result.position_covariance[i][j] = input.covariance[i][j];
    }
    result.valid = true;
}


bool SearchRtkIntegerAmbiguities(const Matrix& ambiguities, const Matrix& covariance,
    int max_search_iterations, rtk_integer_candidates_t& candidates)
{
    candidates = rtk_integer_candidates_t{};
    try
    {
        const int n = static_cast<int>(ambiguities.size());
        if (!ValidMatrix(ambiguities, n, 1) || !ValidMatrix(covariance, n, n) || max_search_iterations < 1)
            throw runtime_error("Invalid integer search input");
        vector<double> scale;
        const Matrix cholesky = FactorCovariance(covariance, scale);
        vector<double> q(n * n), lower(n * n), d(n), z(n * n), transformed(n);
        for (int i = 0; i < n; ++i)
        {
            if (fabs(ambiguities[i][0]) > 1e12) throw runtime_error("Ambiguity exceeds integer precision limit");
            z[i + i * n] = 1.0;
            for (int j = 0; j < n; ++j) q[i + j * n] = 0.5 * covariance[i][j] + 0.5 * covariance[j][i];
        }
        LdFactor(n, q, lower, d);
        Reduce(n, lower, d, z);
        Matrix z_transpose = zeros(n, n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
            {
                z_transpose[i][j] = z[j + i * n];
                transformed[i] += z_transpose[i][j] * ambiguities[j][0];
            }
        for (double value : transformed)
            if (!isfinite(value) || fabs(value) > 1e12) throw runtime_error("Reduced ambiguity exceeds integer precision limit");
        const Matrix reduced_candidates = Search(n, lower, d, transformed, max_search_iterations, candidates.search_iterations);
        const Matrix inverse_z = inverse(z_transpose);
        if (!ValidMatrix(inverse_z, n, n)) throw runtime_error("LAMBDA inverse transform failed");
        Matrix integers = mul(inverse_z, reduced_candidates);
        for (auto& row : integers)
            for (double& v : row)
            {
                const double rounded = RoundInteger(v);
                if (!isfinite(v) || fabs(v) > 1e12 || fabs(v - rounded) > 1e-5)
                    throw runtime_error("LAMBDA inverse transform lost integer precision");
                v = rounded;
            }
        const Matrix recovered = mul(z_transpose, integers);
        for (int i = 0; i < n; ++i)
            for (int c = 0; c < 2; ++c)
                if (!isfinite(recovered[i][c]) || fabs(recovered[i][c] - reduced_candidates[i][c]) > 1e-5)
                    throw runtime_error("LAMBDA integer back-transformation is inconsistent");
        Matrix differences = zeros(n, 2);
        bool distinct = false;
        for (int i = 0; i < n; ++i)
        {
            distinct = distinct || integers[i][0] != integers[i][1];
            for (int c = 0; c < 2; ++c) differences[i][c] = ambiguities[i][0] - integers[i][c];
        }
        if (!distinct) throw runtime_error("Integer candidates are identical");
        // 回到原模糊度坐标系重算二次型，避免把降相关坐标或矩阵转置错误带入Ratio。
        const Matrix weighted = SolveCovariance(cholesky, scale, differences);
        for (int c = 0; c < 2; ++c)
        {
            double norm = 0.0;
            for (int i = 0; i < n; ++i) norm += differences[i][c] * weighted[i][c];
            if (!isfinite(norm) || norm < 0.0) throw runtime_error("Invalid candidate squared norm");
            candidates.squared_norms[c] = norm;
        }
        if (candidates.squared_norms[0] > candidates.squared_norms[1])
        {
            swap(candidates.squared_norms[0], candidates.squared_norms[1]);
            for (auto& row : integers) swap(row[0], row[1]);
        }
        if (candidates.squared_norms[1] <= 0.0) throw runtime_error("Degenerate candidate distances");
        candidates.integers = std::move(integers);
        candidates.valid = true;
        return true;
    }
    catch (const exception& error)
    {
        const int iterations = candidates.search_iterations;
        candidates = rtk_integer_candidates_t{};
        candidates.search_iterations = iterations;
        candidates.error_message = error.what();
        return false;
    }
}

bool ResolveRtkAmbiguities(const rtk_solution_t& input,
    const rtk_ambiguity_options_t& options, rtk_ambiguity_solution_t& result)
{
    result = rtk_ambiguity_solution_t{};
    result.ratio_threshold = options.ratio_threshold;
    try
    {
        const int nx = input.parameter_count, n = nx - 3;
        if (!input.valid || n < 1 || !ValidMatrix(input.x, nx, 1) ||
            !ValidMatrix(input.covariance, nx, nx)) throw runtime_error("Invalid float solution");
        for (int i = 0; i < 3; ++i)
            if (!isfinite(input.rover_xyz[i]) || !isfinite(input.base_xyz[i]) ||
                !isfinite(input.baseline_xyz[i]) || !isfinite(input.baseline_enu[i]))
                throw runtime_error("Nonfinite float coordinates");
        vector<double> scale;
        FactorCovariance(input.covariance, scale);
        CopyFloatPosition(input, result);
        if (!isfinite(options.ratio_threshold) || options.ratio_threshold <= 1.0 || options.max_search_iterations < 1)
            throw runtime_error("Invalid ambiguity resolution options");
        if (!ValidMatrix(input.ambiguity_covariance, n, n) ||
            !ValidMatrix(input.position_ambiguity_covariance, 3, n) || static_cast<int>(input.ambiguities.size()) != n)
            throw runtime_error("Ambiguity covariance blocks or mapping dimensions disagree");
        set<tuple<int, int, int>> mapping;
        result.float_ambiguities = zeros(n, 1);
        result.ambiguities = input.ambiguities;
        for (int i = 0; i < n; ++i)
        {
            const auto& item = input.ambiguities[i];
            if (item.column != i + 3 || item.sat < 1 || item.sat > MAXSAT ||
                item.reference_sat < 1 || item.reference_sat > MAXSAT || item.sat == item.reference_sat ||
                (item.sys != SYS_GPS && item.sys != SYS_CMP) || item.frequency < 0 || item.frequency > 1 ||
                !isfinite(item.wavelength) || item.wavelength <= 0.0 ||
                !mapping.emplace(item.sat, item.reference_sat, item.frequency).second)
                throw runtime_error("Invalid or duplicate ambiguity column mapping");
            result.float_ambiguities[i][0] = input.x[i + 3][0];
            for (int j = 0; j < nx; ++j)
            {
                const double block = j < 3 ? input.position_ambiguity_covariance[j][i] : input.ambiguity_covariance[i][j - 3];
                if (fabs(block - input.covariance[i + 3][j]) / scale[i + 3] / scale[j] > 1e-9)
                    throw runtime_error("QNN/QXN disagree with full covariance");
            }
        }
        for (const auto& row : input.rows)
        {
            // 奇偶校验未通过的相位仍可参与浮点估计，但此时不作整周约束。
            if (row.type == RtkObservationType::Phase && (row.lli & LLI_HALFC))
            {
                result.error_message = "Unresolved half-cycle flag; retained FLOAT";
                return true;
            }
        }
        if (!SearchRtkIntegerAmbiguities(result.float_ambiguities, input.ambiguity_covariance,
            options.max_search_iterations, result.candidates)) throw runtime_error(result.candidates.error_message);
        result.ratio = result.candidates.squared_norms[0] == 0.0 ? numeric_limits<double>::infinity() :
            result.candidates.squared_norms[1] / result.candidates.squared_norms[0];
        if (!(result.ratio > options.ratio_threshold)) return true;

        vector<double> ambiguity_scale;
        const Matrix lower = FactorCovariance(input.ambiguity_covariance, ambiguity_scale);
        Matrix delta = zeros(n, 1);
        for (int i = 0; i < n; ++i) delta[i][0] = result.candidates.integers[i][0] - result.float_ambiguities[i][0];
        // 条件均值：X_fixed=X_float+QXN*QNN^-1*(N_integer-N_float)。
        const Matrix correction = mul(input.position_ambiguity_covariance, SolveCovariance(lower, ambiguity_scale, delta));
        Matrix fixed_covariance = sub(result.position_covariance, mul(input.position_ambiguity_covariance,
            SolveCovariance(lower, ambiguity_scale, transpose(input.position_ambiguity_covariance))));
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < i; ++j) fixed_covariance[i][j] = fixed_covariance[j][i] =
                0.5 * fixed_covariance[i][j] + 0.5 * fixed_covariance[j][i];
        vector<double> position_scale;
        FactorCovariance(fixed_covariance, position_scale);
        XYZ base{input.base_xyz[0], input.base_xyz[1], input.base_xyz[2]};
        if (hypot(base.X, base.Y) < 1.0) throw runtime_error("Invalid base coordinate for ENU rotation");
        double fixed_rover[3], baseline[3];
        for (int i = 0; i < 3; ++i)
        {
            fixed_rover[i] = input.rover_xyz[i] + correction[i][0];
            baseline[i] = fixed_rover[i] - input.base_xyz[i];
            if (!isfinite(fixed_rover[i]) || !isfinite(baseline[i])) throw runtime_error("Nonfinite fixed coordinates");
        }
        unique_ptr<BLH> blh(XYZtoBLH(base, 6378137.0, 1.0 / 298.257223563));
        unique_ptr<ENU> enu(XYZtoENU(blh.get(), baseline));
        if (!isfinite(enu->E) || !isfinite(enu->N) || !isfinite(enu->U)) throw runtime_error("Nonfinite fixed ENU");
        for (int i = 0; i < 3; ++i)
        {
            result.rover_xyz[i] = fixed_rover[i];
            result.baseline_xyz[i] = baseline[i];
        }
        result.baseline_enu[0] = enu->E;
        result.baseline_enu[1] = enu->N;
        result.baseline_enu[2] = enu->U;
        result.position_covariance = std::move(fixed_covariance);
        result.fixed = true;
        return true;
    }
    catch (const exception& error)
    {
        result.error_message = error.what();
        return false;
    }
}

bool WriteRtkAmbiguitySolution(ostream& out, const rtk_ambiguity_solution_t& solution)
{
    if (!out.good() || !solution.valid) return false;
    const auto flags = out.flags();
    const auto precision = out.precision();
    out << defaultfloat << setprecision(15);
    out << "RTK_AR week=" << solution.time.week << " sow=" << solution.time.sec
        << " status=" << (solution.fixed ? "FIX" : "FLOAT") << " ratio=" << solution.ratio
        << " threshold=" << solution.ratio_threshold << " search_ok=" << solution.candidates.valid
        << " search_iterations=" << solution.candidates.search_iterations << '\n';
    if (!solution.error_message.empty()) out << "AR_MESSAGE " << solution.error_message << '\n';
    out << "ROVER_XYZ_m " << solution.rover_xyz[0] << ' ' << solution.rover_xyz[1] << ' ' << solution.rover_xyz[2] << '\n';
    out << "BASELINE_XYZ_m " << solution.baseline_xyz[0] << ' ' << solution.baseline_xyz[1] << ' ' << solution.baseline_xyz[2] << '\n';
    out << "BASELINE_ENU_m " << solution.baseline_enu[0] << ' ' << solution.baseline_enu[1] << ' ' << solution.baseline_enu[2] << '\n';
    out << "POSITION_COV_m2 3 3\n";
    for (const auto& row : solution.position_covariance)
    {
        for (double value : row) out << value << ' ';
        out << '\n';
    }
    if (solution.candidates.valid)
    {
        out << "SQUARED_NORMS " << solution.candidates.squared_norms[0] << ' ' << solution.candidates.squared_norms[1] << '\n';
        for (int i = 0; i < static_cast<int>(solution.ambiguities.size()); ++i)
        {
            const auto& item = solution.ambiguities[i];
            out << "AMB column=" << item.column << " sat=" << sat2id(item.sat) << " ref=" << sat2id(item.reference_sat)
                << " frequency=" << item.frequency << " float_cycle=" << solution.float_ambiguities[i][0]
                << " best_cycle=" << solution.candidates.integers[i][0] << " second_cycle=" << solution.candidates.integers[i][1] << '\n';
        }
    }
    out << "END_RTK_AR\n";
    out.flags(flags);
    out.precision(precision);
    return out.good();
}

bool WriteRtkAmbiguitySolutionCsvHeader(ostream& out)
{
    out << "week,sow,status,ratio,search_ok,search_iterations,best_squared_norm,second_squared_norm,"
        "X_m,Y_m,Z_m,dX_m,dY_m,dZ_m,E_m,N_m,U_m,std_X_m,std_Y_m,std_Z_m\n";
    return out.good();
}

bool WriteRtkAmbiguitySolutionCsv(ostream& out, const rtk_ambiguity_solution_t& solution)
{
    if (!out.good() || !solution.valid) return false;
    const auto flags = out.flags();
    const auto precision = out.precision();
    out << defaultfloat << setprecision(15);
    out << solution.time.week << ',' << solution.time.sec << ',' << (solution.fixed ? "FIX" : "FLOAT") << ',';
    if (solution.candidates.valid) out << solution.ratio;
    out << ',' << solution.candidates.valid << ',' << solution.candidates.search_iterations << ',';
    if (solution.candidates.valid) out << solution.candidates.squared_norms[0];
    out << ',';
    if (solution.candidates.valid) out << solution.candidates.squared_norms[1];
    for (double v : solution.rover_xyz) out << ',' << v;
    for (double v : solution.baseline_xyz) out << ',' << v;
    for (double v : solution.baseline_enu) out << ',' << v;
    for (int i = 0; i < 3; ++i) out << ',' << sqrt(solution.position_covariance[i][i]);
    out << '\n';
    out.flags(flags);
    out.precision(precision);
    return out.good();
}
