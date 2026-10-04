#include "rtk_mathematical_model.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
using namespace std;

// 任务二仅构建短基线、单历元、非组合双差函数模型
// 数据流：共视历元 -> 站间单差 -> 分系统选参考星 -> 确定行列映射 -> 构建B、L。

constexpr double DEG_TO_RAD = 3.14159265358979323846 / 180.0;

// 与epoch.common一一对应的站间单差，参考星选择和B/L共用这份数据。
// 单差数组保留原共视数组的下标；未通过检查的元素valid=false，不参与后续构模。
struct sd_obs_t
{
    bool valid = false;
    int sat = 0;
    int sys = SYS_NONE;
    double code[2]{};              // P_rover-P_base，m
    double phase[2]{};             // lambda*(L_rover-L_base)，m
    double geometry = 0.0;        // rho_rover-rho_base，m
    double coefficient[3]{};      // 流动站坐标偏导：-(sat-rover)/rho
    double rover_elevation = 0.0;
    unsigned char lli[2]{};
};

// 一个系统对应一个分组：一颗参考星 + 若干非参考星，GPS与BDS不跨系统作差。
struct system_group_t
{
    int sys = SYS_NONE;
    int reference_index = -1;      // epoch.common/sd数组中的下标
    vector<int> target_indices;
};

// 统一失败出口：清空本次构造的中间结果，防止调用者误用上个历元的B/L。
static bool Fail(rtk_function_model_t& model, const string& message)
{
    model = rtk_function_model_t{};
    model.error_message = message;
    return false;
}

// 检查GPS周与周内秒的基本范围；具体两站是否同步由后续timediff判断。
static bool ValidTime(const gtime_t& time)
{
    return time.week >= 0 && isfinite(time.sec) && time.sec >= 0.0 && time.sec < 604800.0;
}

// 排除空指针、NaN/无穷大以及接近原点的无效坐标；坐标单位为米。
static bool ValidXYZ(const double* xyz)
{
    return xyz != nullptr && isfinite(xyz[0]) &&
        isfinite(xyz[1]) && isfinite(xyz[2]) &&
        hypot(xyz[0], xyz[1], xyz[2]) > 1.0;
}

// 频点下标0/1：GPS为L1/L2，BDS为B1I/B3I；波长lambda=c/f，单位m/cycle。
static double Wavelength(int sys, int frequency)
{
    if (sys == SYS_GPS)
    {
        return Clight / (frequency == 0 ? FREQ_GPS_L1 : FREQ_GPS_L2);
    }
    return Clight / (frequency == 0 ? FREQ_BDS_B1 : FREQ_BDS_B3);
}

// 确认观测码属于当前系统和频点，本模型不混用其他频点或组合观测。
static bool SupportedCode(int sys, int frequency, unsigned char code)
{
    if (sys == SYS_GPS)
    {
        if (frequency == 0)
            return code == CODE_L1C || code == CODE_L1P;
        return code == CODE_L2P || code == CODE_L2C || code == CODE_L2W;
    }
    return code == (frequency == 0 ? CODE_L2I : CODE_L6I);
}

// 一个频点必须在两站同时具备有效伪距、相位和足够SNR，且观测码一致。
// 单频模式只检查下标0；双频模式在调用处分别检查下标0和1。
static bool ValidFrequency(const obsd_t& base, const obsd_t& rover,int sys, int frequency, double snr_mask)
{
    return SupportedCode(sys, frequency, base.code[frequency]) &&
        base.code[frequency] == rover.code[frequency] &&
        isfinite(base.P[frequency]) && base.P[frequency] > 0.0 &&
        isfinite(rover.P[frequency]) && rover.P[frequency] > 0.0 &&
        isfinite(base.L[frequency]) && base.L[frequency] != 0.0 &&
        isfinite(rover.L[frequency]) && rover.L[frequency] != 0.0 &&
        isfinite(base.SNR[frequency]) && base.SNR[frequency] > snr_mask &&
        isfinite(rover.SNR[frequency]) && rover.SNR[frequency] > snr_mask;
}

// 用各站接收机坐标和该站保存的卫星坐标计算几何距离rho，单位m。
static double GeometricRange(const double receiver[3], const double satellite[3])
{
    return hypot(satellite[0] - receiver[0],satellite[1] - receiver[1], satellite[2] - receiver[2]);
}

// GEO分类沿用本项目satpos.cpp的PRN规则；GEO仍可作为非参考星。
static bool IsBdsGeo(int sat)
{
    int prn = 0;
    return satsys(sat, &prn) == SYS_CMP && (prn <= 5 || prn >= 59);
}

// 第一步：逐颗检查共视观测，计算站间单差和流动站坐标偏导。
// 结构/时刻/坐标错误返回false；高度角或观测质量不达标则跳过该颗卫星。
static bool ComputeSingleDifference(const rtk_epoch_t& epoch,const double rover_xyz[3], const rtk_config_t& config,const rtk_model_options_t& options, vector<sd_obs_t>& sd,string& error)
{
    // 不压缩数组，保证sd[i]始终对应epoch.common[i]。
    sd.resize(epoch.ncommon);
    bool seen[MAXSAT + 1]{};
    const double elevation_mask = config.elevation_mask * DEG_TO_RAD;

    for (int i = 0; i < epoch.ncommon; ++i)
    {
        // 两站观测数组顺序可能不同，必须通过各自下标取得同一颗卫星。
        const common_sat_t& item = epoch.common[i];
        if (item.sat < 1 || item.sat > MAXSAT || seen[item.sat] ||item.base_index < 0 || item.base_index >= epoch.base_obs.n ||item.rover_index < 0 || item.rover_index >= epoch.rover_obs.n)
        {
            error = "Invalid or duplicate common satellite/index";
            return false;
        }
        seen[item.sat] = true;
        const obsd_t& base = epoch.base_obs.data[item.base_index];
        const obsd_t& rover = epoch.rover_obs.data[item.rover_index];
        if (base.sat != item.sat || rover.sat != item.sat ||item.base_sat.sat != item.sat || item.rover_sat.sat != item.sat)
        {
            error = "Observation and satellite-state indices do not match";
            return false;
        }

        int prn = 0;
        const int sys = satsys(item.sat, &prn);
        if ((sys != SYS_GPS && sys != SYS_CMP) || (sys == SYS_GPS && !options.use_gps) || (sys == SYS_CMP && !options.use_bds))
        {
            continue;
        }

        if (!ValidTime(base.time) || !ValidTime(rover.time) ||
            fabs(timediff(base.time, epoch.time)) >= config.sync_tolerance ||
            fabs(timediff(rover.time, epoch.rover_obs.data[0].time)) >= config.sync_tolerance ||
            fabs(timediff(base.time, rover.time)) >= config.sync_tolerance)
        {
            error = "Common satellite observation times are not synchronized";
            return false;
        }
        if (!ValidXYZ(item.base_sat.pos) || !ValidXYZ(item.rover_sat.pos))
        {
            error = "Invalid satellite coordinates";
            return false;
        }
        // 复用任务一已保存的各站高度角，不重新实现坐标转换。
        if (!isfinite(item.base_azel[1]) || !isfinite(item.rover_azel[1]) || item.base_azel[1] < elevation_mask || item.rover_azel[1] < elevation_mask)
        {
            continue;

        }

        // 所选频点全部通过检查后，才让该卫星参与本模型。
        bool frequencies_ok = true;
        for (int f = 0; f < options.frequency_count; ++f)
        {
            frequencies_ok = frequencies_ok && ValidFrequency(base, rover, sys, f, config.snr_mask);

        }
        if (!frequencies_ok)
        {
            continue;

        }

        // 基准站坐标在本模型中视为已知量，只估计流动站坐标改正。
        const double base_range = GeometricRange(epoch.base_spp.XYZ, item.base_sat.pos);
        const double rover_range = GeometricRange(rover_xyz, item.rover_sat.pos);
        if (!isfinite(base_range) || !isfinite(rover_range) ||base_range <= 1.0 || rover_range <= 1.0)
        {
            error = "Invalid receiver-to-satellite geometric range";
            return false;
        }
        sd_obs_t& result = sd[i];
        result.valid = true;
        result.sat = item.sat;
        result.sys = sys;
        result.geometry = rover_range - base_range;
        result.rover_elevation = item.rover_azel[1];
        // rho=|sat-rover|，故对rover坐标求导得到负的视线单位向量。
        for (int axis = 0; axis < 3; ++axis)
        {
            result.coefficient[axis] = -(item.rover_sat.pos[axis] - rover_xyz[axis]) / rover_range;

        }
        for (int f = 0; f < options.frequency_count; ++f)
        {
            result.code[f] = rover.P[f] - base.P[f];
            // 原始L以周为单位，乘波长后才能与伪距、几何距离统一相减。
            result.phase[f] = Wavelength(sys, f) * (rover.L[f] - base.L[f]);
            // 合并两站的失锁标志，仅保留信息；这里不进行周跳修复。
            result.lli[f] = base.LLI[f] | rover.LLI[f];
        }
    }
    return true;
}

// 第二步：在同一系统的合格卫星中选流动站高度角最高的一颗作为参考星。默认不允许BDS GEO作参考星，但GEO可以保留为参与双差的非参考星。
static int SelectReferenceSatellite(const vector<int>& indices,const vector<sd_obs_t>& sd, const rtk_model_options_t& options)
{
    int best = -1;
    for (int i : indices)
    {
        if (options.exclude_bds_geo_reference && IsBdsGeo(sd[i].sat))
            continue;
        if (best < 0 || sd[i].rover_elevation > sd[best].rover_elevation)
            best = i;
    }
    return best; // 同高度角时保留已按sat排序的第一颗
}

// 按GPS、BDS顺序组织单差数据，卫星按内部sat编号排序，使输出顺序稳定。target_indices存非参考星下标，reference_index存该系统唯一的参考星下标。
static bool BuildGroups(const vector<sd_obs_t>& sd,const rtk_model_options_t& options, vector<system_group_t>& groups,string& error)
{
    for (int sys : { SYS_GPS, SYS_CMP })
    {
        vector<int> indices;
        for (int i = 0; i < static_cast<int>(sd.size()); ++i)
        {
            if (sd[i].valid && sd[i].sys == sys)
            {
                indices.push_back(i);
            }
        }
            
        sort(indices.begin(), indices.end(), [&](int a, int b) { return sd[a].sat < sd[b].sat; });
        // 一个系统不足两颗星时不能产生双差，省去该系统。
        if (indices.size() < 2)
        {
            continue;

        }
        system_group_t group;
        group.sys = sys;
        group.reference_index = SelectReferenceSatellite(indices, sd, options);
        if (group.reference_index < 0)
        {
            error = "No eligible BDS reference satellite (all candidates are GEO)";
            return false;
        }
        for (int i : indices)
        {
            if (i != group.reference_index)
            {
                group.target_indices.push_back(i);

            }
        }
            
        groups.push_back(group);
    }
    return true;
}

// 伪距行决定三维坐标能否独立估计；相位行另含自由模糊度参数。
// 用双差坐标系数组成三列矩阵，通过带选主元的消元检查秩是否为3。
// 至少三个独立双差方向才能估计dX/dY/dZ；仅满足观测数量不代表几何满秩。
static bool HasThreeDimensionalGeometry(const vector<system_group_t>& groups,const vector<sd_obs_t>& sd)
{
    Matrix geometry;
    for (const system_group_t& group : groups)
    {
        for (int i : group.target_indices)
        {
            vector<double> row(3);
            for (int axis = 0; axis < 3; ++axis)
            {
                row[axis] = sd[i].coefficient[axis] - sd[group.reference_index].coefficient[axis];

            }
            geometry.push_back(row);
        }
    }
        
    // 每找到一个非零主元就增加一个独立方向；容差用于排除数值近零。
    int rank = 0;
    for (int column = 0; column < 3 && rank < static_cast<int>(geometry.size()); ++column)
    {
        int pivot = rank;
        for (int row = rank + 1; row < static_cast<int>(geometry.size()); ++row)
        {
            if (fabs(geometry[row][column]) > fabs(geometry[pivot][column]))
            {
                pivot = row;
            }
        }
                
        if (fabs(geometry[pivot][column]) <= 1e-10)
        {
            continue;
        }
        swap(geometry[rank], geometry[pivot]);
        const double divisor = geometry[rank][column];
        for (int k = column; k < 3; ++k)
        {
            geometry[rank][k] /= divisor;

        }
        for (int row = rank + 1; row < static_cast<int>(geometry.size()); ++row)
        {
            const double factor = geometry[row][column];
            for (int k = column; k < 3; ++k)
                geometry[row][k] -= factor * geometry[rank][k];
        }
        ++rank;
    }
    return rank == 3;
}

// 第三步：先确定每条观测的行、每个模糊度的列，再用同一映射填B和L。
// 每个系统内的行顺序：相位f0、相位f1、伪距f0、伪距f1（单频省略f1）。
// 列0/1/2为dX/dY/dZ；其后为GPS各频点模糊度，再为BDS各频点模糊度。
static void BuildRowAndColumnLayout(const vector<system_group_t>& groups,const vector<sd_obs_t>& sd, int frequency_count, rtk_function_model_t& model)
{
    for (const system_group_t& group : groups)
    {
        const int pair_count = static_cast<int>(group.target_indices.size());
        const int reference_sat = sd[group.reference_index].sat;
        // 跳过前三个坐标参数，以及先前系统已经占用的模糊度列。
        const int ambiguity_start = 3 + static_cast<int>(model.ambiguities.size());
        if (group.sys == SYS_GPS)
        {
            model.gps_reference_sat = reference_sat;
            model.gps_satellite_count = pair_count + 1;
        }
        else
        {
            model.bds_reference_sat = reference_sat;
            model.bds_satellite_count = pair_count + 1;
        }
        model.double_difference_count += pair_count;
        for (int f = 0; f < frequency_count; ++f)
        {
            for (int k = 0; k < pair_count; ++k)
            {
                dd_ambiguity_info_t ambiguity;
                // 每个频点占pair_count列，每颗非参考星对应一个双差模糊度。
                ambiguity.column = ambiguity_start + f * pair_count + k;
                ambiguity.sys = group.sys;
                ambiguity.frequency = f;
                ambiguity.sat = sd[group.target_indices[k]].sat;
                ambiguity.reference_sat = reference_sat;
                ambiguity.wavelength = Wavelength(group.sys, f);
                model.ambiguities.push_back(ambiguity);
            }
        }
            
        for (RtkObservationType type : { RtkObservationType::Phase, RtkObservationType::Code })
            for (int f = 0; f < frequency_count; ++f)
                for (int k = 0; k < pair_count; ++k)
                {
                    const int target = group.target_indices[k];
                    dd_row_info_t row;
                    row.sys = group.sys;
                    row.frequency = f;
                    row.sat = sd[target].sat;
                    row.reference_sat = reference_sat;
                    row.type = type;
                    row.ambiguity_column = type == RtkObservationType::Phase ? ambiguity_start + f * pair_count + k : -1;
                    row.common_index = target;
                    row.reference_common_index = group.reference_index;
                    row.wavelength = Wavelength(group.sys, f);
                    // 合并两站、非参考星/参考星共四条原始相位的LLI。
                    row.lli = sd[target].lli[f] | sd[group.reference_index].lli[f];
                    model.rows.push_back(row);
                }
    }
    model.observation_count = static_cast<int>(model.rows.size());
    model.parameter_count = 3 + static_cast<int>(model.ambiguities.size());
}

// 第四步：填设计矩阵B（观测数nv行、未知参数数nx列）。
// 每行前三列是双差几何距离对流动站XYZ的偏导；相位行另有一个+lambda。
// 伪距不含模糊度，因此伪距行在第3列及其后全部为0（列号从0开始）。
static void BuildDesignMatrix(const vector<sd_obs_t>& sd, rtk_function_model_t& model)
{
    model.B = zeros(model.observation_count, model.parameter_count);
    for (int r = 0; r < model.observation_count; ++r)
    {
        const dd_row_info_t& row = model.rows[r];
        for (int axis = 0; axis < 3; ++axis)
        {
            model.B[r][axis] = sd[row.common_index].coefficient[axis] -
                sd[row.reference_common_index].coefficient[axis];
        }
            
        // 相位方程含+lambda*N；其他模糊度列保持初始化的0。
        if (row.type == RtkObservationType::Phase)
        {
            model.B[r][row.ambiguity_column] = row.wavelength;
        }
    }
}

// 第五步：填L=O-C，即双差观测量减双差几何距离，大小为nv×1。
// 相位已乘波长转换为米；不扣除模糊度初值，因为未知数是绝对双差N。
static void BuildResidualVector(const vector<sd_obs_t>& sd, rtk_function_model_t& model)
{
    model.L = zeros(model.observation_count, 1);
    for (int r = 0; r < model.observation_count; ++r)
    {
        const dd_row_info_t& row = model.rows[r];
        const sd_obs_t& satellite = sd[row.common_index];
        const sd_obs_t& reference = sd[row.reference_common_index];
        const int f = row.frequency;
        // 对站间单差再作星间差：目标星单差减参考星单差。
        const double observed = row.type == RtkObservationType::Phase ?
            satellite.phase[f] - reference.phase[f] : satellite.code[f] - reference.code[f];
        // 双差几何项也使用同样顺序，保证B与L的符号约定一致。
        model.L[r][0] = observed - (satellite.geometry - reference.geometry);
    }
}

// 构模和导出前的自检：维数、数值、行列映射及模糊度系数必须一致。
// 设双差卫星对数为D、频点数为F：nv=2*F*D，nx=3+F*D。
static bool ValidModel(const rtk_function_model_t& model)
{
    if (!model.error_message.empty() || model.double_difference_count < 3 ||
        (model.frequency_count != 1 && model.frequency_count != 2) ||
        model.observation_count != 2 * model.frequency_count * model.double_difference_count ||
        model.parameter_count != 3 + model.frequency_count * model.double_difference_count ||
        model.observation_count < model.parameter_count ||
        static_cast<int>(model.B.size()) != model.observation_count ||
        static_cast<int>(model.L.size()) != model.observation_count ||
        static_cast<int>(model.rows.size()) != model.observation_count ||
        static_cast<int>(model.ambiguities.size()) != model.parameter_count - 3)
    {
        return false;

    }
    for (int c = 3; c < model.parameter_count; ++c)
    {
        const dd_ambiguity_info_t& ambiguity = model.ambiguities[c - 3];
        if (ambiguity.column != c || !isfinite(ambiguity.wavelength) || ambiguity.wavelength <= 0.0)
        {
            return false;

        }
    }
    for (int r = 0; r < model.observation_count; ++r)
    {
        const dd_row_info_t& row = model.rows[r];
        if (static_cast<int>(model.B[r].size()) != model.parameter_count ||
            model.L[r].size() != 1 || !isfinite(model.L[r][0]) ||
            row.frequency < 0 || row.frequency >= model.frequency_count ||
            row.sat == row.reference_sat || row.common_index < 0 || row.reference_common_index < 0)
        {
            return false;

        }
        if (row.type == RtkObservationType::Phase)
        {
            if (row.ambiguity_column < 3 || row.ambiguity_column >= model.parameter_count)
            {
                return false;

            }
            const dd_ambiguity_info_t& ambiguity = model.ambiguities[row.ambiguity_column - 3];
            if (ambiguity.sys != row.sys || ambiguity.frequency != row.frequency ||
                ambiguity.sat != row.sat || ambiguity.reference_sat != row.reference_sat ||
                ambiguity.wavelength != row.wavelength)
            {
                return false;

            }
        }
        else if (row.type != RtkObservationType::Code || row.ambiguity_column != -1)
        {
            return false;

        }
        for (int c = 0; c < model.parameter_count; ++c)
        {
            if (!isfinite(model.B[r][c]) || (c >= 3 && model.B[r][c] != (c == row.ambiguity_column ? row.wavelength : 0.0)))
            {
                return false;

            }
        }
            
    }
    return true;
}

// 以下两个小函数仅用于输出可读标签，不参与数值计算。
static const char* SystemName(int sys)
{
    return sys == SYS_GPS ? "GPS" : "BDS";
}

static const char* FrequencyName(int sys, int frequency)
{
    if (sys == SYS_GPS)
    {
        return frequency == 0 ? "L1" : "L2";

    }
    return frequency == 0 ? "B1I" : "B3I";
}

// 保留原四参数接口,默认双系统双频模式。
bool BuildRtkFunctionModel(const rtk_epoch_t& epoch, const double rover_xyz[3],const rtk_config_t& config, rtk_function_model_t& model)
{
    return BuildRtkFunctionModel(epoch, rover_xyz, config, rtk_model_options_t{}, model);
}

// 对外入口：epoch是任务一准备好的输入，rover_xyz是本次线性化展开点。
// 本函数不修改输入epoch；成功返回true，失败返回false并填写model.error_message。
bool BuildRtkFunctionModel(const rtk_epoch_t& epoch, const double rover_xyz[3],
    const rtk_config_t& config, const rtk_model_options_t& options, rtk_function_model_t& model)
{
    // 0. 清空输出并检查配置、数量、SPP状态及历元同步，避免非法访问。
    model = rtk_function_model_t{};
    if ((!options.use_gps && !options.use_bds) ||
        (options.frequency_count != 1 && options.frequency_count != 2) ||
        !isfinite(config.sync_tolerance) || config.sync_tolerance <= 0.0 ||
        !isfinite(config.elevation_mask) || config.elevation_mask < 0.0 || config.elevation_mask >= 90.0 ||
        !isfinite(config.snr_mask) || config.snr_mask < 0.0)
    {
        return Fail(model, "Invalid system/frequency/quality configuration");

    }
    if (epoch.ncommon <= 0 || epoch.ncommon > MAXSAT ||
        epoch.base_obs.n <= 0 || epoch.base_obs.n > MAXOBS ||
        epoch.rover_obs.n <= 0 || epoch.rover_obs.n > MAXOBS)
    {
        return Fail(model, "Invalid observation or common satellite count");

    }
    if (epoch.base_spp.stat != 1 || epoch.rover_spp.stat != 1 ||
        !ValidXYZ(epoch.base_spp.XYZ) || !ValidXYZ(rover_xyz))
    {
        return Fail(model, "SPP/linearization coordinates are not valid");

    }
    if (!ValidTime(epoch.time) || !ValidTime(epoch.base_obs.data[0].time) ||
        !ValidTime(epoch.rover_obs.data[0].time) ||
        fabs(timediff(epoch.base_obs.data[0].time, epoch.time)) >= config.sync_tolerance ||
        fabs(timediff(epoch.base_obs.data[0].time, epoch.rover_obs.data[0].time)) >= config.sync_tolerance)
    {
        return Fail(model, "Epoch times are not synchronized");

    }

    // 1. 形成站间单差；2. 分系统选择参考星，并检查三维几何是否可解。
    vector<sd_obs_t> sd;
    string error;
    if (!ComputeSingleDifference(epoch, rover_xyz, config, options, sd, error))
        return Fail(model, error);
    vector<system_group_t> groups;
    if (!BuildGroups(sd, options, groups, error))
        return Fail(model, error);
    if (!HasThreeDimensionalGeometry(groups, sd))
        return Fail(model, "Too few independent double differences for three coordinates");

    // 3. 保存本历元信息；4. 先建立映射，再按映射填充B和L并自检。
    model.time = epoch.time;
    model.time_diff = timediff(epoch.base_obs.data[0].time, epoch.rover_obs.data[0].time);
    model.frequency_count = options.frequency_count;
    copy_n(epoch.base_spp.XYZ, 3, model.base_xyz);
    copy_n(rover_xyz, 3, model.rover_xyz);
    BuildRowAndColumnLayout(groups, sd, options.frequency_count, model);
    BuildDesignMatrix(sd, model);
    BuildResidualVector(sd, model);
    if (!ValidModel(model))
        return Fail(model, "Function model dimensions or numeric values are not consistent");
    return true;
}

// 导出一整个历元的模型，out可以是文件流，也可以是控制台流。
// 先输出参考星和行列含义，再输出B/L数值，便于逐行检查；不输出任务三权阵。
bool WriteRtkFunctionModel(ostream& out, const rtk_function_model_t& model)
{
    if (!out.good() || !ValidModel(model))
    {
        return false;

    }
    // 临时使用15位有效数字；结束时恢复调用者原有的流格式。
    const auto old_flags = out.flags();
    const auto old_precision = out.precision();
    out << defaultfloat << setprecision(15);
    out << "FUNCTION_MODEL week=" << model.time.week << " sow=" << model.time.sec
        << " dt_s=" << model.time_diff << " nv=" << model.observation_count
        << " nx=" << model.parameter_count << " dof=" << model.observation_count - model.parameter_count
        << " frequencies=" << model.frequency_count << " dd_pairs=" << model.double_difference_count << '\n';
    out << "CONVENTION SD=rover-base DD=sat-reference L_unit=m N_unit=cycle\n";
    out << "BASE_XYZ_m " << model.base_xyz[0] << ' ' << model.base_xyz[1] << ' ' << model.base_xyz[2] << '\n';
    out << "ROVER_LINEARIZATION_XYZ_m " << model.rover_xyz[0] << ' ' << model.rover_xyz[1] << ' ' << model.rover_xyz[2] << '\n';
    out << "REFERENCE GPS=" << (model.gps_reference_sat ? sat2id(model.gps_reference_sat) : "NONE")
        << " BDS=" << (model.bds_reference_sat ? sat2id(model.bds_reference_sat) : "NONE")
        << " gps_sats=" << model.gps_satellite_count << " bds_sats=" << model.bds_satellite_count << '\n';
    // PARAM解释每一列未知数；ROW解释每一行来自哪颗星、哪种观测。
    out << "PARAM 0 dX_m\nPARAM 1 dY_m\nPARAM 2 dZ_m\n";
    for (const dd_ambiguity_info_t& ambiguity : model.ambiguities)
    {
        out << "PARAM " << ambiguity.column << " N_DD_cycle " << SystemName(ambiguity.sys)
            << ' ' << FrequencyName(ambiguity.sys, ambiguity.frequency) << " sat=" << sat2id(ambiguity.sat)
            << " ref=" << sat2id(ambiguity.reference_sat) << " wavelength_m=" << ambiguity.wavelength << '\n';
    }
        
    for (int r = 0; r < model.observation_count; ++r)
    {
        const dd_row_info_t& row = model.rows[r];
        out << "ROW " << r << ' ' << SystemName(row.sys) << ' ' << FrequencyName(row.sys, row.frequency)
            << ' ' << (row.type == RtkObservationType::Phase ? "PHASE" : "CODE")
            << " sat=" << sat2id(row.sat) << " ref=" << sat2id(row.reference_sat)
            << " ambiguity_col=" << row.ambiguity_column << " lli=" << static_cast<int>(row.lli) << '\n';
    }
    // 矩阵头记录行列数，随后逐行写数值；END标记本历元模型结束。
    out << "B " << model.observation_count << ' ' << model.parameter_count << '\n';
    for (const auto& row : model.B)
    {
        for (int c = 0; c < model.parameter_count; ++c)
            out << (c == 0 ? "" : " ") << row[c];
        out << '\n';
    }
    out << "L " << model.observation_count << " 1\n";
    for (const auto& row : model.L)
    {
        out << row[0] << '\n';

    }
    out << "END_FUNCTION_MODEL\n";
    out.flags(old_flags);
    out.precision(old_precision);
    return out.good();
}
