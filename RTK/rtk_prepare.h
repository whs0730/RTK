#pragma once

#include "decode.h"
#include "obs.h"

// RTK任务一配置参数
struct rtk_config_t
{
    double sync_tolerance = 0.001;  // 时间同步阈值，实际判断使用 < 0.001 s
    double elevation_mask = 10.0;   // 截止高度角，单位：度
    double snr_mask = 30.0;         // 信噪比必须 > 30 dB-Hz
};

// 单个测站的运行状态
struct station_state_t
{
    raw_t raw;                      // 解码状态，内部包含obs和nav

    satpos_t sats[MAXOBS];          // 当前历元卫星位置与钟差
    sol_t spp_sol;                  // 当前历元SPP结果
    sol_t last_spp_sol;             // 上一个成功历元SPP结果

    int valid_sat = 0;
    bool has_last_spp = false;
    bool eof = false;
};

// 一颗共视卫星的数据
struct common_sat_t
{
    int sat = 0;                    // 内部卫星编号

    int base_index = -1;            // 在基准站obs数组中的下标
    int rover_index = -1;           // 在流动站obs数组中的下标

    satpos_t base_sat;              // 基准站对应的卫星状态
    satpos_t rover_sat;             // 流动站对应的卫星状态

    double base_azel[2] = { 0.0 };    // 基准站方位角、高度角
    double rover_azel[2] = { 0.0 };   // 流动站方位角、高度角
};

// 一个已经同步的RTK历元
struct rtk_epoch_t
{
    gtime_t time;                   // 同步历元时间
    double time_diff = 0.0;         // 两站历元时间差

    obs_t base_obs;                 // 基准站当前历元观测值
    obs_t rover_obs;                // 流动站当前历元观测值

    sol_t base_spp;                 // 基准站SPP结果
    sol_t rover_spp;                // 流动站SPP结果

    common_sat_t common[MAXSAT];     // 共视卫星
    int ncommon = 0;
    int common_gps = 0;
    int common_bds = 0;
}; 
// 从一个测站的OEM文件中读取下一个观测历元。
// 读取过程中遇到星历报文时，自动更新station.raw.nav。
bool ReadNextObservationEpoch(
    FILE* fp,
    station_state_t& station
);

// 比较两个测站当前观测历元。
// 返回-1：基准站较早；0：同步；1：流动站较早。
int CompareStationEpochTime(
    const station_state_t& base,
    const station_state_t& rover,
    const rtk_config_t& config,
    double& dt
);

// 计算一个测站的卫星状态并进行SPP。
bool ProcessStationEpoch(
    station_state_t& station
);

// 将两个同步历元保存到rtk_epoch_t，并选取共视卫星。
bool BuildRtkEpoch(
    const station_state_t& base,
    const station_state_t& rover,
    const rtk_config_t& config,
    rtk_epoch_t& epoch
);

// 检查一颗卫星是否满足RTK数据准备要求。
bool PassRtkObservationCheck(
    const obsd_t& base_obs,
    const obsd_t& rover_obs,
    const satpos_t& base_sat,
    const satpos_t& rover_sat,
    const sol_t& base_sol,
    const sol_t& rover_sol,
    const rtk_config_t& config,
    double base_azel[2],
    double rover_azel[2]
);