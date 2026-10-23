#pragma once
#include "spp.h"
#include "decode.h"
#include "obs.h"
#include "station_satpos.h"
// 两站同步与共视筛选共用的门限；时间使用GPST，角度配置与保存值的单位不同。
struct rtk_config_t
{
    double sync_tolerance = 0.001;  // 两站接收时刻差的绝对值必须严格小于此值，s
    double elevation_mask = 10.0;   // 截止高度角，单位：度
    double snr_mask = 30.0;         // 信噪比门限，dB-Hz；观测值必须严格大于此值
};

// 每站独立持有解码历史，不能让两站共用raw_t中的锁定时间、LLI或导航缓存。
// 当前历元失败时只清空当前结果，last_spp_sol仍可作为以后历元的SPP初值。
struct station_state_t
{
    raw_t raw;                      // 解码状态，内部包含obs和nav

    satpos_t sats[MAXOBS];          // 当前历元卫星位置与钟差
    sol_t spp_sol;                  // 当前历元SPP结果
    sol_t last_spp_sol;             // 上一个成功历元SPP结果
    solvel_t speed_sol;       // 当前历元速度和接收机钟速
    
    int valid_sat = 0;

    bool speed_ok = false;
    bool position_ok = false;
    bool has_last_spp = false;
    bool eof = false;
};

// 共视下标指向本历元两站各自的观测数组，不是内部sat编号，也不是压缩后的卫星状态下标。
struct common_sat_t
{
    int sat = 0;                    // 内部卫星编号

    int base_index = -1;            // 在基准站obs数组中的下标
    int rover_index = -1;           // 在流动站obs数组中的下标

    satpos_t base_sat;              // 基准站对应的卫星状态
    satpos_t rover_sat;             // 流动站对应的卫星状态

    double base_azel[2] = { 0.0 };    // 基准站方位角、高度角，rad；以该站SPP坐标计算
    double rover_azel[2] = { 0.0 };   // 流动站方位角、高度角，rad；以该站SPP坐标计算
};

// 保存两站观测和SPP的值副本；common只列出筛选通过的卫星，obs数组本身不压缩。
// 卫星状态采用共同星历计算，SPP仍是各站独立导航缓存得到的结果。
struct rtk_epoch_t
{
    gtime_t time;                   // 基准站接收时刻，GPST周与周内秒
    double time_diff = 0.0;         // 基准站减流动站的接收时刻差，s

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
// 调用前应保证两站都有观测；空观测也返回0，不能据此判定同步。
int CompareStationEpochTime(
    const station_state_t& base,
    const station_state_t& rover,
    const rtk_config_t& config,
    double& dt
);

// 按该站观测计算卫星状态、SPP及多普勒速度；位置成功即返回true，速度允许失败。
// SPP成功后更新该站的last_spp_sol，不使用另一站或已知基站坐标替代SPP。
bool ProcessStationEpoch(
    station_state_t& station
);

// 两站用同一套有效广播星历，分别重算发射时刻卫星状态并筛选共视星。
// 各站原有观测、星历、SPP和卫星状态缓存不修改；没有可用共同星历的卫星不参与RTK。
// 返回true只说明至少有一颗合格共视星；三维双差几何是否可解由函数模型另行检查。
bool BuildRtkEpoch(
    const station_state_t& base,
    const station_state_t& rover,
    const rtk_config_t& config,
    rtk_epoch_t& epoch
);

// 两站都必须具备前两个频点的有效伪距、相位和同码观测，SNR严格大于门限。
// 准备层固定按双频筛选；后续选择单频模型不会恢复在此处被剔除的卫星。
// LLI在此处不作整周固定判据，交给估计与模糊度模块处理；高度角等于门限时可通过。
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
// 按课程参考格式导出共视星列表和两站SPP摘要；BLH角度用度。
bool WritePreparedEpoch(ostream& out,const rtk_epoch_t& epoch);
// 按卫星、频点列出已解码的观测值，不筛共视星。
bool WriteStationObservations(ostream& out, const obs_t& obs);
// 返回是否读到下一历元；写入失败仅置output_failed，调用方需另行检查该标志。
bool ReadAndExportNext(FILE* fp,station_state_t& station,ostream& out,bool& output_failed);
