#pragma once

#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <string>

#include "obs.h"
// 解码模块：根据NovAtel解码观测数据和卫星星历

// NovAtel OEM4 原始数据状态
// 保存当前报文缓冲区、已解码观测值、导航星历和跨历元辅助状态
class raw_t {
public:
    gtime_t time;                               // 当前报文时间
    gtime_t tobs[MAXSAT][NFREQ + NEXOBS];       // 各卫星各频点上一条观测时间
    double lockt[MAXSAT][NFREQ + NEXOBS];       // 锁定时间，用于辅助周跳判断
    uint8_t halfc[MAXSAT][NFREQ + NEXOBS];      // 半周模糊度标记缓存
    obs_t obs;                                  // 当前历元观测值
    nav_t nav;                                  // 当前导航星历
    uint8_t buff[MAXRAWLEN];                    // 当前 OEM4 报文缓冲区
    int nbyte;                                  // 已接收字节数
    int len;                                    // 当前报文总长度
    int ephsat;                                 // 最近更新星历的卫星编号
    int ephset;                                 // 星历更新标记
    char opt[256];                              // 预留选项

    raw_t() {
        nbyte = 0;
        len = 0;
        ephsat = 0;
        ephset = 0;
        memset(buff, 0, sizeof(buff));
        memset(tobs, 0, sizeof(tobs));
        memset(lockt, 0, sizeof(lockt));
        memset(halfc, 0, sizeof(halfc));
        memset(&obs, 0, sizeof(obs_t));
        memset(&nav, 0, sizeof(nav_t));
        opt[0] = '\0';
    }

    // 新报文时间与当前观测历元不同时，清空观测集合，开始缓存新历元。
    void reset_obs() {
        if (obs.n > 0 && fabs(timediff(obs.data[0].time, time)) > 1E-9) {
            obs.n = 0;
        }
    }
};



// 按内部卫星编号查找当前导航数据中的广播星历。
eph_t* find_eph(nav_t* nav, int sat);

// 判断星历槽是否已经解码出可用的轨道参数。
bool has_valid_eph(const eph_t* eph);

// 从小端字节流读取 2 字节无符号整数，供主程序读取报文类型。
uint16_t U2(uint8_t* p);

// 解码一个完整 OEM4 报文，返回值表示是否得到观测/星历等有效数据。
int decode_oem4(raw_t* raw);

// OEM4 同步头检测，用于从字节流中找报文起点。
int sync_oem4(uint8_t* buff, uint8_t data);

// 文件模式输入：逐字节读取并在凑齐完整报文后调用 decode_oem4()。
int input_oem4f(raw_t* raw, FILE* fp);
