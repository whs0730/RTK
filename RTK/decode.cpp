#include "decode.h"
#include <iostream>
// 字节序转换和数据解析工具：NovAtel 二进制字段采用小端存储。
static uint8_t U1(const uint8_t* p) { return *p; }
uint16_t U2(uint8_t* p) { uint16_t u; memcpy(&u, p, 2); return u; }
static uint32_t U4(uint8_t* p) { uint32_t u; memcpy(&u, p, 4); return u; }
static int32_t  I4(uint8_t* p) { int32_t  i; memcpy(&i, p, 4); return i; }
static float    R4(uint8_t* p) { float    r; memcpy(&r, p, 4); return r; }
static double   R8(uint8_t* p) { double   r; memcpy(&r, p, 8); return r; }
// 符号扩展：将指定位宽的无符号值转换为带符号
static int32_t exsign(uint32_t v, int bits) {
    return (int32_t)(v & (1 << (bits - 1)) ? v | (~0u << bits) : v);
}
// CRC32 校验（用于验证OEM4 消息完整性）
static uint32_t rtk_crc32(const uint8_t* buff, int len) {
    uint32_t crc = 0;
    for (int i = 0; i < len; i++) {
        crc ^= buff[i];
        for (int j = 0; j < 8; j++) {
            if (crc & 1) crc = (crc >> 1) ^ 0xEDB88320;
            else crc >>= 1;
        }
    }
    return crc;
}

// 在导航数据表中查找指定卫星的星历。
eph_t* find_eph(nav_t* nav, int sat)
{
    for (int i = 0; i < MAXSAT; i++) {
        if (nav->eph[i].sat == sat) {
            return &nav->eph[i];
        }
    }
    return nullptr;
}

// 星历是否具备参与卫星位置计算的基本参数。
bool has_valid_eph(const eph_t* eph)
{
    if (!eph) return false;
    if (eph->sat <= 0) return false;
    if (eph->sqrtA <= 0.0) return false;
    return true;
}
//根据卫星号 sat 和观测信号类型 code，返回该观测值对应的载波频率
// 对于 GLONASS 需要使用卫星频率通道（glo_fcn）进行偏移计算
static double sat2freq(int sat, uint8_t code, nav_t* nav)
{
    int sys, prn = 0;
    sys = satsys(sat, &prn);

    if (sys == SYS_GPS) {

        if (code == CODE_L1C || code == CODE_L1P)
            return FREQ_GPS_L1;

        if (code == CODE_L2P || code == CODE_L2C || code == CODE_L2W)
            return FREQ_GPS_L2;

        if (code == CODE_L5Q)
            return 1176.45e6;
    }

    else if (sys == SYS_GLO) {

        int fcn = nav->glo_fcn[prn - 1];
        if (fcn < -7 || fcn > 6) fcn = 0;

        if (code == CODE_L1C)
            return 1602.0e6 + fcn * 562.5e3;

        if (code == CODE_L2C)
            return 1246.0e6 + fcn * 437.5e3;
    }

    else if (sys == SYS_GAL) {

        if (code == CODE_L1C)
            return 1575.42e6;

        if (code == CODE_L5Q)
            return 1176.45e6;
    }

    else if (sys == SYS_CMP) {

        if (code == CODE_L2I)   // B1I
            return FREQ_BDS_B1;

        if (code == CODE_L6I)   // B3I
            return FREQ_BDS_B3;
    }

    return 0.0;
}

// 将观测码映射到观测数组索引（L1->0, L2->1, L5->2 等）
static int code2idx(int sys, int code)
{
    if (sys == SYS_GPS) {

        if (code == CODE_L1C || code == CODE_L1P)
            return 0;

        if (code == CODE_L2P || code == CODE_L2C || code == CODE_L2W)
            return 1;

        if (code == CODE_L5Q)
            return 2;
    }

    else if (sys == SYS_GLO) {

        if (code == CODE_L1C)
            return 0;

        if (code == CODE_L2C)
            return 1;
    }

    else if (sys == SYS_GAL) {

        if (code == CODE_L1C)
            return 0;

        if (code == CODE_L5Q)
            return 2;
    }

    else if (sys == SYS_CMP) {

        if (code == CODE_L2I)   // B1I
            return 0;

        if (code == CODE_L6I)   // B3I
            return 1;
    }

    return -1;
}

// 将 NovAtel 报文中的信号类型编号转换为本程序使用的观测码
static int sig2code(int sys, int sigtype) {
    if (sys == SYS_GPS) {
        switch (sigtype) {
        case 0:  return CODE_L1C;
        case 9:  return CODE_L2P;
        }
    }
    else if (sys == SYS_GLO) {
        switch (sigtype) {
        case 0:  return CODE_L1C;
        case 1:  return CODE_L2C;
        case 5:  return CODE_L2P;
        }
    }
    else if (sys == SYS_GAL) {
        switch (sigtype) {
        case 2:  return CODE_L1C;
        case 12: return CODE_L5Q;
        }
    }
    else if (sys == SYS_CMP) {
        switch (sigtype) {
        case 0:  return CODE_L2I;  // B1I D1
        case 4:  return CODE_L2I;  // B1I D2

        case 2:  return CODE_L6I;  // B3I D1
        case 6:  return CODE_L6I;  // B3I D2
        }
    }
    return 0;
}
// 32 位状态字stat中解析出卫星系统、信号类型、跟踪状态、锁定标志、奇偶校验、相位半周等信息，返回信号观测索引（0/1/2）或 -1
static int decode_track_stat(uint32_t stat, int* sys, int* code, int* track,int* plock, int* clock, int* parity, int* halfc) 
{
    int satsys, sigtype, idx = -1;

    *track = stat & 0x1F; //取最低 5 位，得到跟踪状态
    *plock = (stat >> 10) & 1;//取第 10 位，表示伪距是否锁定
    *parity = (stat >> 11) & 1;//取第 11 位，表示奇偶校验是否通过
    *clock = (stat >> 12) & 1;//取第 12 位，表示载波锁定状态
    satsys = (stat >> 16) & 7;//取第 16~18 位，得到卫星系统编号
    *halfc = (stat >> 28) & 1; //取第 28 位，表示半周模糊度/半周标志
    sigtype = (stat >> 21) & 0x1F;//取第 21~25 位，得到信号类型编号
    //判断系统类型
    switch (satsys) {
    case 0: *sys = SYS_GPS; break;
    case 1: *sys = SYS_GLO; break;
    case 2: *sys = SYS_SBS; break;
    case 3: *sys = SYS_GAL; break;
    case 4: *sys = SYS_CMP; break;
    case 5: *sys = SYS_QZS; break;
    case 6: *sys = SYS_IRN; break;
    default: return -1;
    }
    //信号类型编号转换为程序内部CODE&&进一步转换为观测数组下标
    if (!(*code = sig2code(*sys, sigtype)) ||
        (idx = code2idx(*sys, *code)) < 0) {
        return -1;
    }
    return idx;
}

// obs 集合中查找或创建对应卫星的观测项，返回索引
static int obsindex(obs_t* obs, gtime_t time, int sat) {
    int i, j;
    if (obs->n >= MAXOBS) return -1;
    //查找到返回下标
    for (i = 0; i < obs->n; i++) {
        if (obs->data[i].sat == sat) return i;
    }
    obs->data[i].time = time;
    obs->data[i].sat = sat;
    for (j = 0; j < NFREQ + NEXOBS; j++) {
        obs->data[i].L[j] = obs->data[i].P[j] = 0.0;
        obs->data[i].D[j] = 0.0;
        obs->data[i].SNR[j] = obs->data[i].LLI[j] = 0;
        obs->data[i].code[j] = CODE_NONE;
    }
    obs->n++;
    return i;
}

// ------------------ 解码函数 ------------------
// 解码 RANGEB 消息（NovAtel Raw Range Binary）(原始观测值报文)
//可以读取伪距 P，载波相位 L，多普勒 D，信噪比 SNR，锁定时间 lockt，跟踪状态字 track status
static int decode_rangeb(raw_t* raw) {
    uint8_t* p = raw->buff + OEM4HLEN;
    int i, index, nobs, prn, sat, sys, code, idx, track, plock, clock, parity, halfc, lli;
    double psr, adr, dop, snr, lockt, tt, freq, glo_bias = 0.0;

    // 报文中第一个字段是卫星数（32位无符号）
    nobs = U4(p);
    if (raw->len < OEM4HLEN + 4 + nobs * 44) {
        cerr << "RANGEB length error" << endl;
        return -1;
    }

    // 如果当前 raw 中缓存的观测时间与新消息时间不同，则重置观测集合
    raw->reset_obs();

    for (i = 0, p += 4; i < nobs; i++, p += 44) {
        //跟踪字中读取系统 sys，观测码 code，频点索引 idx，锁定标志 plock，伪距有效标志 clock，奇偶校验 parity，半周标志 halfc
        if ((idx = decode_track_stat(U4(p + 40), &sys, &code, &track, &plock, &clock, &parity, &halfc)) < 0) {
            continue;
        }
        prn = U2(p);
        if (sys == SYS_GLO) prn -= 37;
        if (!(sat = satno(sys, prn))) continue;

        // GLONASS 奇偶校验未知则跳过（无法可靠解码）
        if (sys == SYS_GLO && !parity) continue;

        psr = R8(p + 4);//伪距 m
        adr = R8(p + 16);//载波相位 周
        dop = R4(p + 28);//多普勒 hz
        snr = R4(p + 32);//信噪比
        lockt = R4(p + 36);//锁定时间，判断周跳

        // 周跳检测：比较 lock time 与上次记录，判断是否发生周跳
        if (raw->tobs[sat - 1][idx].week != 0) {
            //计算当前观测时间和上一条观测时间的时间差
            tt = timediff(raw->time, raw->tobs[sat - 1][idx]);
            lli = (lockt - raw->lockt[sat - 1][idx] + 0.05 <= tt) ? LLI_SLIP : 0;
        }
        else {
            lli = 0;
        }
        if (!parity) lli |= LLI_HALFC;
        if (halfc)   lli |= LLI_HALFA;
        //更新当前卫星当前频点的上一次历元状态
        raw->tobs[sat - 1][idx] = raw->time;
        raw->lockt[sat - 1][idx] = lockt;
        raw->halfc[sat - 1][idx] = halfc;

        // 相位符号调整
        adr = -adr;

        if (!clock) psr = 0.0;
        if (!plock) adr = dop = 0.0;

        // 将解析出的观测值写入 obs 结构
        if ((index = obsindex(&raw->obs, raw->time, sat)) >= 0) {
            raw->obs.data[index].L[idx] = adr;
            raw->obs.data[index].P[idx] = psr;
            raw->obs.data[index].D[idx] = (float)dop;
            raw->obs.data[index].SNR[idx] = snr / SNR_UNIT;
            raw->obs.data[index].LLI[idx] = (unsigned char)lli;
            raw->obs.data[index].code[idx] = (unsigned char)code;
        }
    }
    return 1;
}

// 解码 RANGECMPB 消息(压缩格式，每颗卫星24字节)
static int decode_rangecmpb(raw_t* raw) {
    uint8_t* p = raw->buff + OEM4HLEN;
    int i, nobs, prn, sat, sys, code, idx, track, plock, clock, parity, halfc, lli, index;
    double psr, adr, adr_rolls, lockt, tt, dop, snr, freq, glo_bias = 0.0;

    // 读取观测数并进行长度验证
    nobs = U4(p);
    if (raw->len < OEM4HLEN + 4 + nobs * 24) {
        cerr << "RANGECMPB length error" << endl;
        return -1;
    }

    raw->reset_obs();

    for (i = 0, p += 4; i < nobs; i++, p += 24) {
        if ((idx = decode_track_stat(U4(p), &sys, &code, &track, &plock, &clock, &parity, &halfc)) < 0) {
            continue;
        }
        prn = U1(p + 17);
        if (sys == SYS_GLO) prn -= 37;
        if (!(sat = satno(sys, prn))) continue;
        if (sys == SYS_GLO && !parity) continue;

        dop = exsign(U4(p + 4) & 0xFFFFFFF, 28) / 256.0;//多普勒恢复
        psr = (U4(p + 7) >> 4) / 128.0 + U1(p + 11) * 2097152.0;//伪距恢复

        // 如果可以得到频率，用相位计数和周跳修正恢复完整相位
        if ((freq = sat2freq(sat, code, &raw->nav)) != 0.0) {
            adr = I4(p + 12) / 256.0;
            adr_rolls = (psr * freq / Clight + adr) / MAXVAL;
            adr = -adr + MAXVAL * floor(adr_rolls + (adr_rolls <= 0 ? -0.5 : 0.5));
        }
        else {
            adr = 1e-9;
        }
        lockt = (U4(p + 18) & 0x1FFFFF) / 32.0;

        if (raw->tobs[sat - 1][idx].week != 0) {
            tt = timediff(raw->time, raw->tobs[sat - 1][idx]);
            lli = (lockt - raw->lockt[sat - 1][idx] + 0.05 <= tt) ? LLI_SLIP : 0;
        }
        else {
            lli = 0;
        }
        if (!parity) lli |= LLI_HALFC;
        if (halfc)   lli |= LLI_HALFA;

        raw->tobs[sat - 1][idx] = raw->time;
        raw->lockt[sat - 1][idx] = lockt;
        raw->halfc[sat - 1][idx] = halfc;

        snr = ((U2(p + 20) & 0x3FF) >> 5) + 20.0;
        if (!clock) psr = 0.0;
        if (!plock) adr = dop = 0.0;

        if ((index = obsindex(&raw->obs, raw->time, sat)) >= 0) {
            raw->obs.data[index].L[idx] = adr;
            raw->obs.data[index].P[idx] = psr;
            raw->obs.data[index].D[idx] = (float)dop;
            raw->obs.data[index].SNR[idx] = snr / SNR_UNIT;
            raw->obs.data[index].LLI[idx] = (unsigned char)lli;
            raw->obs.data[index].code[idx] = (unsigned char)code;
        }
    }
    return 1;
}
// 解码 GPS 广播星历报文，并写入 nav.eph[sat-1]。
static int decode_gpsephemb(raw_t* raw)
{
    uint8_t* p = raw->buff + OEM4HLEN;

    int prn = (int)U4(p + 0);      //  GPS 卫星 PRN
    int week = (int)U4(p + 24);     //  toe 对应的 GPS 周
    double toe = R8(p + 32);        // 星历参考时刻 toe，单位秒
    double toc = R8(p + 164);       // 卫星钟参考时刻toc
    double A = R8(p + 40);        // 轨道长半轴 A
    double deln = R8(p + 48);        // 平均角速度改正 Δn
    double M0 = R8(p + 56);        // 平近点角 M0
    double ecc = R8(p + 64);        // 偏心率 e
    double omg = R8(p + 72);        // 近地点角距 ω
    double cuc = R8(p + 80);        // 升交距角改正
    double cus = R8(p + 88);        
    double crc = R8(p + 96);        //轨道半径改正
    double crs = R8(p + 104);       
    double cic = R8(p + 112);       //轨道倾角改正
    double cis = R8(p + 120);       
    double i0 = R8(p + 128);       //轨道倾角i0
    double idot = R8(p + 136);     // 倾角变化率IDOT
    double OMG0 = R8(p + 144);       // 升交点经度Ω0
    double OMGd = R8(p + 152);       // 升交点赤经变化率Ωdot
    
    double af0 = R8(p + 180);       //钟差二次多项式参数 af0/af1/af2
    double af1 = R8(p + 188);   
    double af2 = R8(p + 196);      

    int sat = satno(SYS_GPS, prn);
    if (!sat) return 0;

    eph_t* eph = &raw->nav.eph[sat - 1];
    memset(eph, 0, sizeof(eph_t));//清零旧星历
    //存储当前星历
    eph->sat = sat;
    eph->prn = prn;
    eph->week = week;

    eph->toe = gpst2time(week, toe);
    eph->toc = gpst2time(week, toc);

    eph->sqrtA = sqrt(A);
    eph->A = A;
    eph->e = ecc;
    eph->M0 = M0;
    eph->OMG0 = OMG0;
    eph->omg = omg;
    eph->i0 = i0;

    eph->deln = deln;
    eph->OMGd = OMGd;
    eph->idot = idot;

    eph->cuc = cuc;
    eph->cus = cus;
    eph->crc = crc;
    eph->crs = crs;
    eph->cic = cic;
    eph->cis = cis;

    eph->af0 = af0;
    eph->af1 = af1;
    eph->af2 = af2;

    raw->ephsat = sat;//最新更新的卫星编号
    return 2;
}
// 解码 BDS 广播星历报文。注意原始 toe/toc 属于 BDT。
static int decode_bdsephemerisb(raw_t* raw)
{
    uint8_t* p = raw->buff + OEM4HLEN;

    int prn = (int)U4(p + 0);      
    int week = (int)U4(p + 4);      
    double tgd1 = R8(p + 20);        //BDS 的两个群延迟参数 TGD1、TGD2 s
    double tgd2 = R8(p + 28);        
    double toc = (double)U4(p + 40);  // 钟差参考时刻 toc BDT周内秒
    double af0 = R8(p + 44);          //BDS 卫星钟差参数
    double af1 = R8(p + 52);         
    double af2 = R8(p + 60);         
    double toe = (double)U4(p + 72);  //轨道参考时刻 toe BDT周内秒
    double sqrtA = R8(p + 76);       
    double ecc = R8(p + 84);       // 偏心率
    double omg = R8(p + 92);       
    double deln = R8(p + 100);      
    double M0 = R8(p + 108);      
    double OMG0 = R8(p + 116);      
    double OMGd = R8(p + 124);      
    double i0 = R8(p + 132);    
    double idot = R8(p + 140);      
    double cuc = R8(p + 148);     
    double cus = R8(p + 156);      
    double crc = R8(p + 164);     
    double crs = R8(p + 172);      
    double cic = R8(p + 180);    
    double cis = R8(p + 188);      

    int sat = satno(SYS_CMP, prn);
    if (!sat) return 0;

    eph_t* eph = &raw->nav.eph[sat - 1];
    memset(eph, 0, sizeof(eph_t));//清零旧星历

    eph->sat = sat;
    eph->prn = prn;
    eph->week = week;
    // BDS广播星历中的toe/toc属于BDT。toes/tocs保留原始BDT周内秒供BDS轨道升交点经度项使用；toe/toc转为GPST，便于和观测历元GPST求时间差。
    eph->toes = toe;  
    eph->tocs = toc;
    eph->toe = bdt2gpst(bdt2time(week, toe));
    eph->toc = bdt2gpst(bdt2time(week, toc));
    eph->sqrtA = sqrtA;
    eph->A = sqrtA * sqrtA;
    eph->e = ecc;
    eph->M0 = M0;
    eph->OMG0 = OMG0;
    eph->omg = omg;
    eph->i0 = i0;

    eph->deln = deln;
    eph->OMGd = OMGd;
    eph->idot = idot;

    eph->cuc = cuc;
    eph->cus = cus;
    eph->crc = crc;
    eph->crs = crs;
    eph->cic = cic;
    eph->cis = cis;

    eph->af0 = af0;
    eph->af1 = af1;
    eph->af2 = af2;
    eph->tgd[0] = tgd1;
    eph->tgd[1] = tgd2;
    raw->ephsat = sat;
    return 2;
}
// 解码 OEM4 消息头并分发到对应的子解码器
int decode_oem4(raw_t* raw) {
    int type = U2(raw->buff + 4); //读取消息 ID
    //定义 GPS 周、消息格式、时间状态、周内秒
    int week, msg, stat; 
    double tow;

    week = U2(raw->buff + 14);
    msg = (U1(raw->buff + 6) >> 4) & 0x3;
    stat = U1(raw->buff + 13);
    tow = U4(raw->buff + 16) * 0.001;
    //计算当前报文的 CRC32，并和报文尾部保存的 CRC 比较
    if (rtk_crc32(raw->buff, raw->len) != U4(raw->buff + raw->len)) {
        cout << "CRC fail type=" << type
            << " week=" << week
            << " tow=" << tow << endl;
        return -1;
    }

    if (stat == 20 || week == 0) {
        cout << "time/stat skip type=" << type
            << " week=" << week
            << " tow=" << tow
            << " stat=" << stat << endl;
        return 0;
    }

    week = adjgpsweek(week);
    raw->time = gpst2time(week, tow);

    if (msg != 0) {
        cout << "msg skip type=" << type
            << " week=" << week
            << " tow=" << tow
            << " msg=" << msg << endl;
        return 0;
    }

    switch (type) {
    case ID_RANGE:          return decode_rangeb(raw);
    case ID_RANGECMP:       return decode_rangecmpb(raw);
    case ID_GPSEPHEM:       return decode_gpsephemb(raw);
    case ID_BDSEPHEMERIS:   return decode_bdsephemerisb(raw);
    default: return 0;
    }
}

// 同步头检测：用于在字节流中查找OEM4 报文起始标志
int sync_oem4(uint8_t* buff, uint8_t data) {
    buff[0] = buff[1]; buff[1] = buff[2]; buff[2] = data;
    return buff[0] == OEM4SYNC1 && buff[1] == OEM4SYNC2 && buff[2] == OEM4SYNC3;
}
// 文件输入：先搜索同步头，再按报文长度读完一帧 OEM4 数据。
int input_oem4f(raw_t* raw, FILE* fp) {
    int i, data; //data 保存当前读到的一个字节

    if (raw->nbyte == 0) {
        for (i = 0;; i++) {
            if ((data = fgetc(fp)) == EOF) return -2;
            if (sync_oem4(raw->buff, (uint8_t)data)) break;
            if (i >= 4096) return 0;
        }
    }
    if (fread(raw->buff + 3, 7, 1, fp) < 1) return -2;
    raw->nbyte = 10;

    raw->len = U2(raw->buff + 8) + OEM4HLEN;
    if (raw->len > MAXRAWLEN - 4) {
        raw->nbyte = 0;
        return -1;
    }
    if (fread(raw->buff + 10, raw->len - 6, 1, fp) < 1) return -2;
    raw->nbyte = 0;

    return decode_oem4(raw);
    //2：文件结束或读取不完整
    //1：报文错误，比如 CRC 错、长度错
    //0：读到报文，但不是当前程序需要的消息，或者时间无效
    // 1：成功解码观测值
    // 2：成功解码星历
}
