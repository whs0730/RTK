#include "obs.h"
//根据卫星系统和prn计算卫星号
int satno(int sys, int prn)
{
    if (prn <= 0) return 0;
    switch (sys) {
    case SYS_GPS:
        if (prn < MINPRNGPS || MAXPRNGPS < prn) return 0;
        return prn - MINPRNGPS + 1;
    case SYS_GLO:
        if (prn < MINPRNGLO || MAXPRNGLO < prn) return 0;
        return NSATGPS + prn - MINPRNGLO + 1;
    case SYS_GAL:
        if (prn < MINPRNGAL || MAXPRNGAL < prn) return 0;
        return NSATGPS + NSATGLO + prn - MINPRNGAL + 1;
    case SYS_QZS:
        if (prn < MINPRNQZS || MAXPRNQZS < prn) return 0;
        return NSATGPS + NSATGLO + NSATGAL + prn - MINPRNQZS + 1;
    case SYS_CMP:
        if (prn < MINPRNCMP || MAXPRNCMP < prn) return 0;
        return NSATGPS + NSATGLO + NSATGAL + NSATQZS + prn - MINPRNCMP + 1;
    case SYS_IRN:
        if (prn < MINPRNIRN || MAXPRNIRN < prn) return 0;
        return NSATGPS + NSATGLO + NSATGAL + NSATQZS + NSATCMP + prn - MINPRNIRN + 1;
    case SYS_LEO:
        if (prn < MINPRNLEO || MAXPRNLEO < prn) return 0;
        return NSATGPS + NSATGLO + NSATGAL + NSATQZS + NSATCMP + NSATIRN +
            prn - MINPRNLEO + 1;
    case SYS_SBS:
        if (prn < MINPRNSBS || MAXPRNSBS < prn) return 0;
        return NSATGPS + NSATGLO + NSATGAL + NSATQZS + NSATCMP + NSATIRN + NSATLEO +
            prn - MINPRNSBS + 1;
    }
    return 0;
}
// 根据内部 sat 编号反推出卫星系统和系统内 PRN。
int satsys(int sat, int* prn)
{
    if (sat <= 0) return SYS_NONE;

    if (sat <= NSATGPS) {
        if (prn) *prn = sat;
        return SYS_GPS;
    }
    else if (sat <= NSATGPS + NSATGLO) {
        if (prn) *prn = sat - NSATGPS;
        return SYS_GLO;
    }
    else if (sat <= NSATGPS + NSATGLO + NSATGAL) {
        if (prn) *prn = sat - NSATGPS - NSATGLO;
        return SYS_GAL;
    }
    else if (sat <= NSATGPS + NSATGLO + NSATGAL + NSATQZS) {
        if (prn) *prn = sat - NSATGPS - NSATGLO - NSATGAL;
        return SYS_QZS;
    }
    else if (sat <= NSATGPS + NSATGLO + NSATGAL + NSATQZS + NSATCMP) {
        if (prn) *prn = sat - NSATGPS - NSATGLO - NSATGAL - NSATQZS;
        return SYS_CMP;
    }
    else if (sat <= NSATGPS + NSATGLO + NSATGAL + NSATQZS + NSATCMP + NSATIRN) {
        if (prn) *prn = sat - NSATGPS - NSATGLO - NSATGAL - NSATQZS - NSATCMP;
        return SYS_IRN;
    }
    return SYS_NONE;
}
// 将内部卫星编号格式化为常见卫星号，例如 G25、C13。
string sat2id(int sat)
{
    int prn = 0;
    int sys = satsys(sat, &prn);
    char buff[16] = { 0 };

    switch (sys) {
    case SYS_GPS:
        sprintf_s(buff, "G%02d", prn);
        break;
    case SYS_GLO:
        sprintf_s(buff, "R%02d", prn);
        break;
    case SYS_GAL:
        sprintf_s(buff, "E%02d", prn);
        break;
    case SYS_QZS:
        sprintf_s(buff, "J%02d", prn);
        break;
    case SYS_CMP:
        sprintf_s(buff, "C%02d", prn);
        break;
    case SYS_IRN:
        sprintf_s(buff, "I%02d", prn);
        break;
    default:
        sprintf_s(buff, "%03d", sat);
        break;
    }
    return string(buff);
}
// 判断卫星位置是否有效。
bool IsValidSatpos(satpos_t* sat)
{
    if (sat == nullptr) {
        return false;
    }

    if (sat->pos[0] == 0.0 && sat->pos[1] == 0.0 && sat->pos[2] == 0.0) {
        return false;
    }

    return true;
}

