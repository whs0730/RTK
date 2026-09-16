#pragma once

// ====================== NovAtel OEM 报文与观测常量 ======================
#define OEM4SYNC1       0xAA        // OEM4/OEM6 二进制报文同步头第 1 字节
#define OEM4SYNC2       0x44        // OEM4/OEM6 二进制报文同步头第 2 字节
#define OEM4SYNC3       0x12        // OEM4/OEM6 二进制报文同步头第 3 字节
#define OEM3SYNC1       0xAA
#define OEM3SYNC2       0x44
#define OEM3SYNC3       0x11
#define OEM4HLEN        28          // OEM4/OEM6 报文头长度，解码时用于定位数据体
#define OEM3HLEN        12
#define MAXRAWLEN       8192        // 单次原始报文/网口接收缓冲区最大字节数
#define MAXOBS          512         // 单历元最多保存的观测值数量
#define MAXSAT          220         // 内部可编号的最大卫星数量
#define NFREQ           3           // 基本频点数，P/L/D/SNR/code 数组的主频点部分
#define NEXOBS          3           // 扩展观测槽数量，和 NFREQ 一起决定每星观测数组长度
#define CODE_NONE       0           // 未识别或未赋值的观测码
#define LLI_SLIP        1           // 周跳标记
#define LLI_HALFC       2           // 半周模糊标记：奇偶校验异常
#define LLI_HALFA       4           // 半周模糊标记：跟踪状态中的 half-cycle 标志
#define SNR_UNIT        1.0         // SNR 存储比例，当前直接按 dB-Hz 保存
#define WL1             0.1902936727984
#define WL2             0.2442102134246
#define Clight          299792458.0 // 真空光速，伪距、钟差和多普勒换算使用
#define MAXVAL          8388608.0   // RANGECMP 载波相位压缩值的整周展开尺度
#define OFF_FRQNO       -7

// ====================== 卫星系统标识 ======================
#define SYS_NONE    0x00            // 无效或未知卫星系统
#define SYS_GPS     0x01            // GPS
#define SYS_SBS     0x02            // SBAS
#define SYS_GLO     0x04            // GLONASS
#define SYS_GAL     0x08            // Galileo
#define SYS_QZS     0x10            // QZSS
#define SYS_CMP     0x20            // BDS/BeiDou
#define SYS_IRN     0x40            // NavIC/IRNSS
#define SYS_LEO     0x80            // LEO 扩展系统
#define SYS_ALL     0xFF

// ====================== 观测码类型 ======================
#define CODE_L1C        1           // GPS/GLO/GAL L1 C/A 或 E1 类信号
#define CODE_L1P        2           // GPS L1 P(Y)
#define CODE_L1W        3
#define CODE_L1Y        4
#define CODE_L1M        5
#define CODE_L1N        6
#define CODE_L1S        7
#define CODE_L1L        8           // GPS L1C pilot,当前只做解码映射
#define CODE_L1X        9
#define CODE_L2I        40          // BDS B1I 信号，本工程映射到观测数组第 1 频点
#define CODE_L2C        10          // GPS/GLO L2C
#define CODE_L2D        11
#define CODE_L2S        12          // GPS L2C pilot，当前只做解码映射
#define CODE_L2L        13
#define CODE_L2X        14
#define CODE_L2P        15          // GPS/GLO L2 P(Y)
#define CODE_L2W        16          // GPS L2 P(Y) 加密跟踪码
#define CODE_L2Y        17
#define CODE_L2M        18
#define CODE_L2N        19
#define CODE_L5I        20
#define CODE_L5Q        21          // GPS/Galileo L5/E5a Q，当前只做解码映射
#define CODE_L5X        22
#define CODE_L7I        23
#define CODE_L7Q        24
#define CODE_L7X        25
#define CODE_L6A        26
#define CODE_L6B        27
#define CODE_L6C        28
#define CODE_L6X        29
#define CODE_L6Z        30
#define CODE_L6S        31
#define CODE_L6L        32
#define CODE_L6I        60          // BDS B3I 信号，本工程映射到观测数组第 2 频点
#define CODE_L8I        33
#define CODE_L8Q        34
#define CODE_L8X        35

// ====================== NovAtel 消息 ID ======================
#define ID_RANGECMP     140         // RANGECMPB 压缩观测值，实时流主要观测来源
#define ID_RANGE        43          // RANGEB 普通观测值
#define ID_RAWEPHEM     41
#define ID_IONUTC       8
#define ID_GLOEPHEMERIS 723
#define ID_GALEPHEMERIS 1122
#define ID_GPSEPHEM     7           // GPS 广播星历
#define ID_BDSEPHEMERIS 1696        // BDS 广播星历
#define ID_QZSSRAWEPHEM 1331
#define ID_NAVICEPHEMERIS 2123
#define ID_RGEB         15
#define ID_RGED         65
#define ID_REPB         14
#define ID_FRMB         54
#define ID_IONB         16
#define ID_UTCB         17

// ====================== 各系统 PRN 范围与内部编号偏移 ======================
#define MINPRNGPS   1               // GPS 最小 PRN
#define MAXPRNGPS   32              // GPS 最大 PRN
#define NSATGPS     (MAXPRNGPS-MINPRNGPS+1) // GPS 内部编号数量
#define MINPRNGLO   1               // GLONASS 最小槽号/PRN
#define MAXPRNGLO   27              // GLONASS 最大槽号/PRN
#define NSATGLO     (MAXPRNGLO-MINPRNGLO+1) // GLONASS 内部编号数量
#define MINPRNGAL   1               // Galileo 最小 PRN
#define MAXPRNGAL   36              // Galileo 最大 PRN
#define NSATGAL     (MAXPRNGAL-MINPRNGAL+1) // Galileo 内部编号数量
#define MINPRNQZS   193             // QZSS 最小 PRN
#define MAXPRNQZS   202             // QZSS 最大 PRN
#define MINPRNQZS_S 183
#define MAXPRNQZS_S 191
#define NSATQZS     (MAXPRNQZS-MINPRNQZS+1) // QZSS 内部编号数量
#define MINPRNCMP   1               // BDS 最小 PRN
#define MAXPRNCMP   63              // BDS 最大 PRN，覆盖 BDS-2/BDS-3 常见编号
#define NSATCMP     (MAXPRNCMP-MINPRNCMP+1) // BDS 内部编号数量
#define MINPRNIRN   1               // NavIC/IRNSS 最小 PRN
#define MAXPRNIRN   14              // NavIC/IRNSS 最大 PRN
#define NSATIRN     (MAXPRNIRN-MINPRNIRN+1) // NavIC/IRNSS 内部编号数量
#define MINPRNLEO   1               // LEO 扩展系统最小编号
#define MAXPRNLEO   10              // LEO 扩展系统最大编号
#define NSATLEO     (MAXPRNLEO-MINPRNLEO+1) // LEO 内部编号数量
#define MINPRNSBS   120             // SBAS 最小 PRN
#define MAXPRNSBS   158             // SBAS 最大 PRN
#define NSATSBS     (MAXPRNSBS-MINPRNSBS+1) // SBAS 内部编号数量

// ====================== 常用频率与地球自转角速度 ======================
#define FREQ_GPS_L1 1575.42e6       // GPS L1 频率，单位 Hz
#define FREQ_GPS_L2 1227.60e6       // GPS L2 频率，单位 Hz

#define FREQ_BDS_B1 1561.098e6      // BDS B1I 频率，单位 Hz
#define FREQ_BDS_B2 1207.14e6
#define FREQ_BDS_B3 1268.52e6       // BDS B3I 频率，单位 Hz

#define OMGed_GPS 7.2921151467e-5   // GPS 坐标计算使用的地球自转角速度，rad/s
#define OMGed_BDS 7.2921150e-5      // BDS 坐标计算使用的地球自转角速度，rad/s
