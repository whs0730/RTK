#pragma once
#include "obs.h"
#include "spp.h"
#include "satpos.h"
#include "decode.h"
#include "error_correction.h"
//星历最大允许年龄
int CaculateSatellitePositions(obsd_t* obs, int n, nav_t* nav, satpos_t* sats);