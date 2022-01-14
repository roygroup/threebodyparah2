#ifndef REGION_H
#define REGION_H

#include <map>
#include <vector>

#include "utils.h"
#include "jpoint.h"

// constants that determine the boundaries of $ (R, s, \varphi) coordinate space

namespace regions
{

constexpr double RMIN      = 2.20;
constexpr double R_SEXTRAP = 3.25;
constexpr double RMAX      = 6.25;

constexpr double SMIN      = 1.00;
constexpr double SMAX      = 3.85;

constexpr double UCHANGE   = 0.00001/HALFPI;

constexpr double S_EX_0 = 3.55;
constexpr double S_EX_1 = 3.85;
constexpr double R_EX_0 = 2.20;
constexpr double R_EX_1 = 2.25;

enum class Region {
    AxilrodTellerMuto,
    R_extrapolation,
    s_extrapolation,
    R_and_s_extrapolation,
    Bounded,
};

std::string region_name(const Region&);
unsigned int number_of_raw_JacobiPoints(const Region&);
Region which_region(const JacobiPoint&);
std::vector<JacobiPoint> get_raw_jpoints(const JacobiPoint&, const Region&);

}

#endif
