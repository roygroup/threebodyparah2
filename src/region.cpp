#include <iostream>

#include "utils.h"
#include "region.h"
#include "trapezoid.h"

namespace regions
{

// UNUSED
std::string region_name(const Region& reg) {
    if (reg == Region::AxilrodTellerMuto) {
        return "AxilrodTellerMuto    ";
    }
    else if (reg == Region::Bounded) {
        return "Bounded              ";
    }
    else if (reg == Region::s_extrapolation) {
        return "s_extrapolation      ";
    }
    else if (reg == Region::R_extrapolation) {
        return "R_extrapolation      ";
    }
    else if (reg == Region::R_and_s_extrapolation) {
        return "R_and_s_extrapolation";
    }
    else {
        std::cerr << "ERROR in " << __func__;
        std::cerr << ": unreachable; the JacobiPoint should have fallen into";
        std::cerr << " one of the previous cases";
        std::cerr << std::endl;
        exit(EXIT_FAILURE);
    }
}

unsigned int number_of_raw_JacobiPoints(const Region& reg) {
    unsigned int N_raw_JacobiPoints;
    switch (reg)
    {
    case Region::AxilrodTellerMuto:
        N_raw_JacobiPoints = 0;
        break;
    case Region::R_extrapolation:
        N_raw_JacobiPoints = 3;
        break;
    case Region::s_extrapolation:
        N_raw_JacobiPoints = 2;
        break;
    case Region::R_and_s_extrapolation:
        N_raw_JacobiPoints = 6;
        break;
    case Region::Bounded:
        N_raw_JacobiPoints = 1;
        break;
    }

    return N_raw_JacobiPoints;
}

Region which_region(const JacobiPoint& jpoint) {
    double R = jpoint.get_R();
    double s = jpoint.get_s();
    double u = jpoint.get_u();

    Region reg;
    if        (R >= R_SEXTRAP && s >  SMAX) {
        reg = Region::AxilrodTellerMuto;
    } else if (R >= RMAX      && s >= SMIN) {
        reg = Region::AxilrodTellerMuto;
    } else if (R <  RMIN      && s >= SMIN && s <= SMAX) {
        reg = Region::R_extrapolation;
    } else if (R <  RMIN      && s > SMAX) {
        reg = Region::R_and_s_extrapolation;
    } else if (R <  R_SEXTRAP && R >= RMIN && s > SMAX) {
        if (u >= UCHANGE) {
            reg = Region::s_extrapolation;
        } else {
            reg = Region::Bounded;
        }
    } else if (RMIN <= R && R < RMAX && SMIN <= s && s <= SMAX) {
        reg = Region::Bounded;
    } else {
        std::cerr << "ERROR in " << __func__;
        std::cerr << ": unreachable; the JacobiPoint should have fallen into";
        std::cerr << " one of the previous cases";
        std::cerr << std::endl;
        exit(EXIT_FAILURE);
    }

    return reg;
}

// --- RAW JACOBIPOINT ------------------------------------------------------------

std::vector<JacobiPoint> get_raw_jpoints(const JacobiPoint& jpoint, const Region& reg) {
    std::vector<JacobiPoint> raw_jpoints;
    raw_jpoints.reserve(regions::number_of_raw_JacobiPoints(reg));

    if (reg == Region::AxilrodTellerMuto) {
        // no RKHS calculations involved, do nothing
    }
    else if (reg == Region::Bounded) {
        raw_jpoints.push_back(jpoint);
    }
    else if (reg == Region::s_extrapolation) {
        double R, s, u;
        std::tie(R, s, u) = jpoint.unpack();

        raw_jpoints.push_back( JacobiPoint::from_jacobicoords(R, regions::S_EX_0, u) );
        raw_jpoints.push_back( JacobiPoint::from_jacobicoords(R, regions::S_EX_1, u) );
    }
    else if (reg == Region::R_extrapolation) {
        double R, s, u;
        std::tie(R, s, u) = jpoint.unpack();

        raw_jpoints.push_back( JacobiPoint::from_jacobicoords(regions::R_EX_0, s, u) );
        raw_jpoints.push_back( JacobiPoint::from_jacobicoords(regions::R_EX_1, s, u) );
        raw_jpoints.push_back( JacobiPoint::from_jacobicoords(R              , s, u) );
    }
    else if (reg == Region::R_and_s_extrapolation) {
        double R, s, u;
        std::tie(R, s, u) = jpoint.unpack();

        raw_jpoints.push_back( JacobiPoint::from_jacobicoords(regions::R_EX_0, regions::S_EX_0, u) );
        raw_jpoints.push_back( JacobiPoint::from_jacobicoords(regions::R_EX_0, regions::S_EX_1, u) );
        raw_jpoints.push_back( JacobiPoint::from_jacobicoords(regions::R_EX_1, regions::S_EX_0, u) );
        raw_jpoints.push_back( JacobiPoint::from_jacobicoords(regions::R_EX_1, regions::S_EX_1, u) );
        raw_jpoints.push_back( JacobiPoint::from_jacobicoords(R              , regions::S_EX_0, u) );
        raw_jpoints.push_back( JacobiPoint::from_jacobicoords(R              , regions::S_EX_1, u) );
    }
    else {
        std::cerr << "ERROR in " << __func__;
        std::cerr << ": unreachable; the JacobiPoint should have fallen into";
        std::cerr << " one of the previous cases";
        std::cerr << std::endl;
        exit(EXIT_FAILURE);
    }
    
    return raw_jpoints;
}

}
