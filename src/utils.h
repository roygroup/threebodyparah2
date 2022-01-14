#ifndef UTILS_H
#define UTILS_H

#include <cmath>

#include "jpoint.h"

constexpr double HALFPI = M_PI / 2.0;

namespace utils
{
    double cosine_transition(double, double, double);
}

namespace atm
{
    constexpr double C9     = 34336.220013464925;

    double wterm(const JacobiPoint&);
    double wterm(double, double);
    double tterm(const JacobiPoint&);
    double tterm(double, double);
    double fterm(const JacobiPoint&);

    double norm(const JacobiPoint&);
    double invnorm(const JacobiPoint&);

    double energy(const JacobiPoint&);
}

#endif
