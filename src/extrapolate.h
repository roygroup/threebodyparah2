#ifndef EXTRAPOLATE_H
#define EXTRAPOLATE_H

#include <functional>
#include <tuple>

#include "utils.h"

namespace extrapolate
{

class RVariableExtrapolator
{
private:
    const double c_min = 6.0;   // lower linear extrapolation limit
    const double c_max = 8.0;   // upper linear extrapolation limit
    const double abs_eng_min = 1.0e-8;

    const double R0 = regions::R_EX_0;
    const double R1 = regions::R_EX_1;

    const double e0, e1;  // RKHS energies at the 0th and 1st extrapolation points
    double c_lin;   // linear slope
    double dR;      // difference between the 0th and 1st extrapolaiton R-values

    double linear(double) const;
    double exponential(double, double) const;
public:
    RVariableExtrapolator(double, double, double, double);
    double extrapolate(double) const;
};

double energy_AxilrodTellerMuto(const JacobiPoint&);
double energy_Bounded(const JacobiPoint&, double);

double                     dfds_nonvanish(double, double);
std::tuple<double, double> get_a_and_b(double, double, double, double, double, double, unsigned int);
double energy_s_extrapolation(const JacobiPoint&, double, double, double, double);

double energy_R_extrapolation(const JacobiPoint&, double, double, double, double, double);
double energy_R_and_s_extrapolation(const JacobiPoint&, double, double, double, double, const std::vector<double>&);

// --- GENERAL EXTRAPOLATION ENERGY -----------------------------------------------

double processed_energy(const JacobiPoint&, const regions::Region&, const std::vector<double>&, std::size_t);

}

#endif
