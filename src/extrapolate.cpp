#include <functional>
#include <iostream>
#include <tuple>
#include <cmath>

#include "utils.h"
#include "region.h"
#include "extrapolate.h"
#include "trapezoid.h"

namespace extrapolate {

// --- EXTRAPOLATION FOR Region::AxilrodTellerMuto --------------------------------

double energy_AxilrodTellerMuto(const JacobiPoint& jpoint) {
    return atm::energy(jpoint);
}

// --- EXTRAPOLATION FOR Region::Bounded ------------------------------------------

double energy_Bounded(const JacobiPoint& jpoint, double energy_raw) {
    double R = jpoint.get_R();
    double s = jpoint.get_s();
    
    double final_energy;
    if (R >= 5.35) {
        double energy_atm = atm::energy(jpoint);
        double weight_atm = utils::cosine_transition(R, 5.35, 6.25);
        double weight_raw = 1.0 - weight_atm;
        final_energy = weight_atm*energy_atm + weight_raw*energy_raw;
    }
    else if (R > 3.25 && s > 3.5) {
        double energy_atm = atm::energy(jpoint);
        double weight_atm = utils::cosine_transition(s, 3.5, 3.85);
        double weight_raw = 1.0 - weight_atm;
        final_energy = weight_atm*energy_atm + weight_raw*energy_raw;
    }
    else {
        final_energy = energy_raw;
    }

    return final_energy;
}

// --- EXTRAPOLATION FOR Region::s_extrapolation ----------------------------------

double dfds_nonvanish(double s, double u) {
    double cosu = cos(HALFPI*u);
    double w  = atm::wterm(s, u);
    double w2 = w*w;

    double denom_ = atm::tterm(s, u);
    double denom  = denom_*denom_;
    double numer  = 3.0 + 2.0*w - (1.0 + 4.0*cosu*cosu)*w2;
    double coeff  = 6.0 * w / s;

    return coeff*numer/denom;
}

std::tuple<double, double> get_a_and_b(double R, double u, double s0, double s1, double e0, double e1, unsigned int N) {
    JacobiPoint jpoint0 = JacobiPoint::from_jacobicoords(R, s0, u);
    JacobiPoint jpoint1 = JacobiPoint::from_jacobicoords(R, s1, u);

    double U0 = atm::norm(jpoint0) * e0;
    double U1 = atm::norm(jpoint1) * e1;
    double f0 = atm::fterm(jpoint0);
    double f1 = atm::fterm(jpoint1);

    std::function<double(double)> dfds_nonv = [=](double s_) {
        return dfds_nonvanish(s_, u);
    };

    auto trapintegrator = TrapezoidIntegrator1D(dfds_nonv);
    
    double integF10 = trapintegrator.integrate(s0, s1, N);
    
    double sinu = sin(HALFPI*u);
    double inv_sinu_sq = 1.0 / (sinu*sinu);

    double a = inv_sinu_sq * (U0*f1 - U1*f0) / integF10;
    double b = inv_sinu_sq * (U1 - U0) / integF10;

    return std::make_tuple(a, b);
}

double energy_s_extrapolation(const JacobiPoint& jpoint, double s0, double s1, double e0, double e1) {
    double R = jpoint.get_R();
    double u = jpoint.get_u();
    
    double a, b;
    std::tie(a, b) = get_a_and_b(R, u, s0, s1, e0, e1, 1000);
    
    double energy_mod = (a + b*atm::fterm(jpoint)) * atm::invnorm(jpoint);
    
    double energy_final;
    if (R > 3.15) {
        double energy_atm = atm::energy(jpoint);
        double weight_atm = utils::cosine_transition(R, 3.15, 3.25);
        double weight_mod = 1.0 - weight_atm;
        
        energy_final = weight_atm*energy_atm + weight_mod*energy_mod;
    }
    else {
        energy_final = energy_mod;
    }

    return energy_final;
}


// --- EXTRAPOLATION FOR Region::R_extrapolation ----------------------------------

RVariableExtrapolator::RVariableExtrapolator(double e0, double e1, double R0, double R1)
    : e0(e0), e1(e1), R0(R0), R1(R1)
{
    dR = R1 - R0;
    c_lin = (e1 - e0)/dR;
}

double RVariableExtrapolator::linear(double R) const {
    return e0 + c_lin*(R - R0);
}

double RVariableExtrapolator::exponential(double c_exp, double R) const {
    return e0 * exp(-c_exp*(R - R0));
}

double RVariableExtrapolator::extrapolate(double R) const {
    double extrapolated_energy;
    if (e0*e1 > 0.0) {
        // hackish way of checking if e0 and e1 have the same sign
        double a0 = std::max(abs_eng_min, fabs(e0));
        double a1 = std::max(abs_eng_min, fabs(e1));

        if (a1 >= a0) {
            extrapolated_energy = linear(R);
        }
        else {
            double c_exp = - (1.0/dR) * log(a1/a0);

            if      (c_exp <= c_min) {
                extrapolated_energy = exponential(c_exp, R);
            }
            else if (c_exp >= c_max) {
                extrapolated_energy = linear(R);
            }
            else {
                double weight_lin = utils::cosine_transition(c_exp, c_min, c_max);
                double energy_lin = linear(R);
                double weight_exp = 1.0 - weight_lin;
                double energy_exp = exponential(c_exp, R);

                extrapolated_energy = weight_lin*energy_lin + weight_exp*energy_exp;
            }
        }
    }
    else {
        extrapolated_energy = linear(R);
    }

    return extrapolated_energy;
}

double energy_R_extrapolation(const JacobiPoint& jpoint, double R0, double R1, double e0, double e1, double e) {
    double R, s, u;
    std::tie(R, s, u) = jpoint.unpack();

    RVariableExtrapolator rextrapolator = RVariableExtrapolator(e0, e1, R0, R1);
    
    double weight_pes = utils::cosine_transition(R, 2.0, 2.2);
    double energy_pes = e;
    
    double weight_ext = 1.0 - weight_pes;
    double energy_ext = rextrapolator.extrapolate(R);

    return weight_pes*energy_pes + weight_ext*energy_ext;
}

// --- EXTRAPOLATION FOR Region::R_and_s_extrapolation ----------------------------

double energy_R_and_s_extrapolation(const JacobiPoint& jpoint, double R0, double R1, double s0, double s1, const std::vector<double>& raw_energies) {
    double R, s, u;
    std::tie(R, s, u) = jpoint.unpack();

    double e0, e1, e2, e3, e4, e5;
    e0 = raw_energies[0];
    e1 = raw_energies[1];
    e2 = raw_energies[2];
    e3 = raw_energies[3];
    e4 = raw_energies[4];
    e5 = raw_energies[5];
    
    // extrapolation above \varphi_cf
    u = std::max(u, regions::UCHANGE);

    JacobiPoint jpoint01 = JacobiPoint::from_jacobicoords(R0, s, u);
    JacobiPoint jpoint23 = JacobiPoint::from_jacobicoords(R1, s, u);
    JacobiPoint jpoint45 = JacobiPoint::from_jacobicoords(R1, s, u);

    double e01 = energy_s_extrapolation(jpoint01, s0, s1, e0, e1);
    double e23 = energy_s_extrapolation(jpoint23, s0, s1, e2, e3);
    double e45 = energy_s_extrapolation(jpoint45, s0, s1, e4, e5);

    return energy_R_extrapolation(jpoint, R0, R1, e01, e23, e45);
}

// --- GENERAL EXTRAPOLATION ENERGY -----------------------------------------------

double processed_energy(const JacobiPoint& jpoint, const regions::Region& reg, const std::vector<double>& raw_energies, std::size_t i_raw_start) {
    double proc_energy;
    if (reg == regions::Region::AxilrodTellerMuto) {
        proc_energy = energy_AxilrodTellerMuto(jpoint);
        //std::cout << "(R, s, u) = (" << jpoint.get_R() << ", " << jpoint.get_s() << ", " << jpoint.get_u() << ") : ";
        //std::cout << "AXILROD" << std::endl;
    }
    else if (reg == regions::Region::Bounded) {
        double e0 = raw_energies[i_raw_start];
        proc_energy = energy_Bounded(jpoint, e0);
        //std::cout << "(R, s, u) = (" << jpoint.get_R() << ", " << jpoint.get_s() << ", " << jpoint.get_u() << ") : ";
        //std::cout << "BOUNDED" << std::endl;
    }
    else if (reg == regions::Region::s_extrapolation) {
        double e0 = raw_energies[i_raw_start];
        double e1 = raw_energies[i_raw_start + 1];
        proc_energy = energy_s_extrapolation(jpoint, regions::S_EX_0, regions::S_EX_1, e0, e1);
        //std::cout << "(R, s, u) = (" << jpoint.get_R() << ", " << jpoint.get_s() << ", " << jpoint.get_u() << ") : ";
        //std::cout << "s-EXTRAPOLATION" << std::endl;
    }
    else if (reg == regions::Region::R_extrapolation) {
        double e0 = raw_energies[i_raw_start];
        double e1 = raw_energies[i_raw_start + 1];
        double e  = raw_energies[i_raw_start + 2];
        proc_energy = energy_R_extrapolation(jpoint, regions::R_EX_0, regions::R_EX_1, e0, e1, e);
        //std::cout << "(R, s, u) = (" << jpoint.get_R() << ", " << jpoint.get_s() << ", " << jpoint.get_u() << ") : ";
        //std::cout << "R-EXTRAPOLATION" << std::endl;
    }
    else if (reg == regions::Region::R_and_s_extrapolation) {
        std::vector<double> Rs_raw_energies = {raw_energies.begin() + i_raw_start, raw_energies.begin() + i_raw_start + 6};
        proc_energy = energy_R_and_s_extrapolation(jpoint, regions::R_EX_0, regions::R_EX_1, regions::S_EX_0, regions::S_EX_1, Rs_raw_energies);
        //std::cout << "(R, s, u) = (" << jpoint.get_R() << ", " << jpoint.get_s() << ", " << jpoint.get_u() << ") : ";
        //std::cout << "R-and-s-EXTRAPOLATION" << std::endl;
    }
    else {
        std::cerr << "ERROR in " << __func__;
        std::cerr << ": unreachable; the JacobiPoint should have fallen into";
        std::cerr << " one of the previous cases";
        std::cerr << std::endl;
        exit(EXIT_FAILURE);
    }

    return proc_energy;
}

}
