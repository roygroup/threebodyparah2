#include <algorithm>
#include <vector>
#include <cmath>
#include <tuple>

#include <iostream>

#include "jpoint.h"

JacobiPoint JacobiPoint::from_unsorted_pairdistances(double rAB, double rBC, double rAC) {
    // sort the distances so r12 <= r23 <= r13
    std::vector<double> distances {rAB, rAC, rBC};
    std::sort(distances.begin(), distances.end());
    double r12 = distances[0];
    double r23 = distances[1];
    double r13 = distances[2];
    
    // convert r12, r23, r13 to the rescaled Jacobi coordinates
    double _R = r12;
    double _r = sqrt( 0.5 * (r13*r13 + r23*r23 - 0.5*r12*r12) );
    double _cost = (r13*r13 - r23*r23) / (2.0*r12*_r);
    
    if      (_cost > 1.0)
        _cost = 1.0;
    else if (_cost < 0.0)
        _cost = 0.0;

    double _t = acos(_cost);

    // get the scaled coordinates
    double _rmin = 0.5 * _R * ( _cost + sqrt(3.0 + _cost*_cost) );
    double _s    = _r / _rmin;
    double _u = _t / (M_PI/2.0);

    // floating point imprecision sometimes puts s slightly below 1;
    if (_s < 1.0)
        _s = 1.0;
    
    return JacobiPoint(_R, _s, _u);
}

JacobiPoint JacobiPoint::from_jacobicoords(double _R, double _s, double _u) {
    return JacobiPoint(_R, _s, _u);
}

// ----------------------------------------------------------------------------

// constructor directly from Jacobi coordinates
JacobiPoint::JacobiPoint(double _R, double _s, double _u) {
    R = _R;
    s = _s;
    u = _u;
}

// copy constructor
JacobiPoint::JacobiPoint(const JacobiPoint& other) {
    R = other.get_R();
    s = other.get_s();
    u = other.get_u();
}

// assignment operator
JacobiPoint& JacobiPoint::operator=(const JacobiPoint& other) {
    if (this == &other) {
        return *this;
    }
    R = other.get_R();
    s = other.get_s();
    u = other.get_u();
    return *this;
}

// destructor
JacobiPoint::~JacobiPoint() {}

// getters for the Jacobi coordinates
double JacobiPoint::get_R() const { return R; }
double JacobiPoint::get_s() const { return s; }
double JacobiPoint::get_u() const { return u; }

std::tuple<double, double, double> JacobiPoint::unpack() const {
    return std::make_tuple(R, s, u);
}
