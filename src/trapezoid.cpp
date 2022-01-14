#include <functional>
#include <cassert>
#include <iostream>

#include "trapezoid.h"

TrapezoidIntegrator1D::TrapezoidIntegrator1D(const std::function<double(double)>& func)
    : func(func)
{
}


double TrapezoidIntegrator1D::integrate(double xmin, double xmax, unsigned int size) const {
    if (size <= 1) {
        std::cerr << "ERROR in " << __func__;
        std::cerr << ": need at least 2 points to perform trapezoidal intergation" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    if (xmax <= xmin) {
        std::cerr << "ERROR in " << __func__;
        std::cerr << ": left endpoint greater than right endpoint" << std::endl;
        exit(EXIT_FAILURE);
    }

    double dx = (xmax - xmin) / (double)(size - 1);
    
    double endpoint_contrib = 0.5 * (func(xmin) + func(xmax));

    double midpoint_contrib = 0.0;
    for (std::size_t i = 1; i < size - 1; ++i) {
        double x = xmin + i*dx;
        midpoint_contrib += func(x);
    }

    double total_area = dx * (endpoint_contrib + midpoint_contrib);

    return total_area;
}

double TrapezoidIntegrator1D::eval(double x) const {
    return func(x);
}

