#ifndef TRAPEZOID_H
#define TRAPEZOID_H

#include <functional>

class TrapezoidIntegrator1D
{
private:
    std::function<double(double)> func;
public:
    TrapezoidIntegrator1D(const std::function<double(double)>&);
    
    double integrate(double, double, unsigned int) const;
    double eval(double) const;
};


#endif
