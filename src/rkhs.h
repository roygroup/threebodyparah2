#ifndef RKHS_H
#define RKHS_H

#include <vector>
#include <string>

#include "jpoint.h"

namespace rkhs
{
    std::vector<double> get_raw_energies(std::string&, const std::vector<JacobiPoint>&);
}

#endif
