#include "rkhs.h"
#include "jpoint.h"

// FORTRAN FILES
extern "C" {
    void calculate_energies_(char*, double*, double*, double*, double*, int*);
}

namespace rkhs
{
    std::vector<double> get_raw_energies(std::string& kernelfile, const std::vector<JacobiPoint>& raw_jpoints) {
    	// call the Fortran code to calculate the interpolated three-body energies using RKHS
    	
    	// unpack the JacobiPoint coordinates into three separate vectors
        int N = raw_jpoints.size();
    	std::vector<double> R_arr; R_arr.reserve(N);
    	std::vector<double> s_arr; s_arr.reserve(N);
    	std::vector<double> u_arr; u_arr.reserve(N);
    	for (std::size_t i = 0; i < N; ++i) {
    		R_arr.push_back(raw_jpoints[i].get_R());
    		s_arr.push_back(raw_jpoints[i].get_s());
    		u_arr.push_back(raw_jpoints[i].get_u());
    	}
    	
    	// calculate the RKHS energies using the fortran code
    	std::vector<double> raw_energies(N);
        char* c_kernelfile = (char*)kernelfile.c_str();
        calculate_energies_(c_kernelfile, R_arr.data(), s_arr.data(), u_arr.data(), raw_energies.data(), &N);
    	
    	return raw_energies;
    }
}
