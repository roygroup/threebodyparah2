// VERSION 1.0
// --------------------------------------------------------------------------------
// MIT License
// 
// Copyright (c) 2021 Alexander Ibrahim
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>

#include "rkhs.h"
#include "jpoint.h"
#include "region.h"
#include "extrapolate.h"

std::string usage_message(std::string& executable) {
	// creates a message that tells the user how to use the executable
	std::stringstream message;

	message << "Usage: " << executable << " KERNELFILE INPUTFILE > OUTPUTFILE        \n";
	message << "                                                                     \n";
	message << "  KERNELFILE: name of .kernel file                                   \n";
	message << "                                                                     \n";
	message << "  INPUTFILE: File of three whitespace-delimited columns where        \n";
	message << "                                                                     \n";
	message << "    column 1: distance between hydrogen molecules 1 and 2            \n";
	message << "    column 2: distance between hydrogen molecules 2 and 3            \n";
	message << "    column 3: distance between hydrogen molecules 1 and 3            \n";
	message << "                                                                     \n";
	message << "    All distances are in Angstroms.                                  \n";
	message << "                                                                     \n";
	message << "  OUTPUTFILE: Single-column file of three-body parahydrogen energies.\n";
	message << "                                                                     \n";
	message << "    All energies are in cm^{-1}.                                     \n";
	
	return message.str();
}
unsigned int total_number_raw_jpoints(const std::vector<regions::Region>& region_list) {
    unsigned int N = 0;
    for (auto reg : region_list) {
        N += regions::number_of_raw_JacobiPoints(reg);
    }

    return N;
}

unsigned int number_of_lines(const std::string& filename) {
	// count the number of lines in `filename`
	unsigned int nlines = 0;
	
	std::string line;
	std::ifstream fin(filename);
	while (std::getline(fin, line))
		++nlines;
	
	fin.close();
	
	return nlines;
}

std::vector<JacobiPoint> get_input_jpoints(const std::string& inputfilename) {
    // prepare the vector of JacobiPoint instances
    int nlines = number_of_lines(inputfilename);
    std::vector<JacobiPoint> jpoints;
    jpoints.reserve(nlines);
    
    // read the pair distances from the three columns of the input file
    std::string line;
    std::stringstream ssvalue;
    std::ifstream fin(inputfilename);
    
    for (int i = 0; i < nlines; ++i) {
        // every line of the input file contains pair distances {r12, r23, r13}
        std::getline(fin, line);
        ssvalue << line;
        
        double r12, r23, r13;
        try {
            ssvalue >> r12;
            ssvalue >> r23;
            ssvalue >> r13;
        }
        catch (...) {
            std::cerr << "ERROR in " << __func__;
            std::cerr << ": unable to read pair distance values from " << inputfilename << std::endl;
            exit(EXIT_FAILURE);
        }
        
        // create a new JacobiPoint instance
        jpoints.push_back(JacobiPoint::from_unsorted_pairdistances(r12, r23, r13));
        
        // reset ssvalue stringstream
        ssvalue.str(std::string());
        ssvalue.clear();
    }

    fin.close();

    return jpoints;
}

// --------------------------------------------------------------------------------

std::vector<double> get_threebody_energies(std::string& kernelfile, const std::vector<JacobiPoint>& jpoint_list) {
    // find the Region of (R, s, \varphi) coordinate space each point belongs to
    std::vector<regions::Region> region_list;
    region_list.reserve(jpoint_list.size());
    for (auto jpoint : jpoint_list) {
        regions::Region reg = regions::which_region(jpoint);
        region_list.push_back(reg);
    }
    
    // use the actual JacobiPoint instances to get the `raw` JacobiPoint instances needed to get the energies
    std::vector<JacobiPoint> raw_jpoint_list;
    raw_jpoint_list.reserve(total_number_raw_jpoints(region_list));
    for (std::size_t i = 0; i < jpoint_list.size(); ++i) {
        JacobiPoint     jpoint = jpoint_list[i];
        regions::Region reg    = region_list[i];

        std::vector<JacobiPoint> jpoints_raw = regions::get_raw_jpoints(jpoint, reg);
        raw_jpoint_list.insert(raw_jpoint_list.end(), jpoints_raw.begin(), jpoints_raw.end());
    }
    
    // need the `raw` energies from the RKHS interpolation method
    std::vector<double> raw_rkhs_energies = rkhs::get_raw_energies(kernelfile, raw_jpoint_list);
    
    // the use `raw` energies to calculate the final threebody energies
    std::size_t i_raw = 0;

    std::vector<double> processed_energies;
    processed_energies.reserve(jpoint_list.size());
    for (std::size_t i = 0; i < jpoint_list.size(); ++i) {
        regions::Region reg = region_list[i];
        JacobiPoint jpoint  = jpoint_list[i];

        processed_energies.push_back(extrapolate::processed_energy(jpoint, reg, raw_rkhs_energies, i_raw));

        std::size_t N_raw_energies_for_region = regions::number_of_raw_JacobiPoints(reg);
        i_raw += N_raw_energies_for_region;
    }
    
    return processed_energies;
}

int main(int argc, char** argv) {
    // parse the command line arguments
    std::string kernelfile, inputfile;
    
    if (argc != 3) {
        std::string executable = argv[0];
        std::cerr << usage_message(executable) << std::endl;
        exit(EXIT_FAILURE);
    }
    else {
        try {
            kernelfile = argv[1];
            inputfile  = argv[2];
        }
        catch (...) {
            std::cerr << "ERROR in " << __func__;
            std::cerr << ": unable to read in the .kernel file and the input file" << std::endl;
    
            std::string executable = argv[0];
            std::cerr << usage_message(executable) << std::endl;
    
            exit(EXIT_FAILURE);
        }
    }
    
    std::vector<JacobiPoint> jpoints = get_input_jpoints(inputfile);
    std::vector<double> energies = get_threebody_energies(kernelfile, jpoints);

    // set output precision and type before outputting energies
    std::cout.precision(12);
    std::cout << std::scientific;
    for (auto energy : energies) {
        std::cout << energy << std::endl;
    }

    return 0;
}
