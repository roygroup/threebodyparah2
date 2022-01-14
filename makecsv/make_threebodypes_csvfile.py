'''
    Include the short-range R-coordinate phantom points.
    The purpose is only to make sure the short-range extrapolation is done smoothly.
'''

import os
import itertools
import numpy as np
import callpes as cp
import rst_data as rst
import threebody_potential_types as tbp
import pes_regions as pr

# ---------------------------------------------------------------------------------

def get_threebody_energies_no_short_range_mix(rstpoint_list, kerneltag):
    rkhspot = tbp.RKHSPotential(kerneltag, 'threebodyparah2')

    # each region needs a diff number of triplets to calculate an energy within it
    regions = cp.get_regions(rstpoint_list)
    rstpoints = []
    for region in regions:
        rstpoints.extend(region.get_rstpoints())
    
    raw_energies = rkhspot.from_RstPoint_list(rstpoints)
    
    # the "raw" energies all lie within the input data mesh region, and
    # will be used in extrapolations to get the final energies
    final_energies = get_extrapolated_energies(regions, raw_energies)

    return final_energies

def get_extrapolated_energies(regions, energies):
    extrapolated_energies = []

    lower_cap = 0
    upper_cap = 0
    for region in regions:
        upper_cap += region.get_capacity()

        region.set_energies(energies[lower_cap:upper_cap])
        region.calculate_final_energy_no_mix()

        lower_cap += region.get_capacity()

        final_energy = region.get_final_energy()
        extrapolated_energies.append(final_energy)

    return extrapolated_energies

# ---------------------------------------------------------------------------------

def read_energy_on_nth_line(n, csvfile):
    with open(csvfile, 'r') as fin:
        for (i, line) in enumerate(fin):
            if i == n:
                tokens = line.split(',')
                energy = float(tokens[-1])

    return energy

def get_R0_R1_energies(i_s, i_t, Ns, Nt, csvfile):
    i_R0_index =         Nt*i_s + i_t
    i_R1_index = Nt*Ns + Nt*i_s + i_t

    energy_R0 = read_energy_on_nth_line(i_R0_index, csvfile)
    energy_R1 = read_energy_on_nth_line(i_R1_index, csvfile)

    return energy_R0, energy_R1

def get_extrapolated_energy(R, i_s, i_t, Ns, Nt, csvfile):
    e0, e1 = get_R0_R1_energies(i_s, i_t, Ns, Nt, csvfile)
    R0 = 2.20
    R1 = 2.25

    extrapolator = pr.ShortRangeExtrapolator(e0, e1, R0, R1)
    return extrapolator.extrap(R)


# ---------------------------------------------------------------------------------

def main():
    short_Rdata = np.array([2.00, 2.05, 2.10, 2.15])
    sdata       = rst.get_sdata()
    tdata       = rst.get_tdata()
    Ns = len(sdata)
    Nt = len(tdata)

    csv_line_tmpl = "{0:0.15f},   {1:0.15f},   {2:0.15f},   {3: 0.15f}\n"
    input_csvfile  = os.path.join('..', 'csvfiles', 'RsmodAVTZ.csv')
    output_csvfile = os.path.join('.', 'threebodypes.csv')
    
    with open(output_csvfile, 'w') as fout:
        for R in short_Rdata:
            for (i_s, s) in enumerate(sdata):
                for (i_t, t) in enumerate(tdata):
                    energy = get_extrapolated_energy(R, i_s, i_t, Ns, Nt, input_csvfile)
                    t_norm = t / (np.pi/2.0)
                    line = csv_line_tmpl.format(R, s, t_norm, energy)
                    fout.write(line)

        with open(input_csvfile, 'r') as fin:
            for line in fin:
                fout.write(line)

if __name__ == '__main__':
    main()




