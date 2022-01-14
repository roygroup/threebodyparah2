'''
    Apply the long-range R-coordinate adjustments to the AVTZ data.
    NOTE: the script `make_smodAVTZ_csvfile.py` already performed many of the long-range
          R-coordinate adjustments; this script performs the remaining ones, and corrects
          for imprecise adjustments made by `make_smodAVTZ_csvfile.py`

    [s-index] [t-index] [index for R_A] [decay constant] ...[indices of R-values to apply the decay to]

    If ([index for R_A] == -1) and ([decay constant] == 0.0), then the indices to apply the
    decay to are being set directly to the ATM potential.
'''

import os
import numpy as np
import rst_data as rst
import threebody_potential_types as tbp

def rkhs_along_R(Rdata, s, t, kerneltag):
    rkhs = tbp.RKHSPotential(kerneltag, 'threebodyparah2')
    rstpoints = [tbp.RstPoint(R, s, t) for R in Rdata]

    return rkhs.from_RstPoint_list(rstpoints)

def atm_along_R(Rdata, s, t):
    atm = tbp.AxilrodTellerMuto()
    rstpoints = [tbp.RstPoint(R, s, t) for R in Rdata]

    return atm.from_RstPoint_list(rstpoints)

class NewEnergyData:
    def __init__(self, line, offset):
        # index offset accounts for extra lines at start of the csvfile
        self.offset = offset

        tokens = line.split()
        self.i_s       = int(tokens[0])
        self.i_t       = int(tokens[1])
        self.i_R0      = int(tokens[2])
        self.alpha     = float(tokens[3])
        self.Rindices_ = [int(iR) for iR in tokens[4:]]

        s = rst.get_sdata()[self.i_s]
        t = rst.get_tdata()[self.i_t]
        Rdata = rst.get_Rdata()

        self.indices_and_new_energies = []
        
        kernelpath  = os.path.join('..', 'kernels', 'AVTZdata.kernel')
        ccsdt_point = rkhs_along_R(Rdata, s, t, kernelpath)
        atm_point   = atm_along_R(Rdata, s, t)

        if self.i_R0 == -1:
            for i_R in self.Rindices_:
                atm_eng = atm_point[i_R]
                self.indices_and_new_energies.append((self.get_index(i_R), atm_eng))
        else:
            R0 = Rdata[self.i_R0]
            ccsdt_eng0 = ccsdt_point[self.i_R0]
            atm_eng0   = atm_point[self.i_R0]
            ratio0     = np.abs((atm_eng0 - ccsdt_eng0)/ccsdt_eng0)
            
            self.indices_and_new_energies.append((self.get_index(self.i_R0), atm_eng0))
    
            for i_R in self.Rindices_:
                R    = Rdata[i_R]
                frac = self.decay_function(R, R0, self.alpha)
                ccsdt_eng  = ccsdt_point[i_R]
                atm_eng    = atm_point[i_R]
                new_energy = ccsdt_eng * (1.0 + ratio0*frac)
                self.indices_and_new_energies.append((self.get_index(i_R), new_energy))

    def decay_function(self, R, R0, alpha):
        return np.exp(-0.5*alpha*(R - R0)**2)
            
    def get_index(self, i_R):
        i_Rstep = 323 * i_R
        i_sstep = 19  * self.i_s
        i_tstep = 1   * self.i_t
        return i_Rstep + i_sstep + i_tstep + self.offset
    
    def get_indices_and_new_energies(self):
        return self.indices_and_new_energies

def replace_csvline_energy(csvline, new_energy):
    tokens = csvline.split(',   ')
    tmpl_line    = "{0:s},   {1:s},   {2:s},   {3: 0.15f}\n"
    return tmpl_line.format(tokens[0], tokens[1], tokens[2], new_energy)

# --------------------------------------------------------------------------------

if __name__ == '__main__':
    neweng_filename = "long_R_adjustments.txt"
    offset = 0

    # read in the new energy
    new_energy_datas = []
    with open(neweng_filename, 'r') as fin:
        for line in fin:
            ned = NewEnergyData(line, offset)
            new_energy_datas.append(ned)
    
    # read in the csv file lines
    csvfilename = os.path.join('..', 'csvfiles', 'smodAVTZ.csv')
    csvlines = open(csvfilename, 'r').readlines()

    for ned in new_energy_datas:
        for (index, new_energy) in ned.get_indices_and_new_energies():
            old_csvline = csvlines[index]
            new_csvline = replace_csvline_energy(old_csvline, new_energy)
            csvlines[index] = new_csvline
    
    new_csvfilename = os.path.join('..', 'csvfiles', 'RsmodAVTZ.csv')
    with open(new_csvfilename, 'w') as fout:
        for line in csvlines:
            fout.write(line)

