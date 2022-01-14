'''
    Apply the long-range s-coordinate adjustments to the AVTZ data.
    NOTE: this performs both the long-range s-coordinate adjustments, and part of the
          long-range R-coordinate adjustments, because at large enough R-values, all
          the adjustments for the MATM potential converge to adjustments for the ATM potential

    The file `sprime_list.txt` contains the data for the adjustments.
    The adjustments fall into 3 categories:
        > [R-value] [t-value] [s-prime]
            > the chosen value of s-prime determines a(R, t) and b(R, t) for the MATM potential
            > all mesh values of s > s-prime will be replaced by the MATM potential
        > [R-value] [t-value] 0 [s-swap]
            > all mesh values of s < s-swap will be determines by the ab initio data
            > all mesh values of s > s-prime will be replace by the ATM potential
        > [R-value] [t-value] 0 1.0
            > all mesh values will be replaced by the ATM potential
'''

import os
import itertools
import numpy as np
import rst_data as rst
import threebody_potential_types as tbp

# ---------------------------------------------------------------------------------

def normalized_rkhs(kerneltag, R, sdata, t):
    rkhs = tbp.RKHSPotential(kerneltag, 'threebodyparah2')
    rstpoints = [tbp.RstPoint(R, s, t) for s in sdata]
    
    atm = tbp.AxilrodTellerMuto()
    normdata = np.array([atm.norm(rstp) for rstp in rstpoints])

    return normdata * rkhs.from_RstPoint_list(rstpoints)

# ---------------------------------------------------------------------------------

class MATMReplacementData:
    eps = 1.0e-5

    def __init__(self, matm_line, sdata, ccsdt_engs):
        self.elements = [float(elem) for elem in matm_line.split()]

        if len(self.elements) == 3:
            self.elements = self.elements + [self.elements[-1]]

        self.sdata = sdata
        self.ccsdt = ccsdt_engs
        self.R      = self.elements[0]
        self.t      = self.elements[1]
        self.s_matm = self.elements[2]
        self.s_swap = self.elements[3]
        
        if (np.abs(self.t) < self.eps):
            self.t = self.eps
        if self.t > np.pi/2.0:
            self.t = np.pi/2.0
        self.rstpoints = [tbp.RstPoint(self.R, s, self.t) for s in self.sdata]

    def get_weighted_energy(self, matm_eng, ccsdt_eng, s, s_swap):
        if s >= s_swap:
            return matm_eng
        else:
            return ccsdt_eng

    # -----------------------------------------------------------------------------
    # CASE 1:  R  t   sMATM     sMATM        : switch from CC to MATM at sMATM
    def is_case1(self):
        return (np.abs(self.elements[2]) > self.eps) and (self.elements[3] > 1.0)
    
    def case1_energies(self, ccsdt_engs, sdata, R, t, s_matm):
        # get the MATM energies
        fine_sdata = np.linspace(1.0, 3.85, 500)
        kernelpath = os.path.join('..', 'kernels', 'AVTZdata.kernel')
        fine_norm_ccsdt = normalized_rkhs(kernelpath, R, fine_sdata, t)
        matm      = tbp.ModifiedAxilrodTellerMuto(fine_norm_ccsdt, fine_sdata, t, s_matm)
        matm_engs = matm.from_RstPoint_list(self.rstpoints)

        # combine the CCSD(T) and MATM energies
        wgt_engs = np.empty(len(sdata), dtype = float)
        for i, (s, ccsdt_eng, matm_eng) in enumerate(zip(sdata, ccsdt_engs, matm_engs)):
            wgt_engs[i] = self.get_weighted_energy(matm_eng, ccsdt_eng, s, s_matm)

        return wgt_engs

    # -----------------------------------------------------------------------------
    # CASE 2:  R  t   0.0     s[CC_TO_ATM]   : switch from CC to ATM at s[CC_TO_ATM]
    def is_case2(self):
        return (np.abs(self.elements[2]) < self.eps) and (self.elements[3] > 1.0)

    def case2_energies(self, ccsdt_engs, sdata, s_swap):
        # get the ATM energies
        atm      = tbp.AxilrodTellerMuto()
        atm_engs = atm.from_RstPoint_list(self.rstpoints)
        
        # combine the CCSD(T) and MATM energies
        wgt_engs = np.empty(len(sdata), dtype = float)
        for i, (s, ccsdt_eng, atm_eng) in enumerate(zip(sdata, ccsdt_engs, atm_engs)):
            wgt_engs[i] = self.get_weighted_energy(atm_eng, ccsdt_eng, s, s_swap)

        return wgt_engs

    # -----------------------------------------------------------------------------
    # CASE 3:  R  t   0.0     1.0            : 100 percent ATM
    def is_case3(self):
        return (np.abs(self.elements[2]) < self.eps) and (np.abs(self.elements[3] - 1.0) < self.eps)

    def case3_energies(self):
        # get the ATM energies
        atm      = tbp.AxilrodTellerMuto()
        atm_engs = atm.from_RstPoint_list(self.rstpoints)

        return atm_engs
    
    def get_modified_energies(self):
        if self.is_case1():
            return self.case1_energies(self.ccsdt, self.sdata, self.R, self.t, self.s_matm)
        if self.is_case2():
            return self.case2_energies(self.ccsdt, self.sdata, self.s_swap)
        elif self.is_case3():
            return self.case3_energies()
        else:
            assert False, "unreachable"

# ---------------------------------------------------------------------------------

class CSVEnergyManager:
    def __init__(self):
        self.Rdata = rst.get_Rdata()
        self.sdata = rst.get_sdata()
        self.tdata = rst.get_tdata()

        self.Rindices = np.arange(len(self.Rdata))
        self.sindices = np.arange(len(self.sdata))
        self.tindices = np.arange(len(self.tdata))

    def extract_energy(self, line):
        values = line.split(',')
        energy = float(values[3])
    
        return energy
    
    def read_csv_energies(self, csvfilename):
        energies = np.empty((len(self.Rdata), len(self.sdata), len(self.tdata)), dtype = float)

        with open(csvfilename, 'r') as fin:
            for (i_R, i_s, i_t) in itertools.product(self.Rindices, self.sindices, self.tindices):
                line = fin.readline()
                energies[i_R][i_s][i_t] = self.extract_energy(line)

        return energies

    def write_csv_energies(self, energies, csvfilename):
        assert energies.shape == (len(self.Rdata), len(self.sdata), len(self.tdata))
        
        FMT_line_eng = "{0:0.15f},   {1:0.15f},   {2:0.15f},   {3: 0.15f}\n"

        with open(csvfilename, 'w') as fout:
            for (i_R, i_s, i_t) in itertools.product(self.Rindices, self.sindices, self.tindices):
                R = self.Rdata[i_R]
                s = self.sdata[i_s]
                t = self.tdata[i_t] / (np.pi/2.0)
                e = energies[i_R][i_s][i_t]
                fout.write(FMT_line_eng.format(R, s, t, e))

def get_s_energies(energies, i_R, i_t):
    s_energies = np.empty(energies.shape[1], dtype = float)
    for i_s in range(energies.shape[1]):
        s_energies[i_s] = energies[i_R][i_s][i_t]

    return s_energies


# ---------------------------------------------------------------------------------

if __name__ == '__main__':
    csvfilename_ccsdt = os.path.join("..", "csvfiles", "AVTZdata.csv")
    csvfilename_smod  = os.path.join("..", "csvfiles", "smodAVTZ.csv")

    csvemanager = CSVEnergyManager()

    ccsdt_energies = csvemanager.read_csv_energies(csvfilename_ccsdt)
    new_energies   = np.empty(ccsdt_energies.shape, dtype = float)

    Rdata = rst.get_Rdata()
    sdata = rst.get_sdata()
    tdata = rst.get_tdata()
    
    with open('sprime_list.txt', 'r') as fin:
        for i_R in np.arange(len(Rdata)):
            for i_t in np.arange(len(tdata)):
                matm_line  = fin.readline()
                ccsdt_engs = ccsdt_energies[i_R,:,i_t]
                matmdata   = MATMReplacementData(matm_line, sdata, ccsdt_engs)
                mod_energies = matmdata.get_modified_energies()
                for i_s in np.arange(len(sdata)):
                    new_energies[i_R][i_s][i_t] = mod_energies[i_s]
    
    csvemanager.write_csv_energies(new_energies, csvfilename_smod)
