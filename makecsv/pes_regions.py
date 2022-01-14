import numpy as np
import threebody_potential_types as tbp
from scipy import integrate

def cos_transition_func(x, xmin, xmax):
    assert xmin < xmax

    if   x <= xmin:
        return 0.0
    elif x >= xmax:
        return 1.0
    else:
        k = (x - xmin)/(xmax - xmin)
        return 0.5*(1.0 - np.cos(np.pi*k))

# ---------------------------------------------------------

class Region:
    def __init__(self, capacity):
        self.rstpoints = []
        self.energies  = []
        self.capacity  = capacity
        self.final_energy = None

    def get_rstpoints(self):
        return self.rstpoints

    def get_capacity(self):
        return self.capacity

    def set_energies(self, energies):
        self.energies.extend(energies)

    def get_final_energy(self):
        return self.final_energy

# ---------------------------------------------------------

class BoundedRegion(Region):
    def __init__(self):
        super().__init__(1)

    def set_rstpoints(self, rstpoint):
        assert isinstance(rstpoint, tbp.RstPoint)
        self.rstpoints = [rstpoint]

    def calculate_final_energy(self):
        e_pes = self.energies[0]
        
        R, s, t = self.rstpoints[0].unpack()
        if R >= 5.35:
            e_atm = tbp.AxilrodTellerMuto().from_RstPoint(self.rstpoints[0])
            w_atm = cos_transition_func(R, 5.35, 6.25)
            w_pes = 1.0 - w_atm
            self.final_energy = w_atm*e_atm + w_pes*e_pes
        elif R > 3.25 and s > 3.5:
            e_atm = tbp.AxilrodTellerMuto().from_RstPoint(self.rstpoints[0])
            w_atm = cos_transition_func(s, 3.5, 3.85)
            w_pes = 1.0 - w_atm
            self.final_energy = w_atm*e_atm + w_pes*e_pes
        else:
            self.final_energy = e_pes

# ---------------------------------------------------------

class ATMRegion(Region):
    def __init__(self):
        super().__init__(0)

    def set_rstpoints(self, rstpoint):
        assert isinstance(rstpoint, tbp.RstPoint)
        self.rstpoints = [rstpoint]

    def calculate_final_energy(self):
        atm = tbp.AxilrodTellerMuto()
        self.final_energy = atm.from_RstPoint(self.rstpoints[0])

    def get_rstpoints(self):
        return []

# ---------------------------------------------------------------------------------

def dfds_nonvanish(s, t):
    cost  = np.cos(t)
    cos2t = np.cos(2.0*t)
    W  = ( s*(cost + np.sqrt(3.0 + cost**2)) )**2
    W2 = W**2

    numer = 3.0 + 2.0*W - (1.0 + 4.0*cost**2)*W2
    denom = (1.0 - 2.0*W*cos2t + W2)**2
    coeff = 6*W/s

    return coeff*numer/denom

def get_a_and_b(R, t, s0, s1, e0, e1):
    rstpoint0 = tbp.RstPoint(R, s0, t)
    rstpoint1 = tbp.RstPoint(R, s1, t)

    atm = tbp.AxilrodTellerMuto()
    u0  = atm.norm(rstpoint0) * e0
    u1  = atm.norm(rstpoint1) * e1
    f0  = atm.fterm(rstpoint0)
    f1  = atm.fterm(rstpoint1)

    integF10, _ = integrate.quad(dfds_nonvanish, s0, s1, args = (t,))
    inv_sint_sq = 1.0 / (np.sin(t)**2)

    a = inv_sint_sq * (u0*f1 - u1*f0)/integF10
    b = inv_sint_sq * (u1 - u0)/integF10

    return a, b

def get_modatm_energy(R, s, t, s0, s1, e0, e1):
    a, b = get_a_and_b(R, t, s0, s1, e0, e1)

    rstp = tbp.RstPoint(R, s, t)
    atm  = tbp.AxilrodTellerMuto()

    fatm = atm.fterm(rstp)
    modatm_energy = (a + b*fatm) * atm.invnorm(rstp)

    return modatm_energy

# ---------------------------------------------------------------------------------

class ShortRangeExtrapolator:
    def __init__(self, e0, e1, R0, R1):
        self.e0 = e0
        self.e1 = e1
        self.R0 = R0
        self.dR = R1 - R0

        self.c_min = 6.0
        self.c_max = 8.0

        self.c_lin = (e1 - e0)/self.dR

    def lin_extrap(self, R):
        return self.e0 + self.c_lin*(R - self.R0)

    def exp_extrap(self, c_exp, R):
        return self.e0 * np.exp(-c_exp*(R - self.R0))
    
    def shift_abs_eng(self, a0):
        a0_min = 1.0e-8
        if a0 < a0_min:
            return a0_min
        else:
            return a0

    def extrap(self, R):
        if np.sign(self.e0) != np.sign(self.e1):
            return self.lin_extrap(R)
        else:
            a0 = self.shift_abs_eng(np.abs(self.e0))
            a1 = self.shift_abs_eng(np.abs(self.e1))
            if a1 >= a0:
                return self.lin_extrap(R)
            else:
                c_exp = -(1.0/self.dR) * np.log(a1/a0)

                if   c_exp <= self.c_min:
                    return self.exp_extrap(c_exp, R)
                elif c_exp >= self.c_max:
                    return self.lin_extrap(R)
                else:
                    w_lin = cos_transition_func(c_exp, self.c_min, self.c_max)
                    e_lin = self.lin_extrap(R)
                    w_exp = 1.0 - w_lin
                    e_exp = self.exp_extrap(c_exp, R)
                    return w_lin*e_lin + w_exp*e_exp

# ---------------------------------------------------------------------------------

class SExtrapRegion(Region):
    def __init__(self):
        super().__init__(2)
        
        # lower and upper values for the s-extrapolation
        self.s_lower = 3.55
        self.s_upper = 3.85
    
    def set_rstpoints(self, rstpoint):
        assert isinstance(rstpoint, tbp.RstPoint)
        
        self.rstpoint0 = rstpoint

        R, s, t = self.rstpoint0.unpack()
        self.rstpoints = [
            tbp.RstPoint(R, self.s_lower, t),
            tbp.RstPoint(R, self.s_upper, t),
        ]

    def calculate_final_energy(self):
        R, s, t = self.rstpoint0.unpack()
        
        s0, s1 = self.s_lower, self.s_upper
        e0, e1 = self.energies[0], self.energies[1]
        e_modatm = get_modatm_energy(R, s, t, s0, s1, e0, e1)

        if R > 3.15:
            atm = tbp.AxilrodTellerMuto()
            e_atm = atm.from_RstPoint(self.rstpoint0)
            w_atm = cos_transition_func(R, 3.15, 3.25)
            w_modatm = 1.0 - w_atm

            self.final_energy = w_atm*e_atm + w_modatm*e_modatm
        else:
            self.final_energy = e_modatm

# ---------------------------------------------------------------------------------

class RExtrapRegion(Region):
    def __init__(self):
        super().__init__(3)

        # the lower and upper values for the R-extrapolation
        self.Rmin = 2.20
        self.dR   = 0.05
        
    def set_rstpoints(self, rstpoint):
        self.rstpoint0 = rstpoint

        R0 = self.Rmin
        R1 = self.Rmin + self.dR
    
        R, s, t = self.rstpoint0.unpack()
        self.rstpoints = [
            tbp.RstPoint(R0, s, t),
            tbp.RstPoint(R1, s, t),
            tbp.RstPoint(R, s, t),
        ]

    def calculate_final_energy(self):
        R, s, t    = self.rstpoint0.unpack()
        e0, e1, e2 = self.energies

        extrapolator = ShortRangeExtrapolator(e0, e1, self.Rmin, self.Rmin + self.dR)
        
        w_pes = cos_transition_func(R, 2.0, 2.2)
        e_pes = e2

        w_exp = 1.0 - w_pes
        e_exp = extrapolator.extrap(R)

        self.final_energy = w_pes*e_pes + w_exp*e_exp

#    def calculate_final_energy_no_mix(self):
#        R, s, t    = self.rstpoint0.unpack()
#        e0, e1, _ = self.energies
#
#        extrapolator      = ShortRangeExtrapolator(e0, e1, self.Rmin, self.Rmin + self.dR)
#        self.final_energy = extrapolator.extrap(R)

# ---------------------------------------------------------------------------------

class SandRExtrapRegion(Region):
    def __init__(self):
        super().__init__(6)
        
        # the lower and upper values for the s-extrapolation
        self.s_lower = 3.55
        self.s_upper = 3.85

        # the lower and upper values for the R-extrapolation
        self.Rmin = 2.20
        self.dR   = 0.05

        self.R0 = self.Rmin
        self.R1 = self.Rmin + self.dR

    def set_rstpoints(self, rstpoint):
        assert isinstance(rstpoint, tbp.RstPoint)

        self.rstpoint0 = rstpoint

        R, s, t = self.rstpoint0.unpack()
        self.rstpoints = [
            tbp.RstPoint(self.R0, self.s_lower, t),
            tbp.RstPoint(self.R0, self.s_upper, t),
            tbp.RstPoint(self.R1, self.s_lower, t),
            tbp.RstPoint(self.R1, self.s_upper, t),
            tbp.RstPoint(R      , self.s_lower, t),
            tbp.RstPoint(R      , self.s_upper, t),
        ]

    def calculate_final_energy(self):
        R, s, t = self.rstpoint0.unpack()
        
        # perform the extrapolation along the s-coordinate
        s0 = self.s_lower
        s1 = self.s_upper

        if t < 1.0e-5:
            t = 1.0e-5

        s_eng0 = get_modatm_energy(self.R0, s, t, s0, s1, self.energies[0], self.energies[1])
        s_eng1 = get_modatm_energy(self.R1, s, t, s0, s1, self.energies[2], self.energies[3])
        s_eng  = get_modatm_energy(R      , s, t, s0, s1, self.energies[4], self.energies[5])
        
        # now perform the extrapolation along the R-coordinate
        extrapolator = ShortRangeExtrapolator(s_eng0, s_eng1, self.R0, self.R1)
        
        w_pes = cos_transition_func(R, 2.0, 2.2)
        e_pes = s_eng

        w_exp = 1.0 - w_pes
        e_exp = extrapolator.extrap(R)

        self.final_energy = w_pes*e_pes + w_exp*e_exp

#    def calculate_final_energy_no_mix(self):
#        R, s, t = self.rstpoint0.unpack()
#
#        if t < 1.0e-5:
#            t = 1.0e-5
#        
#        # perform the extrapolation along the s-coordinate
#        s0 = self.s_lower
#        s1 = self.s_upper
#
#        s_eng0 = get_modatm_energy(self.R0, s, t, s0, s1, self.energies[0], self.energies[1])
#        s_eng1 = get_modatm_energy(self.R1, s, t, s0, s1, self.energies[2], self.energies[3])
#        
#        # now perform the extrapolation along the R-coordinate
#        extrapolator      = ShortRangeExtrapolator(s_eng0, s_eng1, self.R0, self.R1)
#        self.final_energy = extrapolator.extrap(R)
