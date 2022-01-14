import numpy as np
from scipy.interpolate import interp1d
import subprocess
import os

# -----------------------------------------------------------------------------

class Cartesian3D:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
    
    def dist_from_origin(self):
        return np.sqrt(self.x**2 + self.y**2 + self.z**2)

    def unpack(self):
        return self.x, self.y, self.z

    def __sub__(self, other):
        return Cartesian3D(self.x - other.x, self.y - other.y, self.z - other.z)

    def __str__(self):
        return "(x, y, z) = ({:0.6f}, {:0.6f}, {:0.6f})".format(self.x, self.y, self.z)

# -----------------------------------------------------------------------------

class TriangleSideLengths:
    '''
        Represents three para-H2 molecules in terms of the side lengths of the
        triangle that it composes, with the condition that R12 <= R23 <= R13.
    '''
    def __init__(self, R12, R23, R13):
        self.R12, self.R23, self.R13 = sorted([R12, R23, R13])

    def unpack(self):
        return self.R12, self.R23, self.R13

    def __str__(self):
        return "(R12, R23, R13) = ({:0.6f}, {:0.6f}, {:0.6f})".format(self.R12, self.R23, self.R13)

    def to_RstPoint(self):
        R12_sq = self.R12**2
        R23_sq = self.R23**2
        R13_sq = self.R13**2

        R = self.R12
        r = np.sqrt( 0.5 * (R13_sq + R23_sq - 0.5*R12_sq) )
        
        # floating point imprecision sometimes puts _t outside of [0, 1]; if so, arccos then fails
        _t = (R13_sq - R23_sq) / (2.0*self.R12*r)
        _t = self.adjust_arccos_argument(_t)
        t  = np.arccos(_t)
        
        cost = np.cos(t)
        rmin = 0.5 * R * (cost + np.sqrt(3.0 + cost**2))
        s    = r / rmin

        # floating point imprecision sometimes puts s slightly below 1; if so, construction of RstPoint fails
        if s < 1.0:
            s = 1.0

        return RstPoint(R, s, t)
    
    def adjust_arccos_argument(self, _t):
        if   (_t > 1.0):
            return 1.0
        elif (_t < 0.0):
            return 0.0
        else:
            return _t

    @staticmethod
    def get_arrays_of_pairdists(tsl_list):
        size     = len(tsl_list)
        r12_data = np.empty(size, dtype = float)
        r23_data = np.empty(size, dtype = float)
        r13_data = np.empty(size, dtype = float)
        
        for (i, tsl) in enumerate(tsl_list):
            r12, r23, r13 = tsl.unpack()
            r12_data[i] = r12
            r23_data[i] = r23
            r13_data[i] = r13
        
        return r12_data, r23_data, r13_data

# -----------------------------------------------------------------------------

class RstPoint:
    '''
        Represents three para-H2 molecules in terms of the (R, s, varphi) rescaled
        Jacobi coordinate. The angle `varphi` is represented using `t`.
    '''
    def __init__(self, R, s, t):
        assert R >  0.0, R
        assert s >= 1.0, s
        assert np.pi/2.0 >= t >= 0.0, t

        self.R = R
        self.s = s
        self.t = t

    def unpack(self):
        return self.R, self.s, self.t

    def get_r(self):
        cost = np.cos(self.t)
        return self.s * 0.5 * self.R * (np.sqrt(3.0 + cost**2) + cost)

    def to_TriangleSideLengths(self):
        cost = np.cos(self.t)
        r    = self.get_r()
        
        right_angle_term = r**2 + 0.25 * self.R**2
        correction_term  = self.R * r * cost

        R12 = self.R
        R23 = np.sqrt(right_angle_term - correction_term)
        R13 = np.sqrt(right_angle_term + correction_term)

        return TriangleSideLengths(R12, R23, R13)

    def to_Triplet(self):
        r      = self.get_r()
        point1 = Cartesian3D(-self.R/2.0, 0.0, 0.0)
        point2 = Cartesian3D( self.R/2.0, 0.0, 0.0)
        point3 = Cartesian3D( r*np.cos(self.t), r*np.sin(self.t), 0.0)

        return Triplet(point1, point2, point3)

    def __str__(self):
        return "(R, s, t) = ({:0.6f}, {:0.6f}, {:0.6f})".format(self.R, self.s, self.t)

# -----------------------------------------------------------------------------

class Triplet:
    '''
        Represents three para-H2 molecules in terms of the 3D Cartesian coordinates
        of each of centre of mass.
    '''
    def __init__(self, point1, point2, point3):
        self.point1 = point1
        self.point2 = point2
        self.point3 = point3
        self.eps    = 1.0e-8

    def to_TriangleSideLengths(self):
        dist12 = self.pair_distance(self.point2, self.point1)
        dist13 = self.pair_distance(self.point3, self.point1)
        dist23 = self.pair_distance(self.point3, self.point2)

        r12, r23, r13 = self.distance_correct(dist12, dist13, dist23)

        return TriangleSideLengths(r12, r23, r13)

    def pair_distance(self, p1, p2):
        dp = p1 - p2
        return np.sqrt(dp.x**2 + dp.y**2 + dp.z**2)
    
    def distance_correct(self, r12, r23, r13):
        # in case floating point errors don't follow the triangle side length requirements
        dist12 = r12
        dist23 = r23
        dist13 = r13
        
        if   np.abs(r12 + r13 - r23) < self.eps:
            dist23 -= self.eps
        elif np.abs(r13 + r23 - r12) < self.eps:
            dist12 -= self.eps
        elif np.abs(r23 + r12 - r13) < self.eps:
            dist13 -= self.eps
        
        return dist12, dist23, dist13
    
    @staticmethod
    def get_arrays_of_pairdists(triplets):
        tsl_list = [trip.to_TriangleSideLengths() for trip in triplets]
        return TriangleSideLengths.get_arrays_of_pairdists(tsl_list)

# -----------------------------------------------------------------------------

class RKHSPotential:
    def __init__(self, kerneltag, rkhs_exec):
        self.kerneltag = kerneltag
        self.rkhs_exec = rkhs_exec
        self.tsl_dist_file = "_triplet_distances.dat"
        self.trip_eng_file = "_triplet_energies.dat"
        self.tsl_tmpl      = "{0:0.12f}   {1:0.12f}   {2:0.12f}\n"

        self.cmd = "./{} {} {} > {}".format(self.rkhs_exec,
                                            self.kerneltag,
                                            self.tsl_dist_file,
                                            self.trip_eng_file)

    def from_TriangleSideLengths_list(self, tsl_list):
        r12data, r23data, r13data = TriangleSideLengths.get_arrays_of_pairdists(tsl_list)

        with open(self.tsl_dist_file, 'w') as fout:
            for (r12, r23, r13) in zip(r12data, r23data, r13data):
                fout.write(self.tsl_tmpl.format(r12, r23, r13))
        
        subprocess.call(self.cmd, shell = True)

        energies = np.loadtxt(self.trip_eng_file, usecols = (0,), unpack = True)

        subprocess.call(['rm', self.tsl_dist_file])
        subprocess.call(['rm', self.trip_eng_file])

        return energies

    def from_RstPoint_list(self, rstpoints):
        tsl_list = [rstp.to_TriangleSideLengths() for rstp in rstpoints]
        return self.from_TriangleSideLengths_list(tsl_list)

    def from_Triplet_list(self, triplets):
        tsl_list = [trip.to_TriangleSideLengths() for trip in triplets]
        return self.from_TriangleSideLengths_list(tsl_list)

# -----------------------------------------------------------------------------

class AxilrodTellerMuto:
    C9 = 34336.220013464925
    def __init__(self):
        pass

    def norm(self, rstpoint):
        assert isinstance(rstpoint, RstPoint)

        R, s, t = rstpoint.unpack()
        Tf = self.Tterm(rstpoint)

        return (R**9/64.0) * Tf**(3.0/2.0) / self.C9

    def invnorm(self, rstpoint):
        assert isinstance(rstpoint, RstPoint)

        R, s, t = rstpoint.unpack()
        Tf = self.Tterm(rstpoint)

        return self.C9 * (64.0/R**9) / Tf**(3.0/2.0)

    def Wterm(self, rstpoint):
        assert isinstance(rstpoint, RstPoint)

        R, s, t = rstpoint.unpack()

        cost = np.cos(t)
        W    = ( s*(cost + np.sqrt(3.0 + cost**2)) )**2
        
        return W

    def Tterm(self, rstpoint):
        assert isinstance(rstpoint, RstPoint)

        R, s, t = rstpoint.unpack()

        cos2t = np.cos(2.0*t)
        W     = self.Wterm(rstpoint)

        return 1.0 - 2.0*cos2t*W + W**2

    def fterm(self, rstpoint):
        assert isinstance(rstpoint, RstPoint)

        R, s, t = rstpoint.unpack()

        cost  = np.cos(t)
        cos2t = np.cos(2.0*t)
        W     = self.Wterm(rstpoint)

        Tf    = 1.0 - 2.0*cos2t*W + W**2
        numer = (1.0 - W*cost**2) * (1.0 - W)

        return -3.0 * (numer / Tf)

    def from_RstPoint(self, rstpoint):
        assert isinstance(rstpoint, RstPoint)

        R, s, t = rstpoint.unpack()

        Tterm = self.Tterm(rstpoint)
        fterm = self.fterm(rstpoint)

        return self.C9 * (64.0/R**9) * (1.0 + fterm) / (Tterm**(3.0/2.0))

    def from_RstPoint_list(self, rstpoints):
        return np.array([self.from_RstPoint(rstp) for rstp in rstpoints])

    def from_TriangleSideLengths(self, tsl):
        assert isinstance(tsl, TriangleSideLengths)

        R12, R23, R13 = tsl.unpack()

        denom = (R12*R13*R23)**3

        R12_sq = R12**2
        R23_sq = R23**2
        R13_sq = R13**2
            
        cos1_ = (R12_sq + R13_sq - R23_sq)
        cos2_ = (R12_sq + R23_sq - R13_sq)
        cos3_ = (R13_sq + R23_sq - R12_sq)
        cos_denom = 8.0 * R12_sq * R23_sq * R13_sq
        fterm = 3.0 * cos1_*cos2_*cos3_ / cos_denom

        return self.C9 * (1.0 + fterm) / denom

    def from_TriangleSideLengths_list(self, tsl_list):
        return np.array([self.from_TriangleSideLengths(tsl) for tsl in tsl_list])

    def from_Triplet(self, triplet):
        assert isinstance(triplet, Triplet)
        return self.from_TriangleSideLengths(triplet.to_TriangleSideLengths())

# -----------------------------------------------------------------------------

class ModifiedAxilrodTellerMuto(AxilrodTellerMuto):
    def __init__(self, norm_ccsdt_energies, sdata, t, svalue = None):
        assert len(norm_ccsdt_energies) == len(sdata)

        atm = AxilrodTellerMuto()
        rstpoints = [RstPoint(2.2, s, t) for s in sdata]      # R here is a dummy variable
        
        self.ccsdt = norm_ccsdt_energies
        self.sdata = sdata
        self.fdata = np.array([atm.fterm(rstp) for rstp in rstpoints])

        dvds_ = np.gradient(self.ccsdt, self.sdata)
        dfds_ = np.gradient(self.fdata, self.sdata)

        self.dvds_func = interp1d(sdata, dvds_, kind = 'cubic')
        self.dfds_func = interp1d(sdata, dfds_, kind = 'cubic')
        self.v_func    = interp1d(sdata, self.ccsdt, kind = 'cubic')
        self.f_func    = interp1d(sdata, self.fdata, kind = 'cubic')

        self.a = None
        self.b = None
        if svalue != None:
            self.set_svalue(svalue)

    def set_svalue(self, svalue):
        dvds = self.dvds_func(svalue)
        dfds = self.dfds_func(svalue)
        v    = self.v_func(svalue)
        f    = self.f_func(svalue)

        self.b = dvds/dfds
        self.a = v - self.b*f

    def from_RstPoint(self, rstpoint):
        assert isinstance(rstpoint, RstPoint)
        assert self.a != None
        assert self.b != None

        R, s, t = rstpoint.unpack()

        Tterm = self.Tterm(rstpoint)
        fterm = self.fterm(rstpoint)

        return self.C9 * (64.0/R**9) * (self.a + self.b*fterm) / (Tterm**(3.0/2.0))

    def from_RstPoint_list(self, rstpoints):
        return np.array([self.from_RstPoint(rstp) for rstp in rstpoints])
