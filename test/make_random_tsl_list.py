import numpy as np
import threebody_potential_types as tbp

class RandomTriangleGenerator:
    def __init__(self, Rmin, Rmax, smax):
        self.Rmin = Rmin
        self.Rmax = Rmax
        self.smin = 1.0
        self.smax = smax
        self.tmin = 0.0
        self.tmax = np.pi/2.0
    
    def get_triangle(self):
        R = np.random.uniform(self.Rmin, self.Rmax)
        s = np.random.uniform(self.smin, self.smax)
        t = np.random.uniform(self.tmin, self.tmax)

        tsl = tbp.RstPoint(R, s, t).to_TriangleSideLengths()

        return tsl


if __name__ == '__main__':
    line_tmpl = '{0: 0.12f}   {1: 0.12f}   {2: 0.12f}\n'
    filename  = 'tsl_random.dat'

    trianglegen = RandomTriangleGenerator(1.9, 8.0, 5.0)
    Ntriangles  = 256
    
    tsl_list = [trianglegen.get_triangle() for i in range(Ntriangles)]

    with open(filename, 'w') as fout:
        for tsl in tsl_list:
            fout.write(line_tmpl.format(*tsl.unpack()))
