'''
    Some example energy curves using the current PES and the ATM potential.
    Both the `threebodypes.kernel` RKHS kernel and the `threebodyparah2`
    executable must be present in this directory for these plots to work.

    Uncomment and comment the functions in the __main__ section at the bottom
    of this file.
'''

import os
import sys
import numpy as np
import threebody_potential_types as tbp
import matplotlib.pyplot as plt

def interior_angle_to_t(alpha):
    sint = np.sin(alpha) / np.sqrt(1.25 - np.cos(alpha))
    return np.arcsin(sint)

def plot_triangle_as_function_of_R(alpha_in_degrees, Rmin, Rmax):
    s = 1.0
    t = interior_angle_to_t(np.radians(alpha_in_degrees))
    Rdata = np.linspace(Rmin, Rmax)

    # the (R, s, \varphi) points that represent the triangles
    rstpoint_list = [tbp.RstPoint(R, s, t) for R in Rdata]

    # energies from the current PES
    rkhspot = tbp.RKHSPotential('threebodypes.kernel', 'threebodyparah2')
    pes_energies = rkhspot.from_RstPoint_list(rstpoint_list)
    
    # energies from the ATM potential
    atm = tbp.AxilrodTellerMuto()
    atm_energies = atm.from_RstPoint_list(rstpoint_list)

    # plot everything
    fig, ax = plt.subplots()

    ax.set_xlabel(r'$ R \ / \ \mathrm{\AA} $', fontsize = 18)
    ax.set_ylabel(r'$ V_3 \ / \ \mathrm{cm}^{-1} $', fontsize = 18)

    ax.set_xlim(Rdata[0], Rdata[-1])
    ax.axhline(y = 0, color = 'k', ls = ':', alpha = 0.5)

    ax.plot(Rdata, pes_energies, 'C3-')
    ax.plot(Rdata, atm_energies, 'C0--')

    fig.tight_layout()
    plt.show()

def plot_isosceles_triangle(base, height_min, height_max):
    heights = np.linspace(height_min, height_max, 200)
    
    # create list of triangle side lengths
    tsl_list = []
    for h in heights:
        sidelen = np.sqrt(h**2 + 0.25*base**2)
        tsl_list.append(tbp.TriangleSideLengths(base, sidelen, sidelen))
    
    # energies from the current PES
    rkhspot = tbp.RKHSPotential('threebodypes.kernel', 'threebodyparah2')
    pes_energies = rkhspot.from_TriangleSideLengths_list(tsl_list)

    # energies from the ATM potential
    atm = tbp.AxilrodTellerMuto()
    atm_energies = atm.from_TriangleSideLengths_list(tsl_list)

    # plot everything
    fig, ax = plt.subplots()

    ax.set_xlabel(r' height [base = {}]'.format(base) + r' $ \ / \ \mathrm{\AA} $', fontsize = 18)
    ax.set_ylabel(r'$ V_3 \ / \ \mathrm{cm}^{-1} $', fontsize = 18)

    ax.set_xlim(heights[0], heights[-1])
    ax.axhline(y = 0, color = 'k', ls = ':', alpha = 0.5)

    ax.plot(heights, pes_energies, 'C3-')
    ax.plot(heights, atm_energies, 'C0--')

    fig.tight_layout()
    plt.show()

# ---------------------------------------------------------------------------------

if __name__ == '__main__':
    '''
        Plotting the 6 different triangles made by a reference molecule and 2 of
        its 12 nearest neighbours.
    '''
    Rmin = 3.0
    Rmax = 4.5
    #plot_triangle_as_function_of_R(60.0, Rmin, Rmax)
    #plot_triangle_as_function_of_R(90.0, Rmin, Rmax)
    #plot_triangle_as_function_of_R(109.47, Rmin, Rmax)
    #plot_triangle_as_function_of_R(120.0, Rmin, Rmax)
    #plot_triangle_as_function_of_R(146.44, Rmin, Rmax)
    #plot_triangle_as_function_of_R(180.0, Rmin, Rmax)
    
    '''
        Plotting the energy for an isosceles triangle, where two molecules are fixed
        at a distance `base` apart (in Angstroms) and the third molecule moves from
        distance `height_min` to `height_max` (both also in Angstroms).
    '''
    base = 2.2
    height_min = 1.9
    height_max = 8.0
    #plot_isosceles_triangle(base, height_min, height_max)
