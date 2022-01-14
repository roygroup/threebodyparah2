import numpy as np
import threebody_potential_types as tbp
import pes_regions as pr

# ---------------------------------------------------------------------------------

def get_threebody_energies(rstpoint_list, kerneltag):
    rkhspot = tbp.RKHSPotential(kerneltag, 'threebodyparah2')

    # each region needs a diff number of triplets to calculate an energy within it
    regions = get_regions(rstpoint_list)
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
        region.calculate_final_energy()

        lower_cap += region.get_capacity()

        final_energy = region.get_final_energy()
        extrapolated_energies.append(final_energy)

    return extrapolated_energies

def get_regions(rstpoint_list):
    regions = []
    for rstpoint in rstpoint_list:
        region = get_region(rstpoint)
        region.set_rstpoints(rstpoint)
        regions.append(region)
    
    return regions

def get_region(rstpoint):
    assert isinstance(rstpoint, tbp.RstPoint)
    R, s, t = rstpoint.unpack()

    in_Rextrap_region = (R < 2.20) and (1.00 <= s <= 3.85)
    in_sandRextrap_region = (R < 2.20) and (s > 3.85)

    in_sextrap_region = (2.20 <= R < 3.25) and (s > 3.85)
    in_atm1_region = (R >= 3.25) and (s > 3.85)
    in_atm2_region = (R >= 6.25) and (1.00 <= s <= 3.85)

    in_bounded_region = (2.20 <= R < 6.25) and (1.00 <= s <= 3.85)
    
    if   in_atm1_region or in_atm2_region:
        return pr.ATMRegion()
    elif in_bounded_region:
        return pr.BoundedRegion()
    elif in_sextrap_region:
        if t >= 0.00001:
            return pr.SExtrapRegion()
        else:
            return pr.BoundedRegion()
    elif in_Rextrap_region:
        return pr.RExtrapRegion()
    elif in_sandRextrap_region:
        return pr.SandRExtrapRegion()
    else:
        assert False, "unreachable: (R, s, t) = ({}, {}, {})".format(R, s, t)
