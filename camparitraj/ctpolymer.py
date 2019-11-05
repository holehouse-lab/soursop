
def get_overlap_concentration(rg):
    """
    Takes Rg in angstroms and returns overlap concentration in Molar

    """
    
    # rg in meters
    Na=6.023e23

    # get rg in meters
    rg = rg*1e-10
        
    # get volume in m3
    v_m3=((4/3)*np.pi*np.power(rg,3))
    
    # get volume in liters
    v_l = v_m3*1000
    
    # get concentration
    c = (1/(v_l))/Na

    return c
