import numpy as np
from uncertainties import ufloat
from uncertainties import UFloat
from uncertainties import unumpy

################################################################################
# UTILITY METHODS ##############################################################
################################################################################

def safe_div(n,d) :
    if n.size == 0 : return None
    if n.shape != d.shape : return None
    if isinstance(n[0,0],UFloat) is False : return None
    z = np.zeros_like(n)
    z.fill(ufloat(np.nan,np.nan))
    it = np.nditer(d,op_flags=['readonly'],flags=['multi_index','refs_ok']) 
    while not it.finished :
        i = it.multi_index[0]
        j = it.multi_index[1]
        if abs(d[i,j].n) == 0. :
            if abs(n[i,j].n) > 0. : 
                z[i,j] = ufloat(np.nan,np.nan)
            else : 
                err = np.sqrt( n[i,j].s**2. + d[i,j].s**2. )
                z[i,j] = ufloat(0.,err) 
                # z[i,j] = ufloat(np.nan,np.nan) 
                # z[i,j] = ufloat(np.inf,np.inf)
        else :
            z[i,j] = n[i,j] / d[i,j]
        it.iternext()
    return z

def safe_div_(n,d) :
    z = ufloat(np.nan,np.nan)
    if abs(d.n) == 0. :
        if abs(n.n) > 0. : 
            z = ufloat(np.nan,np.nan)
        else : 
            err = np.sqrt( n.s**2. + d.s**2. )
            z = ufloat(0.,err) 
    else :
        z = n / d
    return z
