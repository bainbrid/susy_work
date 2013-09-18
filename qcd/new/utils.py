import numpy as np
import uncertainties as unc
from uncertainties import ufloat
from uncertainties import UFloat
from uncertainties import unumpy
from copy import copy

################################################################################
# UTILITY METHODS ##############################################################
################################################################################

def view(arr) :
    return np.flipud(copy(arr).transpose())

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
            #z[i,j] = ufloat(0.,0.)
            z[i,j] = ufloat(np.nan,np.nan)
#            if abs(n[i,j].n) > 0. : 
#                z[i,j] = ufloat(np.nan,np.nan)
#            else : 
#                err = np.sqrt( n[i,j].s**2. + d[i,j].s**2. )
#                z[i,j] = ufloat(0.,err) 
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

def correlated( input1, input2, operation="+", correlation=1.0 ) :
    corr_matrix = np.array([[1.,correlation],
                            [correlation,1.]])
    (output1,output2) = unc.correlated_values_norm([(input1.n,input1.s),(input2.n,input2.s)], corr_matrix)
    if   operation == "+" : return output1 + output2
    elif operation == "-" : return output1 - output2
    elif operation == "*" : return output1 * output2
    elif operation == "/" : return output1 / output2 if output2 > 0. else unc.ufloat(np.inf,np.inf)
    else : return unc.float(np.nan,np.nan)

def correlated_array( input1, input2, operation="+", correlation=1.0 ) :
    corr_array = np.vectorize( correlated )
    return corr_array( input1, input2, operation, correlation ) 

################################################################################
# TEST #########################################################################
################################################################################

if __name__=="__main__":

    input1 = unc.ufloat(10.,2.)
    input2 = unc.ufloat(10.,2.)
    level = 1.
    
    print "Input:          variable #1 = {:.2f}".format( input1 )
    print "Input:          variable #2 = {:.2f}".format( input2 )
    print "Addition:       uncorrelated= {:.2f}, correlated= {:.2f}".format( input1+input2, correlated(input1,input2,"+",level) )
    print "Subtraction:    uncorrelated= {:.2f}, correlated= {:.2f}".format( input1-input2, correlated(input1,input2,"-",level) )
    print "Multiplication: uncorrelated= {:.2f}, correlated= {:.2f}".format( input1*input2, correlated(input1,input2,"*",level) )
    print "Division:       uncorrelated= {:.2f}, correlated= {:.2f}".format( input1/input2, correlated(input1,input2,"/",level) )
    
    in1 = np.empty((1,5),dtype=UFloat)
    in2 = np.empty((1,5),dtype=UFloat)
    in1.fill(input1)
    in2.fill(input2)
    print "Input array #1, input array #2, subtraction:"
    print in1
    print in2
    print correlated_array(in1,in2,"-",level)
