import uncertainties as unc
import numpy as np
import utils
import math

u = unc.ufloat(1, 0.1, "u variable")
v = unc.ufloat(1, 0.1, "v variable")
s = u + v

#print u,v,s

#cov_matrix = unc.covariance_matrix([u, v, sum_value])
#print cov_matrix

#corr_matrix = unc.correlation_matrix([s,u,v])
#print corr_matrix

#mm = unc.ufloat(1.25, 1.e-9,)
#at = unc.ufloat(0.51, 1.e-9,)
#qcd = mm-(1.5+(at-0.5)*10.)
#qcd = mm-(1.5+(at-0.5)*10.)
#c_matrix = unc.correlation_matrix([qcd,mm,at])
#print c_matrix

c1 = unc.ufloat(1.5, 1.e-9,)
c2 = unc.ufloat(0.5, 1.e-9,)
c3 = unc.ufloat(10., 1.e-9,)
c4 = unc.ufloat(1.06, 1.e-9,)
at = unc.ufloat(0.51, 1.e-9,)
mm=(c1+(at-c2)*c3)
c_matrix = unc.correlation_matrix([mm,at])
print c_matrix

#gaus = 0.75*math.exp(-0.5*pow((mm.n-(1.5+(at.n-0.5)*5.))/1.06,2.))
#exp = math.exp(70.0-120.0*at.n)*580.
#qcd = gaus*exp
#c_matrix = unc.correlation_matrix([qcd,gaus,exp])
#print c_matrix


#(u2, v2, sum2) = unc.correlated_values([u.n, v.n, sum_value.n], cov_matrix)
#print u2, v2, sum2
#print sum2 - (u2+2*v2)
#print unc.covariance_matrix([u2, v2, sum2])
#
#(u3, v3, sum3) = unc.correlated_values_norm( [u, v, sum_value], corr_matrix)
#print u3, v3, sum3
#print sum3 - (u3+2*v3)
#
#print unc.covariance_matrix([u2, v2, sum2])
#correlation_matrix()
