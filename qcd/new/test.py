import uncertainties as unc

u = unc.ufloat(1, 0.1, "u variable")
v = unc.ufloat(10, 0.1, "v variable")
sum_value = u+2*v

print u,v,sum_value
print sum_value - (u+2*v)

cov_matrix = unc.covariance_matrix([u, v, sum_value])
print cov_matrix

corr_matrix = unc.correlation_matrix([u, v, sum_value])
print corr_matrix

(u2, v2, sum2) = unc.correlated_values([u.n, v.n, sum_value.n], cov_matrix)
print u2, v2, sum2
print sum2 - (u2+2*v2)
print unc.covariance_matrix([u2, v2, sum2])

(u3, v3, sum3) = unc.correlated_values_norm( [u, v, sum_value], corr_matrix)
print u3, v3, sum3
print sum3 - (u3+2*v3)

#print unc.covariance_matrix([u2, v2, sum2])
#correlation_matrix()
