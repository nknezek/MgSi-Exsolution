import solubility_library as sol
import numpy as np

s = sol.MgSiSystem()
wtp_i_c = np.array([.01, .01, .85, .13])
mass_c=1.97e24 #kg
wt_i_c = s.wtp_i_c_2_wt_i_c(wtp_i_c, mass_c)
print(wtp_i_c, np.sum(wtp_i_c))
print(wt_i_c, np.sum(wt_i_c))
M_i_c = s.wtp_i_c_2_M_i_c(wtp_i_c, mass_c)
X_i_c = M_i_c/np.sum(M_i_c)
print(M_i_c)
print(X_i_c, np.sum(X_i_c))