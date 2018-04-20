import pandas as pd
import os
import numpy as np

names = np.array(['time', 'r_i', 'T_cmb0', 'X_Mg_0', 'X_Si_0', 'X_O_0', 'MgNumFp', 'MgNumPv', 'X_MgFeO_b', 'X_SiO2_b', 'nu_present', 'deltaT0', 'layer_thickness', 'overturn'])

for name in os.listdir("./"):
    if name[-3:]==".csv":

data = pd.read_csv('run_data_0.csv', names = names)
data['ratio_ri'] = data['r_i']/1220./1e3

tmp1 = np.where(( data['ratio_ri'] <1.1) & ( data['ratio_ri'] > 0.9))

fraction_MgFe_b = 0.8
X_MgFeO_b = 0.16
X_SiO2_b = 0.01
tofile = './valid_results/'
for i in range(0,tmp1[0].shape[0]):
    time, r_i, T_cmb0, X_Mg_0, X_Si_0, X_O_0, MgNumFp, MgNumPv, X_MgFeO_b, X_SiO2_b, nu_present, deltaT0, layer_thickness, overturn = data.loc[tmp1[0][i],['time', 'r_i', 'T_cmb0', 'X_Mg_0', 'X_Si_0', 'X_O_0', 'MgNumFp', 'MgNumPv', 'X_MgFeO_b', 'X_SiO2_b', 'nu_present', 'deltaT0', 'layer_thickness', 'overturn']]
    foldername = "./Tc{:.1f}_XM{:.3f}_XS{:.3f}_XO{:.3f}/".format(T_cmb0, X_Mg_0, X_Si_0, X_O_0)
    os.system('cp -r {fr} {to}'.format(fr=foldername, to=tofile))