import pandas as pd
import os

names = np.array(['time', 'T_cmb0', 'T_um0', 'r_i', 'wpms','wpss','wpos', 'fraction_MgFe_b', 'X_MgFeO_b','X_SiO2_b', 'nu_present','layer_thickness','overturn','X_Mg_0', 'X_Si_0', 'X_O_0'])
data = pd.read_csv('run_data_200.csv', names = names)
data['ratio_ri'] = data['r_i']/1290./1e3

tmp1 = np.where(( data['ratio_ri'] <1.1) & ( data['ratio_ri'] > 0.9)) 

fraction_MgFe_b = 0.8
X_MgFeO_b = 0.16
X_SiO2_b = 0.01
tofile = './valid_results/'
for i in range(0,tmp1[0].shape[0]):
	T_cmb0,wpms,wpss,wpos, fraction_MgFe_b, X_MgFeO_b, X_SiO2_b, nu_present,layer_thickness,overturn = data.loc[tmp1[0][i],['T_cmb0', 'wpms', 'wpss', 'wpos', 'fraction_MgFe_b','X_MgFeO_b','X_SiO2_b','nu_present','layer_thickness','overturn']]
	foldername = "./Tc{Tc:d}_WtMg{wtmg:.3f}_WtSi{wtsi:.3f}_WtO{wto:.3f}_fMb{fmb:.2f}_Xmb{xmb:.2f}_XSb{xsb:.2f}_nu{nu:.2e}_lthck{lthck:.0e}_ovt{ovt:.0e}".format(Tc=T_cmb0, wtmg=wpms,wtsi=wpss,wto=wpos, fmb=fraction_MgFe_b, xmb=X_MgFeO_b, xsb=X_SiO2_b, nu=nu_present,lthck=layer_thickness,ovt=overturn)
	os.system('cp -r {fr} {to}'.format(fr=foldername, to=tofile))