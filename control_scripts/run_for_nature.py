import numpy as np
from shutil import copyfile
import matplotlib.pyplot as plt
import sys, os
import dill
sys.path.append('../')
import mg_si
import csv
import datetime

from mg_si import plot as mplt

print(sys.argv)
layer_thickness = 100 # m
overturn = 600 # Myr

times = np.linspace(0,4568e6*365.25*24*3600,20000)

## background mantle state
MgNumFp = 0.8
MgNumPv = 0.93
X_MgFeO_b = 0.311
X_SiO2_b = 0.015
## Mantle viscosity
pl = mg_si.planet.Custom()
nu_present = 10**21/pl.params.mantle.rho #[m^2/s]

T_cmbs_all = np.linspace(4800,6500,round((6500-4800)/100)+1)
Ntc = 3
iT = int(sys.argv[1])

T_cmbs = T_cmbs_all[Ntc*iT:(iT+1)*Ntc]
X_Mgs = np.linspace(1e-5,.05,round(.05/.005)+1)
X_Sis = np.linspace(1e-5,.05,round(.05/.005)+1)
X_Os = np.linspace(1e-5,.15,round(.15/.005)+1)

basefolder = '../computed_solutions_nature/'
for T_cmb0 in T_cmbs:
	for X_Mg_0 in X_Mgs:
		for X_Si_0 in X_Sis:
			for X_O_0 in X_Os:
				print(T_cmb0, X_Si_0, X_Mg_0, X_O_0)
				pl = mg_si.planet.Custom()
				pl.reactions._set_layer_thickness(layer_thickness)
				pl.reactions._set_overturn_time(overturn)
				deltaT0 = pl.mantle_layer.get_dT0(T_cmb0)
				T_um0 = T_cmb0-deltaT0
				try:
					filepath = basefolder+ "Tc{:.1f}_XM{:.3f}_XS{:.3f}_XO{:.3f}/".format(T_cmb0, X_Mg_0, X_Si_0, X_O_0)
					if not os.path.exists(basefolder):
						os.mkdir(basefolder)
					if os.path.exists(filepath+'data.m'):
						print(filepath+' already computed')
						continue
					if not os.path.exists(filepath):
						os.mkdir(filepath)

					Moles_0 = pl.reactions.compute_Moles_0(X_Mg_0, X_Si_0, X_O_0, T_cmb0)
					x0 = [T_cmb0, T_um0]
					x0 = x0+Moles_0
					pl.params.reactions.Moles_0 = Moles_0

					Mm_b = pl.reactions.mantle.compute_Mm_b(X_MgFeO=X_MgFeO_b, X_SiO2=X_SiO2_b, MgNumFp=MgNumFp, MgNumPv=MgNumPv)
					pl.params.reactions.Mm_b = Mm_b

					T_present = 1350 # [K]
					nu_old =  nu_present/1e3
					T_old = T_um0
					A,nu0 = pl.mantle_layer.find_arrenhius_params(nu_present, T_present, nu_old, T_old, set_values=True)

					solution = pl.integrate(times, x0)
					dill.dump((pl,times,solution), open(filepath+'data.m','wb'))
					mplt.temperature(pl, times, solution, filepath=filepath)
					mplt.coremoles(pl, times, solution, filepath=filepath)
					mplt.composition(pl, times, solution, filepath=filepath)
					plt.close('all')
					time = str(datetime.datetime.now())
					r_i = pl.core_layer.r_i(solution[-1,0], one_off=True)
					csvdata = [time, r_i, T_cmb0, X_Mg_0, X_Si_0, X_O_0, MgNumFp, MgNumPv, X_MgFeO_b, X_SiO2_b, nu_present, deltaT0, layer_thickness, overturn]
					with open(basefolder+'run_data{}.csv'.format(iT), 'a') as f:
						writer = csv.writer(f)
						writer.writerow(csvdata)
					f.close()
					del pl
					del csvdata,writer,f
					print('finished computing '+filepath)
				except:
					try:
						del pl
						time = str(datetime.datetime.now())
						r_i = np.nan
						csvdata = [time, r_i, T_cmb0, X_Mg_0, X_Si_0, X_O_0, MgNumFp, MgNumPv, X_MgFeO_b, X_SiO2_b, nu_present, deltaT0, layer_thickness, overturn]
						with open(basefolder+'run_data{}.csv'.format(iT), 'a') as f:
							writer = csv.writer(f)
							writer.writerow(csvdata)
						f.close()
						del csvdata,writer,f
						print('problem with '+str(csvdata))
					except:
						print("couldn't do anything")



