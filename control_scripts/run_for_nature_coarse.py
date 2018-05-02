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
r_i_real = 1220e3


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

T_min_all = 4900
T_max_all = 6000
T_cmbs_all = np.linspace(T_min_all,T_max_all,round((T_max_all-T_min_all)/100)+1)
Ntc = 3
iT = int(sys.argv[1])

T_cmbs = T_cmbs_all[Ntc*iT:(iT+1)*Ntc]
dMg = .01
minMg = 1e-5
maxMg = .05

dO = .01
minO = 1e-5
maxO = .20

dSi = .01
minSi = 1e-5
maxSi = .15

X_Mgs = np.linspace(minMg,maxMg, round((maxMg-minMg)/dMg)+1)
X_Sis = np.linspace(minSi, maxSi, round((maxSi-minSi)/dSi)+1)
X_Os = np.linspace(minO, maxO, round((maxO-minO)/dO)+1)

basefolder = '/media/nknezek/compute_storage/computed_solutions_nature/'
alldatafile = 'all_parameters.m'

Ntotal = len(T_cmbs)*len(X_Mgs)*len(X_Sis)*len(X_Os)
i = 1
for T_cmb0 in T_cmbs:
	for X_Mg_0 in X_Mgs:
		for X_Si_0 in X_Sis:
			for X_O_0 in X_Os:
				time = str(datetime.datetime.now())
				print('{} - {:.0f}K - {:.3f}Mg {:.3f}Si {:.3f}O - {}/{}'.format(time, T_cmb0, X_Mg_0, X_Si_0, X_O_0,i,Ntotal))
				i += 1
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
						print('already computed')
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

					# plot and store solution info
					solution = pl.integrate(times, x0)
					mplt.temperature(pl, times, solution, filepath=filepath)
					mplt.coremoles(pl, times, solution, filepath=filepath)
					mplt.composition(pl, times, solution, filepath=filepath)
					plt.close('all')
					dill.dump((pl,times,solution), open(filepath+'data.m','wb'))
					r_i = pl.core_layer.r_i(solution[-1,0], one_off=True)
					
					# Store Run Info into csv file
					csvdata = [time, r_i, T_cmb0, X_Mg_0, X_Si_0, X_O_0, MgNumFp, MgNumPv, X_MgFeO_b, X_SiO2_b, nu_present, deltaT0, layer_thickness, overturn]
					with open(basefolder+'run_data{}.csv'.format(iT), 'a') as f:
						writer = csv.writer(f)
						writer.writerow(csvdata)
					f.close()

					# if the inner-core size is within 10% of real inner-core, compute entropy and heat terms
					if np.abs(r_i/r_i_real-1)<.1:
						print('r_i within 10%, computing and storing entropy history')
						t_N, all_parameters = pl.core_layer.compute_all_parameters(times, solution)
						mplt.Q_all(pl, t_N, all_parameters, filepath=filepath)
						mplt.E_all(pl, t_N, all_parameters, filepath=filepath)
						dill.dump((t_N,all_parameters), open(filepath+alldatafile,'wb'))
						plt.close('all')
					del pl
					del csvdata
					print('==== successfully finished computing')
				except:
					try:
						del pl
						r_i = 'nan'
						csvdata = [time, r_i, T_cmb0, X_Mg_0, X_Si_0, X_O_0, MgNumFp, MgNumPv, X_MgFeO_b, X_SiO2_b, nu_present, deltaT0, layer_thickness, overturn]
						with open(basefolder+'run_data{}.csv'.format(iT), 'a') as f:
							writer = csv.writer(f)
							writer.writerow(csvdata)
						f.close()
						print('############## problem with '+str(csvdata)+'\n')
					except:
						print("!!!!!!!!!!!!!!!!!!!!!!!\ncouldn't do anything\n!!!!!!!!!!!!!!!!!\n")



