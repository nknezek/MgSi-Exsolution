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

import signal

class timeout:
    def __init__(self, seconds=1, error_message=None):
        self.seconds = seconds
        if error_message is None:
            self.error_message = 'solution timed out after {}s'.format(seconds)
    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)
    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)
    def __exit__(self, type, value, traceback):
        signal.alarm(0)

r_i_real = 1220e3

layer_thickness = 100 # m
overturn = 1e10 # Myr
times = np.linspace(0,4568e6*365.25*24*3600,40000)

## background mantle state
MgNumFp = 0.8
MgNumPv = 0.93
X_MgFeO_b = 0.311
X_SiO2_b = 0.015
## Mantle viscosity
pl = mg_si.planet.Custom()
nu_present = 10**21/pl.params.mantle.rho #[m^2/s]


T_cmb0 = 5700.0

X_Mgs = [0.010]
X_Sis = [0.120]
X_Os = [0.080]

basefolder = '/Users/nknezek/code/MgSi-Exsolution/comparison_runs/inf_mantle_overturn/'
alldatafile = 'new_parameters.m'

Ntotal = len(X_Mgs)*len(X_Sis)*len(X_Os)
i = 1
for X_Mg_0 in X_Mgs:
    for X_Si_0 in X_Sis:
        for X_O_0 in X_Os:    
            times = np.linspace(0,4568e6*365.25*24*3600,40000)

            time = str(datetime.datetime.now())
            print('{} - {:.0f}K - {:.3f}Mg {:.3f}Si {:.3f}O'.format(time, T_cmb0, X_Mg_0, X_Si_0, X_O_0))
            pl = mg_si.planet.Custom()
            pl.reactions._set_layer_thickness(layer_thickness)
            pl.reactions._set_overturn_time(overturn)
            deltaT0 = pl.mantle_layer.get_dT0(T_cmb0)
            T_um0 = T_cmb0-deltaT0
            try:
                filepath = basefolder+ "Tc{:.1f}_XM{:.3f}_XS{:.3f}_XO{:.3f}/".format(T_cmb0, X_Mg_0, X_Si_0, X_O_0)
                if not os.path.exists(basefolder):
                    os.mkdir(basefolder)
                if os.path.exists(filepath+'r_i.m'):
                    print('r_i already computed')
                    r_i = dill.load(open(filepath+'r_i.m','rb'))
                    continue
                if os.path.exists(filepath+'data.m'):
                    print('already computed')
                    try:
                        pl,times,solution = dill.load(open(filepath+'data.m','rb'))
                        r_i = pl.core_layer.r_i(solution[-1,0], one_off=True)
                        dill.dump(r_i, open(filepath+'r_i.m','wb'))
                        del pl
                        continue
                    except:
                        print('could not load previously computed solution, re-computing now')
                if not os.path.exists(filepath):
                    os.mkdir(filepath)
            except:
                print('!!!!! Problem setting up folders',sys.exc_info()[1])
                del pl
                continue
            try:
                pl.params.reactions.ParamCitationFeO = 'Fischer2015'
                pl.params.reactions.ParamCitationSiO2 = 'Fischer2015'
                pl.params.reactions.ParamCitationMgO = 'Badro2015'

                Moles_0 = pl.reactions.compute_Moles_0(X_Mg_0, X_Si_0, X_O_0, T_cmb0)
                x0 = [T_cmb0, T_um0]
                x0 = x0+Moles_0
                pl.params.reactions.Moles_0 = Moles_0
                Mm_b = pl.reactions.mantle.compute_Mm_b(X_MgFeO=X_MgFeO_b, X_SiO2=X_SiO2_b, MgNumFp=MgNumFp, MgNumPv=MgNumPv)
                pl.params.reactions.Mm_b = Mm_b
            except:
                print('Problem with initial mole 0',sys.exc_info()[1])
                del pl
                continue
            try:
                T_present = 1350 # [K]
                nu_old =  nu_present/1e3
                T_old = T_um0
                A,nu0 = pl.mantle_layer.find_arrenhius_params(nu_present, T_present, nu_old, T_old, set_values=True)
            except:
                print('!!!!! Problem setting viscosity',sys.exc_info()[1])
                del pl
                continue
            # plot and store solution info
            timeoutsecs = 120
            try:
                with timeout(seconds=timeoutsecs):
                    solution = pl.integrate(times, x0)  
            except TimeoutError as err:
                print('!!!!!!',err,'with {} timesteps'.format(len(times)))
                continue
            except:
                print("!!!!!! Unexpected error during solution:",sys.exc_info()[1])
                continue
            try:
                r_i = pl.core_layer.r_i(solution[-1,0], one_off=True)
            except:
                print('!!!!! Problem computing r_i',sys.exc_info()[1])
                del pl
                continue
            try:
                dill.dump(r_i, open(filepath+'r_i.m','wb'))
                dill.dump((pl, times, solution), open(filepath+'data.m','wb'))
            except:
                print('!!!!! Problem saving solution',sys.exc_info()[1])
                del pl
                continue
            try:
                # if the inner-core size is within 10% of real inner-core, compute entropy and heat terms
                if np.abs(r_i/r_i_real-1)<.1:
                    print('r_i {:.0f} km within 10%, computing and storing entropy history'.format(r_i/1e3))
                    t_N, all_parameters = pl.core_layer.compute_all_parameters(times, solution)
            #             mplt.Q_all(pl, t_N, all_parameters, filepath=filepath)
                    mplt.E_all(pl, t_N, all_parameters, filepath=filepath)
                    dill.dump((t_N, all_parameters), open(filepath+alldatafile,'wb'))
                    plt.close('all')
                else:
                    print('r_i {:.0f} km not within 10%, entropies not computed'.format(r_i/1e3))
                del pl
                print('==== successfully finished computing')
                continue
            except:
                print('!!!!! Problem computing entropies',sys.exc_info()[1])
                continue