import numpy as np
import os
import pandas as pd
import dill
import sys
sys.path.append('../')
import mg_si.plot as mplt
import matplotlib.pyplot as plt
datafolder = '../computed_solutions_nature/'
datafile = 'data.m'
alldatafile = 'all_parameters.m'
import datetime

column_names = ['time', 'r_i', 'T_cmb0', 'X_Mg_0', 'X_Si_0', 'X_O_0', 'MgNumFp', 'MgNumPv', 'X_MgFeO_b', 'X_SiO2_b', 'nu_present', 'deltaT0', 'layer_thickness', 'overturn']
df10 = pd.read_csv(datafolder+'ri10p_data.csv', names=column_names)
N = len(df10)
for i,row in df10.iterrows():
    try:
        foldername = "Tc{:.1f}_XM{:.3f}_XS{:.3f}_XO{:.3f}/".format(row['T_cmb0'],row['X_Mg_0'],row['X_Si_0'],row['X_O_0'])
        
        time = str(datetime.datetime.now())
        print(time +' '+foldername+', {}/{}'.format(i+1,N))
        if not os.path.exists(datafolder+foldername):
            continue
        if os.path.exists(datafolder+foldername+alldatafile):
            continue
        pl,times,solution = dill.load(open(datafolder+foldername+datafile,'rb'))
        t_N, all_parameters = pl.core_layer.compute_all_parameters(times, solution)
        mplt.Q_all(pl, t_N, all_parameters, filepath=datafolder+foldername)
        mplt.E_all(pl, t_N, all_parameters, filepath=datafolder+foldername)
        dill.dump((t_N,all_parameters), open(datafolder+foldername+alldatafile,'wb'))
        plt.close('all')
    except:
        try:
            print("!!!!!!!!!!!!!!!! something went wrong with:")
            foldername = "Tc{:.1f}_XM{:.3f}_XS{:.3f}_XO{:.3f}/".format(row['T_cmb0'],row['X_Mg_0'],row['X_Si_0'],row['X_O_0'])
            time = str(datetime.datetime.now())
            print(time +' '+foldername+', {}/{}'.format(i+1,N))
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        except:
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! something really wrong")
