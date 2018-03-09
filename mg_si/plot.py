import numpy as np
import matplotlib.pyplot as plt

from .reactions import MgSi

def temperature(planet, times, solution, filepath='./', savename='temperatures.png', N_approx=1000):
    Nt = len(times)
    di = int((len(times) - 1) // N_approx)
    N = np.min((Nt // di, (Nt - 1) // di))
    solution = solution[::di][:N]
    times = times[::di][:N]

    T_cmb = solution[:, 0]
    T_um = solution[:, 1]
    # M_c, M_m = planet.reactions.unwrap_Moles(solution[:, 2:], return_sum=True, split_coremantle=True)
    t_plt = times / 3.16e7 / 1e9
    # names_c = planet.params.reactions.core.species
    # names_c.append('core')
    # names_m = planet.params.reactions.mantle.species
    # names_m.append('mantle')

    r_i = np.array([planet.core_layer.r_i(x, recompute=True) for x in T_cmb])

    plt.figure()
    plt.plot(t_plt, T_cmb, label='Tc')
    plt.plot(t_plt, T_um, label='Tm')
    plt.plot(t_plt, r_i / 1e3, label='Ri')
    plt.ylim(0, 7000)
    plt.legend(loc=0)
    plt.ylabel('Temperature (C)')
    plt.grid()
    plt.title('Temps, R_i={:.0f}km'.format(r_i[-1] / 1e3))
    if savename:
        plt.savefig(filepath + savename)

def coremoles(planet, times, solution, filepath='./', savename='coremoles.png', N_approx=1000):
    Nt = len(times)
    di = int((len(times) - 1) // N_approx)
    N = np.min((Nt // di, (Nt - 1) // di))
    solution = solution[::di][:N]
    times = times[::di][:N]

    T_cmb = solution[:, 0]
    T_um = solution[:, 1]
    M_c, M_m = planet.reactions.unwrap_Moles(solution[:, 2:], return_sum=True, split_coremantle=True)
    t_plt = times / 3.16e7 / 1e9
    names_c = planet.params.reactions.core.species
    names_c.append('core')
    names_m = planet.params.reactions.mantle.species
    names_m.append('mantle')

    plt.figure(figsize=(13, 4))
    plt.subplot(121)
    plt.title('Core')
    for i, M in enumerate(M_c[:-1]):
        plt.plot(t_plt, M, label=names_c[i])
    plt.legend(loc=0)
    # plt.ylim(-1,1)
    plt.ylabel('Moles')
    plt.xlabel('time (Byr)')
    plt.grid()
    plt.subplot(122)
    plt.title('Core')
    for i, M in enumerate(M_c):
        M = (M / M[0] - 1) * 100
        plt.plot(t_plt, M, label=names_c[i])
    plt.legend(loc=0)
    # plt.ylim(-2,1)
    plt.ylabel('% change from initial Moles')
    plt.xlabel('time (Byr)')
    plt.grid()
    if savename:
        plt.savefig(filepath + savename)

def composition(planet, times, solution, filepath='./', savename='composition.png', N_approx=1000):
    Nt = len(times)
    di = int((len(times) - 1) // N_approx)
    N = np.min((Nt // di, (Nt - 1) // di))
    solution = solution[::di][:N]
    times = times[::di][:N]

    # T_cmb = solution[:, 0]
    # T_um = solution[:, 1]
    M_c, M_m = planet.reactions.unwrap_Moles(solution[:, 2:], return_sum=True, split_coremantle=True)
    t_plt = times / 3.16e7 / 1e9
    names_c = planet.params.reactions.core.species
    names_c.append('core')
    names_m = planet.params.reactions.mantle.species
    names_m.append('mantle')

    plt.figure(figsize=(13, 4))
    plt.subplot(121)
    plt.title('Core Composition')
    M0 = np.zeros(N)
    for i in [0, 1, 3, 2]:
        if i != 0:
            M0 = M1
        M1 = M0 + M_c[i] / M_c[-1]
        plt.fill_between(t_plt, M0, M1, label=names_c[i])
    plt.legend(loc=0)
    plt.ylim(0, .3)
    plt.ylabel('Mole fraction')
    plt.xlabel('time (Byr)')
    plt.grid()
    plt.subplot(122)
    plt.title('Mantle Layer Composition')
    M0 = np.zeros(N)
    for i in [0, 1, 2, 3, 4]:
        if i != 0:
            M0 = M1
        M1 = M0 + M_m[i] / M_m[-1]
        plt.fill_between(t_plt, M0, M1, label=names_m[i])
    plt.legend(loc=0)
    plt.ylim(0, 1)
    plt.ylabel('Mole fraction')
    plt.xlabel('time (Byr)')
    plt.grid()
    if savename:
        plt.savefig(filepath + savename)

def dTdt(planet, times, solution, filepath='./', savename='dTdt.png', N_approx=1000):
    Nt = len(times)
    di = int((len(times) - 1) // N_approx)
    N = np.min((Nt // di, (Nt - 1) // di))
    solution = solution[::di][:N]
    times = times[::di][:N]

    T_cmb = solution[:,0]
    # T_um = solution[:,1]
    # M_c, M_m = planet.reactions.unwrap_Moles(solution[:, 2:], return_sum=True, split_coremantle=True)
    t_plt = times / 3.16e7 / 1e9
    # names_c = planet.params.reactions.core.species
    # names_c.append('core')
    # names_m = planet.params.reactions.mantle.species
    # names_m.append('mantle')
    plt.figure()

    dTdt = np.diff(T_cmb) / np.diff(times)*3.16e7*1e9
    plt.plot(t_plt[:-1], dTdt)
    plt.grid()
    plt.ylabel('dT/dt (C / Byr)')

    plt.xlabel('time (Byr)')
    plt.grid()
    if savename:
        plt.savefig(filepath + savename)

def MgSiOequilibrium(planet, times, solution, filepath='./', savename='MgSiOeq.png', N_approx=1000):
    Nt = len(times)
    di = int((len(times) - 1) // N_approx)
    N = np.min((Nt // di, (Nt - 1) // di))
    solution = solution[::di][:N]
    times = times[::di][:N]

    # T_cmb = solution[:,0]
    # T_um = solution[:,1]
    # M_c, M_m = planet.reactions.unwrap_Moles(solution[:, 2:], return_sum=True, split_coremantle=True)
    t_plt = times / 3.16e7 / 1e9
    # names_c = planet.params.reactions.core.species
    # names_c.append('core')
    # names_m = planet.params.reactions.mantle.species
    # names_m.append('mantle')
    plt.figure()

    M_Mg_eq = np.zeros(N)
    M_Si_eq = np.zeros(N)
    M_O_eq = np.zeros(N)
    for t in range(N):
        M_Mg_eq[t], M_Si_eq[t], M_O_eq[t] = planet.reactions.compute_Moles_eq(Moles=solution[t, 2:], T_cmb=solution[t, 0])
    plt.plot(t_plt, M_Mg_eq, 'r--', label='M_Mg_eq')
    plt.plot(t_plt, solution[:N, 2], 'r-', label='M_Mg')
    plt.plot(t_plt, M_Si_eq, 'g--', label='M_Si_eq')
    plt.plot(t_plt, solution[:N, 3], 'g-', label='M_Si')
    plt.plot(t_plt, M_O_eq, 'b--', label='M_O_eq')
    plt.plot(t_plt, solution[:N, 5], 'b-', label='M_O')

    # plt.plot(t_plt, solution[:N,4], 'k-', label='M_Fe')
    plt.legend()
    plt.title('equlibrium core/mantle values')
    plt.ylabel('Moles')
    plt.grid()

    plt.xlabel('time (Byr)')
    plt.grid()
    if savename:
        plt.savefig(filepath + savename)

def MgFefraction(planet, times, solution, filepath='./', savename='MgFefraction.png', N_approx=1000):
    Nt = len(times)
    di = int((len(times) - 1) // N_approx)
    N = np.min((Nt // di, (Nt - 1) // di))

    M_c, M_m = planet.reactions.unwrap_Moles(solution[:, 2:], return_sum=True, split_coremantle=True)
    dM_Mg_dt = (np.diff(M_c[0]) / np.diff(times))[::di][:N]
    dM_Fe_dt = (np.diff(M_c[2]) / np.diff(times))[::di][:N]

    times = times[::di][:N]
    solution = solution[::di][:N]

    t_plt = times / 3.16e7 / 1e9
    M_c, M_m = planet.reactions.unwrap_Moles(solution[:, 2:], return_sum=True, split_coremantle=True)
    M_MgO = M_m[0]
    M_SiO2 = M_m[1]
    M_FeO = M_m[2]
    M_MgSiO3 = M_m[3]
    M_FeSiO3 = M_m[4]

    plt.figure()
    plt.title('Mg/(Mg+Fe)')
    plt.plot(t_plt, dM_Mg_dt / (dM_Mg_dt + dM_Fe_dt), label='Mg/Fe exsolve from core')
    plt.plot(t_plt, M_MgO / (M_MgO + M_FeO), label='MgO')
    plt.plot(t_plt, M_MgSiO3 / (M_MgSiO3 + M_FeSiO3), '--', label='MgSiO3')
    plt.legend(loc=0)
    plt.ylabel('Mg/(Mg+Fe)')

    plt.xlabel('time (Byr)')
    plt.grid()
    if savename:
        plt.savefig(filepath + savename)

def K_vals(planet, times, solution, filepath='./', savename='K_vals.png', N_approx=1000):
    Nt = len(times)
    di = int((len(times) - 1) // N_approx)
    N = np.min((Nt // di, (Nt - 1) // di))
    times = times[::di][:N]
    solution = solution[::di][:N]

    M_c, M_m = planet.reactions.unwrap_Moles(solution[:, 2:], return_sum=True, split_coremantle=True)
    t_plt = times / 3.16e7 / 1e9

    plt.figure()

    X_Mg = M_c[0] / M_c[4]
    X_Si = M_c[1] / M_c[4]
    X_Fe = M_c[2] / M_c[4]
    X_O = M_c[3] / M_c[4]

    X_MgO = M_m[0] / M_m[5]
    X_SiO2 = M_m[1]/ M_m[5]
    X_FeO = M_m[2]/ M_m[5]
    X_MgSiO3 = M_m[3] / M_m[5]
    X_FeSiO3 = M_m[4] / M_m[5]

    K1 = X_MgO * X_SiO2 / X_MgSiO3
    K2 = X_FeO * X_SiO2 / X_FeSiO3
    K3 = X_FeO * X_MgSiO3 / (X_MgO * X_FeSiO3)
    K4 = X_Mg * X_O / X_MgO
    K5 = X_Fe * X_O / X_FeO
    K6 = X_Si * X_O ** 2 / X_SiO2
    K4_dat, _ = planet.reactions.func_KD_MgO_val(solution[:, 0])
    K5_dat, _ = planet.reactions.func_KD_FeO_val(solution[:, 0])
    K6_dat, _ = planet.reactions.func_KD_SiO2_val(X_Si, X_O, solution[:, 0])

    plt.figure()
    plt.title('KD Values over time')
    plt.plot(t_plt, K1, '--', label='K_MgSiO3 from X')
    plt.plot(t_plt, K2, '--', label='K_FeSiO3 from X')
    plt.plot(t_plt, K3, '--', label='K_MgFe from X')
    plt.plot(t_plt, K4_dat, '-', label='K_Mg eqn')
    plt.plot(t_plt, K4, '--', label='K_Mg from X')
    plt.plot(t_plt, K5_dat, '-', label='K_Fe eqn')
    plt.plot(t_plt, K5, '--', label='K_Fe from X')
    plt.plot(t_plt, K6_dat, '-', label='K_Si eqn')
    plt.plot(t_plt, K6, '--', label='K_Si from X')
    plt.legend(loc=0)

    plt.xlabel('time (Byr)')
    plt.grid()
    if savename:
        plt.savefig(filepath + savename)

def Q_all(planet, times, all_parameters, filepath='./', savename='Q_all.png'):
    t_plt = times / 3.17e7 / 1e9
    all = all_parameters
    plt.figure(figsize=(10, 7))
    plt.plot(t_plt, all.Qgm / 1e12, 'g-', label="Qg MgO")
    plt.plot(t_plt, all.Qlm / 1e12, 'g--', label="Ql MOg")
    plt.plot(t_plt, all.Qgs / 1e12, 'b-', label="Qg SiO2")
    plt.plot(t_plt, all.Qls / 1e12, 'b--', label="Ql SiO2")
    plt.plot(t_plt, all.Qgf / 1e12, 'r-', label="Qg FeO")
    plt.plot(t_plt, all.Qlf / 1e12, 'r--', label="Ql FeO")
    # plt.plot(t_plt, all.Qrc / 1e12, label='Qr core')
    plt.plot(t_plt, all.Qs / 1e12, label="Qs secular cooling")
    plt.plot(t_plt, all.Qg / 1e12, label="Qg i.core")
    plt.plot(t_plt, all.Ql / 1e12, label="Ql i.core")
    plt.plot(t_plt, all.Qcmb / 1e12, label="Qcmb")
    plt.plot(t_plt, all.Qphi / 1e12, 'k--', label="Qphi")
    plt.plot(t_plt, all.Qphi / 1e12, 'k-', label="Qk")
    plt.legend(loc=0)
    plt.grid()
    plt.ylim(0, 70)
    plt.ylabel('Heat (TW)')
    plt.title('Heat Production')
    plt.xlabel('time (Byr)')
    if savename:
        plt.savefig(filepath + savename)

def E_all(planet, times, all_parameters, filepath='./', savename='E_all.png'):
    t_plt = times / 3.17e7 / 1e9
    all = all_parameters
    min_E = 100e6  # [MW/K]
    plt.figure(figsize=(10, 7))
    plt.plot(t_plt, all.Ephi / 1e6, '--k', label="Entropy available for Dynamo", linewidth=2.)
    plt.plot(t_plt, min_E * np.ones_like(t_plt) / 1e6, '--k', label='Min. Ent. for Dynamo', linewidth=10, alpha=0.5)
    plt.plot(t_plt, all.Es / 1e6, label="Specific Heat (cooling)")
    plt.plot(t_plt, all.Egm / 1e6, label="MgO (grav.)")
    plt.plot(t_plt, all.Egs / 1e6, label="SiO2 (grav.)")
    plt.plot(t_plt, all.Egf / 1e6, label="FeO (grav.)")
    plt.plot(t_plt, all.Eg / 1e6, label="IC (grav.)")
    plt.plot(t_plt, all.El / 1e6, label="IC (latent heat)")
    plt.plot(t_plt, all.Er / 1e6, label='Radioactive heating')
    plt.plot(t_plt, all.Ek / 1e6, label="Adiabatic entropy sink")
    # plt.plot(t_plt, DE/1e6, label="E tot")
    plt.ylim(0, 3000)
    plt.xlim(0, 4.568)
    plt.title('Entropy Production')
    plt.ylabel('Entropy (MW/K)')
    plt.xlabel('time (Byr)')
    plt.legend(loc=0)
    plt.grid()
    # plt.savefig('S_lowvisc.png')
    if savename:
        plt.savefig(filepath + savename)
