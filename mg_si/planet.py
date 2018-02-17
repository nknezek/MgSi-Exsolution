import numpy as np
import scipy.integrate as integrate
import mg_si
from mg_si.base import Parameters

class Planet(object):
    def __init__(self, layers, params):
        self.layers = layers
        self.Nlayers = len(layers)

        self.radius = self.layers[-1].outer_radius
        self.volume = 4. / 3. * np.pi * self.radius ** 3
        for layer in self.layers:
            layer.planet = self
        self.params = params

    def integrate(self, x0, xkey=None):
        raise NotImplementedError('must implement an integrate method')

class Stevenson(Planet):
    '''Implements Stevenson 1983 2-layer thermal model

     extends base class Planet
    '''
    def __init__(self, case=1):
        params = Parameters("Stevenson 1983")
        Planet.__init__(self, [mg_si.core.Stevenson(params=params, case=case), mg_si.mantle.Stevenson(params=params)], params)
        self.core_layer = self.layers[0]
        self.mantle_layer = self.layers[1]

    def ODE(self, x, t):
        '''define the ODE for thermal evolution

        :param x:
        :param t:
        :return:
        '''
        T_cmb = x[0]
        T_um = x[1]
        try:
            self.params.Mantle_verbose
            vm = True
        except:
            vm = False
        dTm_dt = self.mantle_layer.energy_balance(T_cmb, T_um, t, verbose=vm)
        cmb_flux = self.mantle_layer.lower_boundary_flux(T_cmb, T_um)
        dTc_dt = self.core_layer.energy_balance(T_cmb, cmb_flux)
        return np.array([dTc_dt, dTm_dt])

    def integrate(self, times, x0, full_output=False):
        solution = integrate.odeint(self.ODE, x0, times, full_output=full_output)
        return solution

class Stevenson_backward(Planet):
    def __init__(self, case=1):
        params = Parameters("Stevenson 1983")
        Planet.__init__(self, [mg_si.core.Stevenson_backwards(params=params, case=case), mg_si.mantle.Stevenson_backwards(params=params)], params)
        self.core_layer = self.layers[0]
        self.mantle_layer = self.layers[1]

    def ODE(self, x, t):
        T_cmb = x[0]
        T_um = x[1]
        dTm_dt = self.mantle_layer.energy_balance(T_cmb, T_um, t)
        cmb_flux = self.mantle_layer.lower_boundary_flux(T_cmb, T_um)
        dTc_dt = self.core_layer.energy_balance(T_cmb, cmb_flux)
        return np.array([-dTc_dt, -dTm_dt])

    def integrate(self, times, x0, full_output=False):
        solution = integrate.odeint(self.ODE, x0, times, full_output=full_output)
        return solution

class NimmoStevenson(Planet):
    def __init__(self, case=1):
        params = Parameters('Stevenson Mantle, Nimmo Core')
        Planet.__init__(self, [mg_si.core.Nimmo(params=params), mg_si.mantle.Stevenson(params=params)], params)
        self.core_layer = self.layers[0]
        self.mantle_layer = self.layers[1]

    def ODE(self, x, t):
        T_cmb = x[0]
        T_um = x[1]
        try:
            self.params.Mantle_verbose
            vm = True
        except:
            vm = False
        dTm_dt = self.mantle_layer.energy_balance(T_cmb, T_um, t, verbose=vm)
        cmb_flux = self.mantle_layer.lower_boundary_flux(T_cmb, T_um)
        dTc_dt = self.core_layer.energy_balance(t, T_cmb, cmb_flux)
        return np.array([dTc_dt, dTm_dt])

    def integrate(self, times, x0, full_output=False):
        solution = integrate.odeint(self.ODE, x0, times, full_output=full_output)
        return solution

class Custom(Planet):
    def __init__(self, case=1):
        params = Parameters('Custom Stevenson Mantle, Custom Nimmo Core')
        Planet.__init__(self, [mg_si.core.Custom(params=params), mg_si.mantle.Custom(params=params)], params)
        self.core_layer = self.layers[0]
        self.mantle_layer = self.layers[1]
        self.reactions = mg_si.reactions.MgSi(params=params)
        self.radiogenics = mg_si.radiogenics.Radiogenics()

    def ODE(self, x, t):
        '''define the ODE for thermal evolution

        :param x:
        :param t:
        :return:
        '''
        T_cmb = x[0]
        T_um = x[1]
        Moles = x[2:]
        try:
            self.params.Mantle_verbose
            vm = True
        except:
            vm = False
        dTm_dt = self.mantle_layer.energy_balance(T_cmb, T_um, t, verbose=vm)
        cmb_flux = self.mantle_layer.lower_boundary_flux(T_cmb, T_um)
        dTc_dt = self.core_layer.energy_balance(t, T_cmb, cmb_flux, Moles)
        dMoles_dt  = self.reactions.dMoles_dt(Moles=Moles, T_cmb=T_cmb, dTc_dt=dTc_dt)
        dM_Mg_dt, dM_Si_dt, dM_Fe_dt, dM_O_dt, dM_c_dt, dM_MgO_dt, dM_SiO2_dt, dM_FeO_dt, dM_MgSiO3_dt, dM_FeSiO3_dt, dM_m_dt = self.reactions.unwrap_Moles(dMoles_dt)

        return np.array([dTc_dt, dTm_dt, dM_Mg_dt, dM_Si_dt, dM_Fe_dt, dM_O_dt, dM_MgO_dt, dM_SiO2_dt, dM_FeO_dt, dM_MgSiO3_dt, dM_FeSiO3_dt])

    def integrate(self, times, x0, full_output=False):
        '''integrate the ODE

        :param times:
        :param x0:
        :param full_output:
        :return:
        '''
        solution = integrate.odeint(self.ODE, x0, times, full_output=full_output)
        return solution
