import numpy as np
from .base import Parameters, Layer


class MantleLayer(Layer):
    def __init__(self, inner_radius, outer_radius, params=None):
        Layer.__init__(self, inner_radius, outer_radius, params)

class Stevenson(MantleLayer):
    def __init__(self, params=None):
        if params is None:
            params = Parameters('Stevenson 1983 for mantle')
        self.params = params
        self.params.mantle = Parameters('Stevenson 1983 for mantle')
        pm = self.params.mantle
        pm.R_p0 = 6371e3  # - [m] from Stevenson 1983 Table II
        pm.R_c0 = 3485e3  # - [m] from Stevenson 1983 pg. 474
        pm.T_s = 293.  # - [K] from Stevenson 1983 Table II
        pm.mu = 1.3 # - [] from Stevenson 1983 pg. 473 and Table II
        pm.alpha = 2e-5 # - [/K] from Stevenson 1983 Table I
        pm.k = 4.0 # - [W/m-K] from Stevenson 1983 Table I
        pm.K = 1e-6 # - [m^2/s] from Stevenson 1983 Table I
        pm.rhoC = 4e6 # - [J/m^3-K] from Stevenson 1983 Table I
        pm.rho = 5000. # - [kg/m^3] -- guess as Stevenson 1983 never explicitly states his assumption for rho or C
        pm.C = pm.rhoC/pm.rho # - [J/K-kg]
        pm.Q_0 = 1.7e-7 # - [W/m^3] from Stevenson 1983 Table I
        pm.lam = 1.38e-17 # - [1/s] from Stevenson 1983 Table I
        pm.A = 5.2e4 # - [K] from Stevenson 1983 Table I
        pm.nu_0 = 4.0e3 # - [m^2/s] from Stevenson 1983 Table I
        pm.Ra_crit = 5e2 # - [] from Stevenson 1983 Table I
        pm.beta = 0.3 # - [] from Stevenson 1983 Table I
        pm.g = 10. # - [m/s^2] from Stevenson 1983 Table II
        pm.Ra_boundary_crit = 2e3 # empirical parameparams

        MantleLayer.__init__(self, pm.R_c0, pm.R_p0, params)

    def adiabat_from_bottom(self, T_magma_ocean_top, distance):
        pm = self.params.mantle
        return T_magma_ocean_top * (1.0 - pm.alpha * pm.g * distance / pm.C)

    def average_mantle_temp(self, T_upper_mantle):
        pm = self.params.mantle
        return T_upper_mantle * pm.mu

    def kinematic_viscosity(self, T_upper_mantle):
        pm = self.params.mantle
        return pm.nu_0 * np.exp(pm.A / T_upper_mantle)

    def find_arrenhius_params(self, nu1, T1, nu2, T2, set_values=True):
        '''
        given a set of two viscosities and temperatures, returns the parameters for Arrehnius visocity relation

        nu(T) = nu0*exp(A/T)

        :param nu1:
        :param T1:
        :param nu2:
        :param T2:
        :return: A, nu0
        '''

        A = np.log(nu1 / nu2) / (1 / T1 - 1 / T2)
        nu0 = nu1*np.exp(-A/T1)
        if set_values:
            pr = self.params.mantle
            pr.A = A
            pr.nu_0 = nu0
        return A, nu0

    def heat_production(self, time):
        '''
        Equation (2) from Stevenson et al 1983
        '''
        pm = self.params.mantle
        return pm.Q_0 * np.exp(-pm.lam * time)

    def lower_mantle_temperature(self, T_upper_mantle):
        '''
        Adiabatic Temperature Increase from the temperature at the base of upper mantle boundary layer to
        the top of the lower boundary layer assuming negligable boundary layer thickness.
        '''
        pm = self.params.mantle
        return T_upper_mantle * (1.0 + pm.alpha * pm.g * self.thickness / pm.C)

    def mantle_rayleigh_number(self, T_mantle_bottom, T_upper_mantle):
        '''
        Equation (19) Stevenson et al 1983
        '''
        pm = self.params.mantle
        nu = self.kinematic_viscosity(T_upper_mantle)
        T_lower_mantle = self.lower_mantle_temperature(T_upper_mantle)
        upper_boundary_delta_T = T_upper_mantle - pm.T_s
        lower_boundary_delta_T = T_mantle_bottom - T_lower_mantle
        # assert upper_boundary_delta_T > 0.0
        # assert lower_boundary_delta_T > 0.0
        delta_T_effective = upper_boundary_delta_T + lower_boundary_delta_T
        return pm.g * pm.alpha * (delta_T_effective) * np.power(self.thickness, 3.) / (nu * pm.K)

    def boundary_layer_thickness(self, Ra):
        '''
        Equation (18) Stevenson et al 1983
        '''
        pm = self.params.mantle
        return self.thickness * np.power(pm.Ra_crit / Ra, pm.beta)

    def upper_boundary_layer_thickness(self, T_mantle_bottom, T_upper_mantle):
        '''
        Use Equations (18,19) from Stevenson et al 1983
        '''
        Ra = self.mantle_rayleigh_number(T_mantle_bottom, T_upper_mantle)
        return self.boundary_layer_thickness(Ra)

    def lower_boundary_layer_thickness(self, T_mantle_bottom, T_upper_mantle):
        '''
        Equations (20,21) Stevenson et al 1983

        '''
        pm = self.params.mantle
        T_lower_mantle = self.lower_mantle_temperature(T_upper_mantle)
        delta_T_lower_boundary_layer = T_mantle_bottom - T_lower_mantle
        average_boundary_layer_temp = T_lower_mantle + delta_T_lower_boundary_layer / 2.
        nu_crit = self.kinematic_viscosity(T_upper_mantle)
        # assert delta_T_lower_boundary_layer > 0.0, "dTlbl={3:.1f} K, T_mb={0:.1f} K, T_lm={1:.1f} K, T_um={2:.1f} K".format(
        #     T_mantle_bottom, T_lower_mantle, T_upper_mantle, delta_T_lower_boundary_layer)
        delta_c = np.power(pm.Ra_boundary_crit * nu_crit * pm.K / (pm.g * pm.alpha * delta_T_lower_boundary_layer),
                           0.333)
        Ra_mantle = self.mantle_rayleigh_number(T_mantle_bottom, T_upper_mantle)
        delta_c_normal = self.boundary_layer_thickness(Ra_mantle)
        return np.minimum(delta_c, delta_c_normal)

    def upper_boundary_flux(self, T_mantle_bottom, T_upper_mantle):
        '''
        Equation (17) from Stevenson et al 1983

        :param T_upper_mantle:
        :param T_mantle_bottom:
        :return:
        '''
        pm = self.params.mantle
        delta_T = T_upper_mantle - pm.T_s
        upper_boundary_layer_thickness = self.upper_boundary_layer_thickness(T_mantle_bottom, T_upper_mantle)
        return pm.k * delta_T / upper_boundary_layer_thickness

    def lower_boundary_flux(self, T_mantle_bottom, T_upper_mantle):
        '''
        Equation (17) from Stevenson et al 1983

        :param T_upper_mantle:
        :param T_mantle_bottom:
        :return:
        '''
        pm = self.params.mantle
        delta_T = T_mantle_bottom - self.lower_mantle_temperature(T_upper_mantle)
        lower_boundary_layer_thickness = self.lower_boundary_layer_thickness(T_mantle_bottom, T_upper_mantle)
        return pm.k * delta_T / lower_boundary_layer_thickness

    def energy_balance(self, T_mantle_bottom, T_upper_mantle, time, verbose=False):
        if verbose:
            print(time/3e7/1e6, T_mantle_bottom, T_upper_mantle)
        pm = self.params.mantle
        mantle_surface_area = self.outer_surface_area
        core_surface_area = self.inner_surface_area

        effective_heat_capacity =  pm.rho * pm.C * pm.mu * self.volume
        internal_heat_energy = self.heat_production(time) * self.volume
        cmb_flux = self.lower_boundary_flux(T_mantle_bottom, T_upper_mantle)
        surface_flux = self.upper_boundary_flux(T_mantle_bottom, T_upper_mantle)
        net_flux_out = mantle_surface_area * surface_flux - core_surface_area * cmb_flux
        dTdt = (internal_heat_energy - net_flux_out) / effective_heat_capacity
        return dTdt

class Stevenson_backwards(Stevenson):
    def __init__(self, params=None):
        Stevenson.__init__(self, params)

    def heat_production(self, time):
        '''
        Equation (2) from Stevenson et al 1983
        '''
        pm = self.params.mantle
        return pm.Q_0 * np.exp(-pm.lam * (pm.time_end-time))

class Driscoll(MantleLayer):
    pass

class Korenaga2005(MantleLayer):
    def __init__(self):
        def __init__(self, params=None):
            if params is None:
                params = Parameters('Korenaga 2005 for mantle')
            self.params = params
            self.params.mantle = Parameters('Korenaga 2005 for mantle')
            pm = self.params.mantle
            pm.R_c0 = 3480e3  # [m] core radius
            pm.R_p0 = 6372e3  # [m] planet radius
            pm.C = 7e27  # [J/K] heat capacity of whole earth (Korenaga 2005, below eq (1))
            pm.Up = 0.3  # Urey ratio at present day (Korenaga 2005)
            pm.E = 300e3  # [J/mol] activation enthalpy in lower mantle
            pm.T_off = 273  # [K] temperature offset for mantle from surface
            pm.R = 8.314  # [J/mol-K] universal gas constant
            pm.Q0 = 36e12  # [W] mantle heat flow at reference mantle temperature T0
            pm.T0 = 1350  # [K] mantle reference temperature at present day
            pm.beta = 0.3  # [-] exponent for heat flow calculations Nu ~ Ra^beta

            # pm.rho_0 =
            # pm.g =
            # pm.D =
            # pm.h =
            # pm.nu_M =
            # pm.nu_L =
            pm.RadC = 200e3  # [m] radius of curvature of plate bending

            pm.CGyr2s = 1e9 * 365.25 * 24 * 3600
            pm.l_238U = 4.9160e-18  # [1/Gyr > 1/s]
            pm.l_235U = 3.1209e-17  # [1/Gyr > 1/s]
            pm.l_232Th = 1.5633e-18  # [1/Gyr > 1/s]
            pm.l_40K = 1.7558e-17  # [1/Gyr > 1/s]

            MantleLayer.__init__(self, pm.R_c0, pm.R_p0, params)

    def nu(self, T):
        '''
        Calculate non-normalized mantle viscosity eq (10) Korenaga 2005

        :param T:
        :param E:
        :return:
        '''
        pm = self.params.mantle
        return np.exp(pm.E / (pm.R * (T + pm.T_off)))

    def Q_traditional(self, T):
        '''
        Calculates mantle heat flow Q [W] from mantle tempearture T [K] eq (9) Korenaga 2005

        :param T:
        :return:
        '''
        C = C0(pm.Q0, pm.T0, pm.E)
        return C * T * (T / nu(T)) ** beta

    def find_Q0p(self):
        '''
        finds the constant 'a' in eq (19) Korenaga 2005

        :return:
        '''
        pm = self.params.mantle
        Q = self.Q_plates(pm.T0)
        pm.a = pm.Q0 / Q
        return pm.a

    def Q_plates(self, T):
        '''
        Calculates surface heat flow scaling with temperature from plate tectonic derviation, eq (19) Korenaga 2005

        :param T:
        :return:
        '''
        Q = pm.a * (pm.C_pe * pm.alpha * pm.rho_0 * pm.g * pm.D * pm.h * T ** 3
                    / (pm.CM_vd * pm.nu_M + pm.CS_vd * pm.nu_L * (pm.h / pm.RadC) ** 3)) ** 0.5
        return Q

    def C0(self):
        '''
        Calculates the constant relating mantle temperature to heat flow Q, eq (9, 10) Korenaga 2005

        :return:
        '''
        pm = self.params.mantle
        C0 = pm.Q0 * pm.T0 ** (-1 - pm.beta) * np.exp(pm.beta * pm.E / (pm.R * (pm.T0 + pm.Toff)))
        return C0

    def heat_production_mantle(self, Hp, t, tp=1.44155e17):
        '''
        Calculates heat production through time given present-day heat production Hp

        :param Hp: heat-production at current day [W/kg, TW, or W/m^3]
        :param t: time [s]
        :param tp: end-time [s] default: 1.44155e17 s (4.568 Byr)
        :return: Hm(t): heat-production at time t, in same units as provided in Hp
        '''
        h_238U = 0.37196  # normalized W/kg Koreneaga 2005
        h_235U = 0.01638  # normalized W/kg Koreneaga 2005
        h_232Th = 0.43028  # normalized W/kg Koreneaga 2005
        h_40K = 0.18137  # normalized W/kg Koreneaga 2005
        pm = self.params.mantle

        return Hp * (self._heat_hp(h_238U, pm.l_238U, t, tp=tp)
                     + self._heat_hp(h_235U, pm.l_235U, t, tp=tp)
                     + self._heat_hp(h_232Th, pm.l_232Th, t, tp=tp)
                     + self._heat_hp(h_40K, pm.l_40K, t, tp=tp))

    def heat_production_core(self, Hp, t, tp=1.44155e17):
        '''
        Calculates heat production through time given present-day heat production

        :param Hp: heat production at present day in [W, W/kg, W/m^3]
        :param t: time [s]
        :param tp: present-day time [s] default 1.44155e17 s (4.568 Byr)
        :return: Hc(t) heat production in core same units as provided in Hp
        '''
        pm = self.params.mantle
        h_40K = 1.  # normalized W/kg
        return Hp * self._heat_hp(h_40K, pm.l_40K, t, tp=tp)

    def _heat_h0(self, h0, lam, t):
        '''
        Calculates heat production through time from heat at time 0

        :param h0:
        :param lam:
        :param t:
        :return:
        '''
        return h0 * np.exp(-lam * t)

    def _heat_hp(self, hp, lam, t, tp=1.44155e17):
        '''
        Calculates heat production through time from heat at present time

        :param hp:
        :param lam:
        :param t:
        :param tp:
        :return:
        '''
        return hp * np.exp(lam * (tp - t))

    def energy_balance(self, time, T_mantle, ):
        pass

class ORourke2016(MantleLayer):
    def __init__(self, params=None):
        if params is None:
            params = Parameters('ORourke 2016 for mantle')
        self.params = params
        self.params.mantle = Parameters('ORourke 2016 for mantle')
        pm = self.params.mantle
        pm.Heff = 300e3  # [J/mol] activation enthalpy in lower mantle
        pm.T_ref = 2500  # [K] reference temperature for lower mantle
        pm.nu_ref = 2e21  # [Pa-s] reference viscosity for lower mantle
        pm.R = 8.314  # [J/mol-K] universal gas constant

        MantleLayer.__init__(self, pm.R_c0, pm.R_p0, params)

    def nu(self, T):
        '''
        Calculates viscosity in Pa-s of lower mantle given average bondary layer temperature

        :param T: average boundary layer temperature [K]
        :return: nu [Pa-s]
        '''

        pm = self.params.mantle
        return pm.nu_ref * np.exp(pm.Heff / pm.R * (1 / T - 1 / pm.T_ref))

    def Q_M(self, T_M):
        pass

    def energy_balance(self, T_M, T_B, T_U):
        pass

class Custom(Stevenson):
    def __init__(self, params=None):
        Stevenson.__init__(self, params=params)
        pm = self.params.mantle
        pm.Hp = 10e12 # [W] current heat production in mantle

    def heat_production(self, time):
        '''
        Overloaded method to use four-component radiogenics for mantle

        :param time: time [s]
        :return: heat production [W/m^3]
        '''
        pm = self.params.mantle
        return self.planet.radiogenics.heat_production_mantle(pm.Hp, time)/self.volume


