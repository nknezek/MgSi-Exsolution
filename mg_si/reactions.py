"""
Created on Sun Jul 17 17:11:58 2016

@author: nknezek
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import scipy.special as sp
from .base import Parameters

class Molar_Calculations():
    def __init__(self, species=None):
        molmass = {
            'Mg': 24.305,
            'Si': 28.0855,
            'Fe': 55.845,
            'O': 15.9994,
            'S': 32.065,
            'Ni': 58.6934,
            'H': 1.00794,
            'C': 12.0107, }
        molmass['MgO'] = molmass['Mg'] + molmass['O']
        molmass['FeO'] = molmass['Fe'] + molmass['O']
        molmass['SiO2'] = molmass['Si'] + 2 * molmass['O']
        molmass['MgSiO3'] = molmass['Mg'] + molmass['Si'] + 3 * molmass['O']
        molmass['FeSiO3'] = molmass['Fe'] + molmass['Si'] + 3 * molmass['O']
        self._molmass = molmass
        self.species = species
        if self.species is not None:
            self.molmass = [molmass[sp] for sp in self.species]
        else:
            self.molmass = None

    def _set_species(self, species):
        self.species = species
        self.molmass = [self._molmass[sp] for sp in species]

    def X2M(self, X, M_tot=None, wt_tot=None):
        '''compute moles given an array of mole fractions

        :param X:
        :return:
        '''
        assert (np.isclose(np.sum(X), 1.))
        if M_tot is None and wt_tot is None:
            raise ValueError('must provide M_tot or wt_tot')
        elif M_tot is not None and wt_tot is not None:
            raise ValueError('provide only one of M_tot or wt_tot')
        elif M_tot is None:
            M_tot = self.Xwtot2Mtot(X, wt_tot)
        return X*M_tot

    def X2wtp(self, X):
        '''computes wt% given array of mole fractions

        :param X:
        :param wt_tot:
        :return:
        '''
        assert (np.isclose(np.sum(X), 1.))
        tmp = X*self.molmass
        return tmp/np.sum(tmp)

    def X2wt(self, X, M_tot=None, wt_tot=None):
        '''computes total wt of each species given array of mole fractions

        :param X:
        :param M_tot:
        :return:
        '''
        assert (np.isclose(np.sum(X), 1.))
        if M_tot is None and wt_tot is None:
            raise ValueError('must provide M_tot or wt_tot')
        elif M_tot is not None and wt_tot is not None:
            raise ValueError('provide only one of M_tot or wt_tot')
        elif M_tot is None:
            M_tot = self.Xwtot2Mtot(X, wt_tot)
        return X*M_tot*self.molmass

    def XMtot2wtot(self, X, M_tot):
        ''' compute total mass given mol fractions and total moles

        :param X:
        :param M_tot:
        :return:
        '''
        assert(np.isclose(np.sum(X), 1.))
        return np.sum(X*self.molmass)*M_tot

    def Xwtot2Mtot(self, X, wt_tot):
        ''' computes total moles given mol fractions and total mass

        :param X:
        :param wt_tot:
        :return:
        '''
        assert (np.isclose(np.sum(X), 1.))
        return wt_tot/np.sum(X*self.molmass)

    def M2X(self, M):
        ''' computes mole fractions given array of moles

        :param M:
        :return:
        '''
        return M/np.sum(M)

    def M2wtp(self, M):
        ''' computes wtp given array of Moles

        :param M:
        :return:
        '''
        wt = M*self.molmass
        return wt/np.sum(wt)

    def dM2dwtp(self, dM, M):
        ''' computes wtp given array of Moles

        :param M:
        :return:
        '''
        dwt = dM*self.molmass
        wt = M*self.molmass
        return dwt/np.sum(wt)

    def M2wt(self, M):
        '''computes total wt of species given array of M

        :param M:
        :return:
        '''
        return M*self.molmass

    def wt2M(self, wt):
        '''compute Moles from wt

        :param wt:
        :return:
        '''
        return wt/self.molmass

    def wt2X(self, wt):
        '''compute mole fractions from wts

        :param wt:
        :return:
        '''
        M = wt/self.molmass
        return M/np.sum(M)

    def wt2wtp(self,wt):
        '''compute wtp from wt

        :param wt:
        :return:
        '''
        return wt/np.sum(wt)

    def wtp2M(self, wtp, M_tot=None, wt_tot=None):
        '''convert an array of wt% values to absolute Moles'''
        assert (np.isclose(np.sum(wtp), 1.))
        if M_tot is None and wt_tot is None:
            raise ValueError('must provide M_tot or wt_tot')
        elif M_tot is not None and wt_tot is not None:
            raise ValueError('provide only one of M_tot or wt_tot')
        elif wt_tot is None:
            wt_tot = self.wtpMtot2wtot(wtp, M_tot)
        return wtp*wt_tot/self.molmass

    def wtp2wt(self, wtp, M_tot=None, wt_tot=None):
        '''convert an array of wt% values to absolute wt

        :param wt:
        :return:
        '''
        assert (np.isclose(np.sum(wtp), 1.))
        if M_tot is None and wt_tot is None:
            raise ValueError('must provide M_tot or wt_tot')
        elif M_tot is not None and wt_tot is not None:
            raise ValueError('provide only one of M_tot or wt_tot')
        elif wt_tot is None:
            wt_tot = self.wtpMtot2wtot(wtp, M_tot)
        return wtp*wt_tot

    def wtp2X(self, wtp):
        '''compute mole fractions given array of wt%

        :param wtp:
        :param wt_tot:
        :return:
        '''
        assert (np.isclose(np.sum(wtp), 1.))
        tmp = wtp/self.molmass
        return tmp/np.sum(tmp)

    def wtpwtot2Mtot(self, wtp, wt_tot):
        '''compute total moles given wt % and total mass'''
        assert (np.isclose(np.sum(wtp), 1.))
        return np.sum(wtp / self.molmass) * wt_tot

    def wtpMtot2wtot(self, wtp, M_tot):
        '''compute total mass given wt % and total Moles

        :param wtp:
        :param M_tot:
        :return:
        '''
        return M_tot * np.sum(wtp) / np.sum(wtp / self.molmass)

class Core_MgSi(Molar_Calculations):
    '''class to compute molar concentrations in core for Mg + Si exsolution

    extends Molar_Calculations
    '''
    def __init__(self, species=None):
        if species is None:
            species = ['Mg','Si','O','Fe']
        super(Core_MgSi, self).__init__(species=species)

class Mantle_MgSi(Molar_Calculations):
    '''class to compute molar concentrations in mantle layer for Mg + Si exsolution

    extends Molar_Calculations
    '''
    def __init__(self, species=None):
        if species is None:
            species = ['MgO','SiO2','FeO','MgSiO3', 'FeSiO3']
        super(Mantle_MgSi, self).__init__(species=species)

    def compute_Mm_b(self, fraction_MgFe, X_MgFeO, X_SiO2, M_tot=None):
        ''' computes background mantle composition given Mg/(Mg+Fe) fraction, mol frac. ferropericlase, and mol frac. SiO2

        :param fraction_MgFe:
        :param X_MgFeO:
        :param X_SiO2:
        :param M_m:
        :return:
        '''
        X_MgO = X_MgFeO*fraction_MgFe
        X_FeO = X_MgFeO*(1-fraction_MgFe)
        X_MgFeSiO3 = (1-X_SiO2-X_MgFeO)
        X_MgSiO3 = X_MgFeSiO3*fraction_MgFe
        X_FeSiO3 = X_MgFeSiO3*(1-fraction_MgFe)
        Xm = np.array([X_MgO, X_SiO2, X_FeO, X_MgSiO3, X_FeSiO3])
        if M_tot is None:
            pr = self.params.reactions
            V_l = pr.V_l
            density = pr.density
            mass_l = density * V_l  # [kg]
        return self.mantle.X2M(Xm, wt_tot=mass_l)

class Reactions_MgSi():
    def __init__(self, params = None):
        if params is None:
            params = Parameters('Mantle reaction layer parameters')
        self.params = params
        try:
            self.params.reactions
        except:
            self.params.reactions = Parameters('Mantle reaction layer parameters')
        pr = self.params.reactions

        Cyr2s = 365.25*24*3600

        pr.thickness = 300 # [m] thickness of layer
        pr.density = 5500 # [kg/m^3] average density of layer in lower mantle

        pr.time_overturn = 800e6*Cyr2s # [s] overturn time of layer
        pr.V_c = 4/3*np.pi*3480e3**3 # [m^3] volume of total core
        pr.V_l = 4/3*np.pi*(3480e3+pr.thickness)**3 - pr.V_c # [m^3] volume of layer
        pr.d = 1. # [-] exponent of layer overturn expression
        pr.tau = pr.time_overturn # [s] constant of layer overturn expression
        pr.P = 135e9 # [Pa] pressure at CMB

        pr.dKMgSiO3_KMgSiO3 = 0.0 # these should be zero as no d/dT
        pr.dKFeSiO3_KFeSiO3 = 0.0 # these should be zero as no d/dT

        self.core = Core_MgSi()
        self.mantle = Mantle_MgSi()

    def _set_layer_thickness(self, thickness):
        pr = self.params.reactions
        pr.thickness = thickness
        pr.V_l = 4/3*np.pi*(3480e3+pr.thickness)**3 - pr.V_c # [m^3] volume of layer

    def dKs_dT(self, T_cmb, Moles):
        dKMgO_dT_KMgO = self.dKMgO_dT_KMgO(T_cmb, Moles)
        dKSiO2_dT_KSiO2 = self.dKSiO2_dT_KSiO2(T_cmb, Moles)
        dKFeO_dT_KFeO = self.dKFeO_dT_KFeO(T_cmb, Moles)
        dKMgSiO3_dT_KMgSiO3 = self.dKMgSiO3_dT_KMgSiO3(T_cmb, Moles)
        dKFeSiO3_dT_KFeSiO3 = self.dKFeSiO3_dT_KFeSiO3(T_cmb, Moles)
        dKs_dT = [dKMgO_dT_KMgO, dKSiO2_dT_KSiO2, dKFeO_dT_KFeO, dKMgSiO3_dT_KMgSiO3, dKFeSiO3_dT_KFeSiO3]
        return dKs_dT

    def dKMgO_dT_KMgO(self, T_cmb, Moles):
        ''' helper function to call Tushar's KD_MgO function

        :param T_cmb:
        :param Moles:
        :return:
        '''
        KMgO, dKMgO_dT = self.func_KD_MgO_val(T_cmb)
        return dKMgO_dT / KMgO

    def dKSiO2_dT_KSiO2(self, T_cmb, Moles):
        '''helper function to call Tushar's func_KD_SiO2

        :param T_cmb:
        :param Moles:
        :return:
        '''
        M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = self.unwrap_Moles(Moles)
        X_Si = M_Si/M_c
        X_O = M_O/M_c
        KSiO2, dKSiO2_dT = self.func_KD_SiO2_val(X_Si, X_O, T_cmb)
        return dKSiO2_dT / KSiO2

    def dKFeO_dT_KFeO(self, T_cmb, Moles):
        ''' helper function to call Tushar's func_KD_FeO'''
        KFeO, dKFeO_dT = self.func_KD_FeO_val(T_cmb)
        return dKFeO_dT / KFeO

    def dKMgSiO3_dT_KMgSiO3(self, T_cmb, Moles):
        '''computes dK_MgSiO3, normally = 0

        :param T_cmb:
        :param Moles:
        :return:
        '''
        pr = self.params.reactions
        dKMgSiO3_KMgSiO3 = pr.dKMgSiO3_KMgSiO3
        return dKMgSiO3_KMgSiO3

    def dKFeSiO3_dT_KFeSiO3(self, T_cmb, Moles):
        '''computes dK_FeSiO3, normally = 0

        :param T_cmb:
        :param Moles:
        :return:
        '''
        pr = self.params.reactions
        dKFeSiO3_KFeSiO3 = pr.dKFeSiO3_KFeSiO3
        return dKFeSiO3_KFeSiO3

    def func_KD_SiO2_val(self, X_Si, X_O, T_inp_base, P_inp_base=139e6, temp_diff_pm=1):
        '''
        K_D for SiO2 value
        Needs the following inputs : T_inp (Temperature in K), P_inp (Pressure in GPa),
        X_Si - mole frac Si in the core
        X_O - mole frac O in the core
        diff_pm = value of temp diff to calculate local deriv
        Note - (all the Mg related terms are zero since no data)
        '''
        P_inp = P_inp_base / 1e6  # convert to GPa
        ### Fit values from Hirose et al. 2017 paper (Eqn 5 in the Supplementary material)
        fit_KD_FeO_a = 0.3009  # (+/- 0.1120)
        fit_KD_FeO_b = 0
        fit_KD_FeO_c = -36.8332  # (+/- 5.5957)
        # Fit values from Rebecca Fischer et al. 2015 (extended Data Table 1 - Hirose 2017)
        fit_KD_Si_a = 1.3  # (+/- 0.3)
        fit_KD_Si_b = -13500  # (+/- 900)
        fit_KD_Si_c = 0
        ### Acitivity coeff fit values - Hirose et al. 2017, extended Data Table 1
        epsf_OO = -9.16  # (+/- 4.27)
        epsf_OSi = 7.73  # (+/- 4.53)
        epsf_SiSi = 0

        def func_KD(T_inp, X_Si=X_Si, X_O=X_O):
            log_KD_Feo = fit_KD_FeO_a + fit_KD_FeO_b / T_inp + fit_KD_FeO_c * P_inp / T_inp
            log_KD_Si = fit_KD_Si_a + fit_KD_Si_b / T_inp + fit_KD_Si_c * P_inp / T_inp

            # Steelmaking Handbook method for correcting gamma used:
            ##  (activity)  log(gamma(T)) = Tr/T*log(gamma(Tr))
            T_ref = 1873  # Kelvin, Hirose 2017

            #### Activity coeff values for gamma_Si (based on Hirose et al. 2017 + Ma 2001 Eqn 23-24 in the paper )
            ## Since the cross term due to Si-Si does not contribute,
            # so the only term that contributes for gamma_Si is the Si-O term (all the Mg related terms are zero since no data)
            sum_v = epsf_OSi * (X_O * (1. + np.log(1. - X_O) / X_O - 1. / (1. - X_O)) - X_O ** 2. * X_Si * (
            1. / (1. - X_Si) + 1. / (1. - X_O) + X_Si / (2. * (1. - X_Si) ** 2.) - 1.)) / 2.303
            log_gam_Si = -(T_ref / T_inp) * (epsf_SiSi * np.log(1. - X_Si) / 2.303 + sum_v)
            del sum_v
            #### Activity coeff values for gamma_O (based on Hirose et al. 2017 + Ma 2001 Eqn 23-24 in the paper )
            # Note - all the Mg related terms are zero since no data
            sum_v = epsf_OSi * (X_Si * (1. + np.log(1. - X_Si) / X_Si - 1. / (1. - X_Si)) - X_Si ** 2. * X_O * (
            1. / (1. - X_O) + 1. / (1. - X_Si) + X_O / (2. * (1. - X_O) ** 2.) - 1.)) / 2.303
            log_gam_O = -(T_ref / T_inp) * (epsf_OO * np.log(1. - X_O) / 2.303 + sum_v)
            del sum_v
            KD_SiO2 = (10. ** log_KD_Si) * ((10. ** log_KD_Feo) ** 2.) / (10. ** (log_gam_Si)) / (10. ** (
            log_gam_O)) ** 2.
            return KD_SiO2

        KD_SiO2 = func_KD(T_inp_base)
        KD_SiO2_T_deriv = (func_KD(T_inp_base + temp_diff_pm) - func_KD(T_inp_base - temp_diff_pm)) / (
        2. * temp_diff_pm)
        return KD_SiO2, KD_SiO2_T_deriv

    def func_KD_FeO_val(self, T_inp, P_inp_base=139e6):
        '''
        K_D for FeO value
        Needs the following inputs : T_inp (Temperature in K), P_inp (Pressure in GPa),
        '''
        P_inp = P_inp_base / 1e6  # convert to GPa
        ### Fit values from Hirose et al. 2017 paper (Eqn 5 in the Supplementary material)
        fit_KD_FeO_a = 0.3009  # (+/- 0.1120)
        fit_KD_FeO_b = 0
        fit_KD_FeO_c = -36.8332  # (+/- 5.5957)
        log_KD_Feo = fit_KD_FeO_a + fit_KD_FeO_b / T_inp + fit_KD_FeO_c * P_inp / T_inp
        KD_FeO = 10. ** log_KD_Feo
        KD_FeO_Tderiv = (KD_FeO) * -1. * fit_KD_FeO_c * P_inp * np.log(10.) / T_inp ** 2.
        return KD_FeO, KD_FeO_Tderiv

    def func_KD_MgO_val(self, T_inp, P_inp_base=139e6):
        '''
        K_D for MgO value
        Needs the following inputs : T_inp (Temperature in K), P_inp (Pressure in GPa),
        '''
        P_inp = P_inp_base / 1e6  # convert to GPa
        ### Fit values from Badro et al. 2015 paper (Eqn 5 in the Supplementary material)
        fit_KD_MgO_a = 1.23  # (+/- 0.7)
        fit_KD_MgO_b = -18816  # (+/- 2600)
        fit_KD_MgO_c = 0
        log_KD_Mgo = fit_KD_MgO_a + fit_KD_MgO_b / T_inp + fit_KD_MgO_c * P_inp / T_inp
        KD_MgO = 10. ** log_KD_Mgo
        KD_MgO_Tderiv = (KD_MgO) * -1. * fit_KD_MgO_b * np.log(10.) / T_inp ** 2.
        return KD_MgO, KD_MgO_Tderiv

    def unwrap_dKs_dT(self, dKs_dT):
        ''' helper function to unwrap dK values from dKs_dT'''
        dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs_dT
        return dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3
    
    def unwrap_Moles(self, Moles, return_sum=True, split_coremantle=False):
        ''' unwrap and compute Moles in core and mantle 
        
        :param Moles: 
        :return: 
        '''
        M_Mg, M_Si, M_Fe, M_O, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3 = list(Moles)
        M_c = np.sum(Moles[:4])
        M_m = np.sum(Moles[4:])
        if return_sum:
            if split_coremantle:
                return Moles[:4], Moles[4:], M_c, M_m
            else:
                return M_Mg, M_Si, M_O, M_Fe, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3,
        else:
            if split_coremantle:
                return Moles[:4], Moles[4:]
            else:
                return list(Moles)

    def erode_term(self, M_i, M_i_0, d=None, tau=None):
        ''' Layer erosion term given current and initial number of moles of species i and initial total moles in the layer

        :param M_i:
        :param M_i_0:
        :param M_m_0:
        :return:
        '''
        pr =self.params.reactions
        if d is None:
            d = pr.d
        if tau is None:
            tau = pr.tau
        return np.sign(M_i_0 - M_i) * M_i_0 / tau * ((np.abs(M_i - M_i_0) / M_i_0 + 1)**d-1)

    def dMi_b(self, Moles=None, dTdt=None):
        '''compute the erosion term incorporated directly into the equations

        :param Moles:
        :param dTdt:
        :return:
        '''
        pr = self.params.reactions
        M_Mg, M_Si, M_O, M_Fe, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = self.unwrap_Moles(Moles)
        M_Mg_b, M_Si_b, M_Fe_b, M_O_b, M_c_b, M_MgO_b, M_SiO2_b, M_FeO_b, M_MgSiO3_b, M_FeSiO3_b, M_m_b = self.unwrap_Moles(pr.Moles_b)
        dM_MgO_dt_b = -self.erode_term(M_MgO, M_MgO_b)/dTdt
        dM_SiO2_dt_b = -self.erode_term(M_SiO2, M_SiO2_b)/dTdt
        dM_FeO_dt_b = -self.erode_term(M_FeO, M_FeO_b)/dTdt
        dM_MgSiO3_dt_b = -self.erode_term(M_MgSiO3, M_MgSiO3_b)/dTdt
        dM_FeSiO3_dt_b = -self.erode_term(M_FeSiO3, M_FeSiO3_b)/dTdt

        # mantle visibility correction
        tau_m = pr.tau/100
        dM_MgO_dt_b += -self.erode_term(M_m, M_m_b, tau=tau_m)/dTdt*M_MgO/M_m
        dM_SiO2_dt_b += -self.erode_term(M_m, M_m_b, tau=tau_m)/dTdt*M_SiO2/M_m
        dM_FeO_dt_b += -self.erode_term(M_m, M_m_b, tau=tau_m)/dTdt*M_FeO/M_m
        dM_MgSiO3_dt_b += -self.erode_term(M_m, M_m_b, tau=tau_m)/dTdt*M_MgSiO3/M_m
        dM_FeSiO3_dt_b += -self.erode_term(M_m, M_m_b, tau=tau_m)/dTdt*M_FeSiO3/M_m
        return [0,0,0,0,0,dM_MgO_dt_b, dM_SiO2_dt_b, dM_FeO_dt_b, dM_MgSiO3_dt_b, dM_FeSiO3_dt_b,0]

    def dMi_dt_erode(self, Moles):
        '''calculate the erosion rate for each molar species in the mantle layer

        :param Moles:
        :return:
        '''
        pr = self.params.reactions

        # calculate the background mantle composition with proper Mg-Fe fraction
        M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = self.unwrap_Moles(Moles)
        fraction_MgFe = (M_MgO/(M_MgO+M_FeO) + M_MgSiO3/(M_MgSiO3+M_FeSiO3))/2
        M_Mg_0, M_Si_0, M_Fe_0, M_O_0, M_c_0, M_MgO_0, M_SiO2_0, M_FeO_0, M_MgSiO3_0, M_FeSiO3_0, M_m_0 = self.unwrap_Moles(pr.Moles_0)
        fraction_MgFe_0 = (M_MgO_0/(M_MgO_0+M_FeO_0) + M_MgSiO3_0/(M_MgSiO3_0+M_FeSiO3_0))/2
        fraction_MgFe_b = (fraction_MgFe_0+fraction_MgFe)/2
        X_MgFeO_0 = (M_MgO_0+M_FeO_0)/M_m_0
        Moles_b = self.compute_Moles_background(fraction_MgFe_b, X_MgFeO_0, M_SiO2_0/M_m_0, M_m_0)
        M_Mg_b, M_Si_b, M_Fe_b, M_O_b, M_c_b, M_MgO_b, M_SiO2_b, M_FeO_b, M_MgSiO3_b, M_FeSiO3_b, M_m_b = self.unwrap_Moles(Moles_b)

        # compute the erosional terms
        dM_MgO_dt_e = self.erode_term(M_MgO, M_MgO_b, M_m_b)
        dM_SiO2_dt_e = self.erode_term(M_SiO2, M_SiO2_b, M_m_b)
        dM_FeO_dt_e = self.erode_term(M_FeO, M_FeO_b, M_m_b)
        dM_MgSiO3_dt_e = self.erode_term(M_MgSiO3, M_MgSiO3_b, M_m_b)
        dM_FeSiO3_dt_e = self.erode_term(M_FeSiO3, M_FeSiO3_b, M_m_b)
        return [dM_MgO_dt_e, dM_SiO2_dt_e, dM_FeO_dt_e, dM_MgSiO3_dt_e, dM_FeSiO3_dt_e]

    def dMoles_dT(self, Moles=None, T_cmb=None, dKs_dT=None, dTdt=None, dMi_b=None):
        '''calcluate the change in Moles vs temperature T for each molar species in the core and mantle

        :param Moles:
        :param T_cmb:
        :param dKs_dT:
        :return:
        '''
        if dKs_dT is None:
            dKs_dT = self.compute_dKs_dT(Moles=Moles, T_cmb=T_cmb)
        if dMi_b is None:
            dMi_b = self.dMi_b(Moles=Moles, dTdt=dTdt)

        # core
        dM_Mg_dT = self.dM_Mg_dTc(Moles, dKs_dT, dMi_b)
        dM_Si_dT = self.dM_Si_dTc(Moles, dKs_dT, dMi_b)
        dM_O_dT = self.dM_O_dTc(Moles, dKs_dT, dMi_b)

        ## Allow Mg, Si, O to come out of core, but not go in.
        # M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = self.unwrap_Moles(Moles)
        # M_Mg_eq, M_Si_eq, M_O_eq = self.compute_Moles_eq(Moles=Moles, T_cmb=T_cmb)
        # dM_Mg_dT = self.logit(M_Mg-M_Mg_eq)*dM_Mg_dT
        # dM_Si_dT = self.logit(M_Si-M_Si_eq)*dM_Si_dT
        # dM_O_dT = self.logit(M_O-M_O_eq)*dM_O_dT

        dM_Fe_dT = self.dM_Fe_dTc(Moles, dKs_dT, dMi_b)
        dM_c_dT = np.sum([dM_Mg_dT, dM_Si_dT, dM_Fe_dT, dM_O_dT])

        # mantle
        dM_MgO_dT = self.dM_MgO_dTc(Moles, dKs_dT, dMi_b)
        dM_SiO2_dT = self.dM_SiO2_dTc(Moles, dKs_dT, dMi_b)
        dM_FeO_dT = self.dM_FeO_dTc(Moles, dKs_dT, dMi_b)
        dM_MgSiO3_dT = self.dM_MgSiO3_dTc(Moles, dKs_dT, dMi_b)
        dM_FeSiO3_dT = self.dM_FeSiO3_dTc(Moles, dKs_dT, dMi_b)
        dM_m_dT = np.sum([dM_MgO_dT, dM_FeO_dT, dM_SiO2_dT, dM_MgSiO3_dT, dM_FeSiO3_dT])
        return [dM_Mg_dT, dM_Si_dT, dM_Fe_dT, dM_O_dT, dM_c_dT, dM_MgO_dT, dM_SiO2_dT, dM_FeO_dT, dM_MgSiO3_dT,
                dM_FeSiO3_dT, dM_m_dT]

    def compute_Moles_eq(self, Moles=None, T_cmb=None):
        M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = self.unwrap_Moles(Moles)
        X_Si = M_Si/M_c
        X_O = M_O/M_c
        K_Si,_ = self.func_KD_SiO2_val(X_Si, X_O, T_inp_base=T_cmb)
        K_Fe, _ = self.func_KD_FeO_val(T_cmb)
        K_Mg, _ = self.func_KD_MgO_val(T_cmb)
        M_Mg_eq = K_Mg*M_MgO*M_c**2 / (M_O*M_m)
        M_Si_eq = K_Si*M_SiO2*M_c**3 / (M_O**2 * M_m)
        M_O_eq = K_Fe*M_FeO*M_c**2 / (M_Fe*M_m)
        return M_Mg_eq, M_Si_eq, M_O_eq

    def dMoles_dt(self, Moles=None, T_cmb=None, dTc_dt=None, dKs_dT=None, dMoles_dT=None):
        '''calculate the change in Moles vs time (t) for each molar species in the core and mantle

        :param Moles:
        :param T_cmb:
        :param dTc_dt:
        :return:
        '''
        if dMoles_dT is None:
            dMoles_dT = self.dMoles_dT(Moles=Moles, T_cmb=T_cmb, dKs_dT=dKs_dT, dTdt=dTc_dt)

        dM_Mg_dT, dM_Si_dT, dM_Fe_dT, dM_O_dT, dM_MgO_dT, dM_SiO2_dT, \
                dM_FeO_dT, dM_MgSiO3_dT, dM_FeSiO3_dT = self.unwrap_Moles(dMoles_dT)
        # core
        dM_Mg_dt = dM_Mg_dT*dTc_dt
        dM_Si_dt = dM_Si_dT*dTc_dt
        dM_O_dt = dM_O_dT*dTc_dt
        dM_Fe_dt = dM_Fe_dT*dTc_dt

        # mantle
        dM_MgO_dt = dM_MgO_dT*dTc_dt
        dM_SiO2_dt = dM_SiO2_dT * dTc_dt
        dM_FeO_dt = dM_FeO_dT * dTc_dt
        dM_MgSiO3_dt = dM_MgSiO3_dT * dTc_dt
        dM_FeSiO3_dt = dM_FeSiO3_dT * dTc_dt
        return [dM_Mg_dt, dM_Si_dt, dM_O_dt, dM_Fe_dt, dM_MgO_dt, dM_SiO2_dt, dM_FeO_dt, dM_MgSiO3_dt, dM_FeSiO3_dt]

    def dM_MgSiO3_dTc(self, Moles, dKs, dMi_b):
        dM_Mg_er, dM_Si_er, dM_Fe_er, dM_O_er, dM_c_er, dM_MgO_er, dM_SiO2_er, dM_FeO_er, dM_MgSiO3_er, dM_FeSiO3_er, dM_m_er = self.unwrap_Moles(
            dMi_b)
        M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = Moles
        dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs
        return M_MgSiO3 * (M_Fe * (M_O * (M_FeO * (M_MgO * (M_SiO2 * (
        4.0 * dM_FeO_er + 8.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 8.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
                                                            -4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                   -4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er)) + M_FeSiO3 * (M_FeO * (
        M_SiO2 * (4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er) + M_m * (-4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er)) + M_MgO * (
                                                                                                         M_SiO2 * (
                                                                                                         4.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 8.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
                                                                                                         4.0 * dM_FeO_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                                         -4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er))) + dKFeO_KFeO * (
                                   M_Mg * M_O * (4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                   M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                                   -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_Si * (M_Mg * (
                                   M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                   M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                                   M_FeO * (2.0 * M_SiO2 + 10.0 * M_m) + M_FeSiO3 * (
                                   -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m))) + M_MgO * M_O * (M_FeO * (
                                   3.0 * M_SiO2 + 6.0 * M_m) + M_FeSiO3 * (-6.0 * M_SiO2 + 15.0 * M_m))) + M_c * (
                                   M_Mg * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                   M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) + M_O * (
                                   M_FeO + M_MgO - M_m) - M_SiO2 * M_m)) + M_MgO * M_O * (
                                   -M_FeO * M_SiO2 - M_FeSiO3 * M_m) + M_Si * (M_Mg * (
                                   M_FeO * (-M_SiO2 - 3.0 * M_m) + M_FeSiO3 * (
                                   4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_O * (
                                   M_FeO + M_FeSiO3)) + M_MgO * (M_FeO * (-2.0 * M_SiO2 - 2.0 * M_m) + M_FeSiO3 * (
                                   2.0 * M_SiO2 - 6.0 * M_m)))))) + M_FeSiO3 * dKFeSiO3_KFeSiO3 * (M_O * (
        M_Fe * M_FeO * M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_Mg * (M_Fe * (
        M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
        -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeO * M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m))) + M_Si * (
                                                                                                   M_Mg * (M_Fe * (
                                                                                                   M_FeO * (
                                                                                                   -M_SiO2 + M_m) + M_MgO * (
                                                                                                   -M_SiO2 + M_m) + M_O * (
                                                                                                   -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_SiO2 * M_m) + M_FeO * (
                                                                                                           M_MgO * (
                                                                                                           -M_SiO2 + M_m) + M_O * (
                                                                                                           -9.0 * M_MgO - 6.0 * M_SiO2 + 15.0 * M_m))) + M_MgO * (
                                                                                                   M_Fe * (M_FeO * (
                                                                                                   -M_SiO2 + M_m) + M_O * (
                                                                                                           -9.0 * M_FeO - 6.0 * M_SiO2 + 15.0 * M_m)) + M_FeO * M_O * (
                                                                                                   -9.0 * M_SiO2 + 9.0 * M_m))) + M_c * (
                                                                                                   M_Mg * (M_Fe * (
                                                                                                   M_FeO * (
                                                                                                   M_SiO2 - M_m) + M_MgO * (
                                                                                                   M_SiO2 - M_m) + M_O * (
                                                                                                   M_FeO + M_MgO - M_m) - M_SiO2 * M_m) + M_FeO * (
                                                                                                           M_MgO * (
                                                                                                           M_SiO2 - M_m) + M_O * (
                                                                                                           M_MgO - M_m))) + M_MgO * (
                                                                                                   M_Fe * (M_FeO * (
                                                                                                   M_SiO2 - M_m) + M_O * (
                                                                                                           M_FeO - M_m)) + M_FeO * M_O * (
                                                                                                   M_SiO2 - M_m)) + M_Si * (
                                                                                                   M_Mg * (M_Fe * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO + M_O + M_SiO2 - 9.0 * M_m) + M_FeO * (
                                                                                                           4.0 * M_MgO + 2.0 * M_SiO2 - 6.0 * M_m)) + M_MgO * (
                                                                                                   M_Fe * (
                                                                                                   4.0 * M_FeO + 2.0 * M_SiO2 - 6.0 * M_m) + M_FeO * (
                                                                                                   4.0 * M_SiO2 - 4.0 * M_m))))) + M_Mg * (
                           M_O * (M_Fe * (M_FeO * (M_SiO2 * (
                           4.0 * dM_FeO_er + 8.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 8.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
                                                   -4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)) + M_MgO * (
                                          M_SiO2 * (
                                          4.0 * dM_FeO_er + 8.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 8.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
                                          -4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                          -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er)) + M_FeO * (
                                  M_MgO * (M_SiO2 * (
                                  4.0 * dM_FeO_er + 8.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 8.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
                                           -4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                  -4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er)) + M_FeSiO3 * (M_FeO * (
                           M_SiO2 * (4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er) + M_m * (
                           -4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er)) + M_MgO * (M_SiO2 * (
                           4.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 8.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
                                                                              4.0 * dM_FeO_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                        -4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er))) + dKMgO_KMgO * (
                           M_Fe * M_O * (-4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                           M_FeO * (4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_Si * (M_Fe * (
                           -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                           M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
                           M_FeO * (M_SiO2 - 4.0 * M_m) + M_FeSiO3 * (
                           9.0 * M_FeO + 9.0 * M_MgO + 4.0 * M_SiO2 - 25.0 * M_m) + M_MgO * (
                           3.0 * M_SiO2 + 6.0 * M_m) - 9.0 * M_SiO2 * M_m)) + M_O * (M_FeO * (
                           M_MgO * (3.0 * M_SiO2 + 6.0 * M_m) - 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (M_FeO * (
                           9.0 * M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                                  3.0 * M_SiO2 + 6.0 * M_m) - 9.0 * M_SiO2 * M_m))) + M_c * (
                           M_Fe * (M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                           M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                                   M_FeSiO3 * (-M_FeO - M_MgO + M_m) + M_SiO2 * (-M_FeO - M_MgO + M_m))) + M_O * (
                           M_FeO * M_SiO2 * (-M_MgO + M_m) + M_FeSiO3 * (
                           M_FeO * (-M_SiO2 + M_m) + M_SiO2 * (-M_MgO + M_m))) + M_Si * (M_Fe * (
                           M_FeO * (-M_SiO2 + M_m) + M_FeSiO3 * (
                           -4.0 * M_FeO - 4.0 * M_MgO - M_SiO2 + 9.0 * M_m) + M_MgO * (
                           -2.0 * M_SiO2 - 2.0 * M_m) + M_O * (
                           -M_FeO - M_FeSiO3 - M_SiO2 + M_m) + 4.0 * M_SiO2 * M_m) + M_FeO * (M_MgO * (
                           -2.0 * M_SiO2 - 2.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (M_FeO * (
                           -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                                                                                          -2.0 * M_SiO2 - 2.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_O * (
                                                                                         M_FeO * (
                                                                                         -M_SiO2 + M_m) + M_FeSiO3 * (
                                                                                         -M_SiO2 + M_m)))))) + M_Si * (
                           M_Fe * (M_FeO * (M_MgO * (
                           M_SiO2 * (dM_FeO_er + 2 * dM_FeSiO3_er + dM_MgO_er + 2 * dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                           -1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                            -dM_MgO_er - dM_MgSiO3_er)) + M_FeSiO3 * (M_FeO * (
                           M_SiO2 * (dM_MgO_er + dM_MgSiO3_er) + M_m * (-dM_MgO_er - dM_MgSiO3_er)) + M_MgO * (
                                                                                      M_SiO2 * (
                                                                                      1.0 * dM_FeSiO3_er + dM_MgO_er + 2.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er) + M_m * (
                                                                                      1.0 * dM_FeO_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                      -dM_MgO_er - dM_MgSiO3_er)) + M_O * (
                                   M_FeO * (M_MgO * (
                                   6.0 * dM_FeO_er + 15.0 * dM_FeSiO3_er + 6.0 * dM_MgO_er + 15.0 * dM_MgSiO3_er + 9.0 * dM_SiO2_er) + M_SiO2 * (
                                            dM_MgO_er + dM_MgSiO3_er) + M_m * (
                                            -4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er)) + M_FeSiO3 * (
                                   M_FeO * (9.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er) + M_MgO * (
                                   6.0 * dM_FeO_er + 15.0 * dM_FeSiO3_er + 15.0 * dM_MgO_er + 24.0 * dM_MgSiO3_er + 9.0 * dM_SiO2_er) + M_SiO2 * (
                                   4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er) + M_m * (
                                   -25.0 * dM_MgO_er - 25.0 * dM_MgSiO3_er)) + M_MgO * (M_SiO2 * (
                                   6.0 * dM_FeO_er + 15.0 * dM_FeSiO3_er + 9.0 * dM_MgO_er + 18.0 * dM_MgSiO3_er + 9.0 * dM_SiO2_er) + M_m * (
                                                                                        -6.0 * dM_FeO_er - 15.0 * dM_FeSiO3_er - 9.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                   -9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er))) + M_Mg * (M_Fe * (M_FeO * (
                           M_SiO2 * (dM_FeO_er + 2 * dM_FeSiO3_er + dM_MgO_er + 2 * dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                           -1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er)) + M_MgO * (M_SiO2 * (
                           dM_FeO_er + 2 * dM_FeSiO3_er + dM_MgO_er + 2 * dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                                                                                                    -1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er)) + M_O * (
                                                                                              M_FeO * (
                                                                                              6.0 * dM_FeO_er + 15.0 * dM_FeSiO3_er + 6.0 * dM_MgO_er + 15.0 * dM_MgSiO3_er + 9.0 * dM_SiO2_er) + M_MgO * (
                                                                                              6.0 * dM_FeO_er + 15.0 * dM_FeSiO3_er + 6.0 * dM_MgO_er + 15.0 * dM_MgSiO3_er + 9.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                              4.0 * dM_FeO_er + 10.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 10.0 * dM_MgSiO3_er + 6.0 * dM_SiO2_er) + M_m * (
                                                                                              -10.0 * dM_FeO_er - 25.0 * dM_FeSiO3_er - 10.0 * dM_MgO_er - 25.0 * dM_MgSiO3_er - 15.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                              -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er)) + M_FeO * (
                                                                                      M_MgO * (M_SiO2 * (
                                                                                      dM_FeO_er + 2 * dM_FeSiO3_er + dM_MgO_er + 2 * dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                                                                                               -1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                      -1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er)) + M_FeSiO3 * (
                                                                                      M_FeO * (M_SiO2 * (
                                                                                      1.0 * dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_m * (
                                                                                               -1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er)) + M_MgO * (
                                                                                      M_SiO2 * (
                                                                                      dM_FeSiO3_er + dM_MgO_er + 2 * dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                                                                                      1.0 * dM_FeO_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                      -1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er)) + M_O * (
                                                                                      M_FeO * (M_MgO * (
                                                                                      6.0 * dM_FeO_er + 15.0 * dM_FeSiO3_er + 6.0 * dM_MgO_er + 15.0 * dM_MgSiO3_er + 9.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                               6.0 * dM_FeO_er + 12.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 10.0 * dM_MgSiO3_er + 6.0 * dM_SiO2_er) + M_m * (
                                                                                               -15.0 * dM_FeSiO3_er - 10.0 * dM_MgO_er - 25.0 * dM_MgSiO3_er - 15.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                                      M_FeO * (
                                                                                      -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er) + M_MgO * (
                                                                                      -3.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er + 6.0 * dM_MgO_er + 15.0 * dM_MgSiO3_er + 9.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                      6.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 10.0 * dM_MgSiO3_er + 6.0 * dM_SiO2_er) + M_m * (
                                                                                      15.0 * dM_FeO_er - 10.0 * dM_MgO_er - 25.0 * dM_MgSiO3_er - 15.0 * dM_SiO2_er)))) + M_O * (
                           M_FeO * (M_MgO * (M_SiO2 * (
                           9.0 * dM_FeO_er + 18.0 * dM_FeSiO3_er + 9.0 * dM_MgO_er + 18.0 * dM_MgSiO3_er + 9.0 * dM_SiO2_er) + M_m * (
                                             -9.0 * dM_FeSiO3_er - 9.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                    -9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er)) + M_FeSiO3 * (M_FeO * (
                           M_SiO2 * (9.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er) + M_m * (
                           -9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er)) + M_MgO * (M_SiO2 * (
                           9.0 * dM_FeSiO3_er + 9.0 * dM_MgO_er + 18.0 * dM_MgSiO3_er + 9.0 * dM_SiO2_er) + M_m * (
                                                                              9.0 * dM_FeO_er - 9.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                          -9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er))) + dKSiO2_KSiO2 * (
                           M_O * (M_Fe * M_MgO * (
                           M_FeO * (-2.0 * M_SiO2 - 4.0 * M_m) + M_FeSiO3 * (4.0 * M_SiO2 - 10.0 * M_m)) + M_Mg * (
                                  M_Fe * (M_FeO * (-2.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                  -2.0 * M_SiO2 - 4.0 * M_m) + 6.0 * M_SiO2 * M_m) + M_FeO * (
                                  M_MgO * (-2.0 * M_SiO2 - 4.0 * M_m) + 6.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                  M_FeO * (-6.0 * M_SiO2 + 6.0 * M_m) + M_MgO * (
                                  -2.0 * M_SiO2 - 4.0 * M_m) + 6.0 * M_SiO2 * M_m))) + M_c * (M_Mg * (M_Fe * (
                           M_FeO * (M_SiO2 + M_m) + M_MgO * (M_SiO2 + M_m) + M_O * (
                           -M_FeO - M_MgO + M_m) - 2.0 * M_SiO2 * M_m) + M_FeO * (M_MgO * (
                           M_SiO2 + M_m) - 2.0 * M_SiO2 * M_m) + M_FeSiO3 * (M_FeO * (
                           2.0 * M_SiO2 - 2.0 * M_m) + M_MgO * (M_SiO2 + M_m) - 2.0 * M_SiO2 * M_m) + M_O * (M_FeO * (
                           -M_MgO + M_m) + M_FeSiO3 * (-M_MgO + M_m))) + M_MgO * (M_Fe * (
                           M_FeO * (M_SiO2 + M_m) + M_FeSiO3 * (-M_SiO2 + 3.0 * M_m) + M_O * (
                           -M_FeO - M_FeSiO3 - M_SiO2 + M_m)) + M_O * (M_FeO * (-M_SiO2 + M_m) + M_FeSiO3 * (
                           -M_SiO2 + M_m)))))) + M_c * (M_Fe * (M_FeO * (M_MgO * (
        M_SiO2 * (-dM_FeO_er - 2.0 * dM_FeSiO3_er - dM_MgO_er - 2.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
        dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_SiO2 * M_m * (dM_MgO_er + dM_MgSiO3_er)) + M_FeSiO3 * (M_FeO * (
        M_SiO2 * (-dM_MgO_er - dM_MgSiO3_er) + M_m * (dM_MgO_er + dM_MgSiO3_er)) + M_MgO * (M_SiO2 * (
        -dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                                            -dM_FeO_er + 1.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                                              dM_MgO_er + dM_MgSiO3_er)) + M_O * (
                                                                M_FeO * (M_MgO * (
                                                                -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_SiO2 * (
                                                                         -dM_MgO_er - dM_MgSiO3_er)) + M_FeSiO3 * (
                                                                M_FeO * (-dM_MgO_er - dM_MgSiO3_er) + M_MgO * (
                                                                -dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                dM_MgO_er + dM_MgSiO3_er)) + M_MgO * (M_SiO2 * (
                                                                -dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                                                      dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                dM_MgO_er + dM_MgSiO3_er))) + M_Mg * (M_Fe * (M_FeO * (
        M_SiO2 * (-dM_FeO_er - 2.0 * dM_FeSiO3_er - dM_MgO_er - 2.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
        dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_MgO * (M_SiO2 * (
        -dM_FeO_er - 2.0 * dM_FeSiO3_er - dM_MgO_er - 2.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
                                                              dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_O * (
                                                                                                              M_FeO * (
                                                                                                              -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_MgO * (
                                                                                                              -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                                                              dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                                              dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er)) + M_FeO * (
                                                                                                      M_MgO * (
                                                                                                      M_SiO2 * (
                                                                                                      -dM_FeO_er - 2.0 * dM_FeSiO3_er - dM_MgO_er - 2.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
                                                                                                      dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                                      dM_MgO_er + 1.0 * dM_MgSiO3_er)) + M_FeSiO3 * (
                                                                                                      M_FeO * (
                                                                                                      M_SiO2 * (
                                                                                                      -dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_m * (
                                                                                                      dM_MgO_er + 1.0 * dM_MgSiO3_er)) + M_MgO * (
                                                                                                      M_SiO2 * (
                                                                                                      -1.0 * dM_FeSiO3_er - dM_MgO_er - 2.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
                                                                                                      -dM_FeO_er + dM_MgSiO3_er + dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                                      dM_MgO_er + 1.0 * dM_MgSiO3_er)) + M_O * (
                                                                                                      M_FeO * (M_MgO * (
                                                                                                      -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                                                               dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_FeSiO3 * (
                                                                                                      M_FeO * (
                                                                                                      dM_FeO_er + dM_FeSiO3_er) + M_MgO * (
                                                                                                      dM_FeO_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                                                      -dM_FeO_er + dM_MgSiO3_er + dM_SiO2_er)))) + M_O * (
                                                        M_FeO * (M_MgO * (M_SiO2 * (
                                                        -dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                          dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                 dM_MgO_er + dM_MgSiO3_er)) + M_FeSiO3 * (M_FeO * (
                                                        M_SiO2 * (-dM_MgO_er - dM_MgSiO3_er) + M_m * (
                                                        dM_MgO_er + dM_MgSiO3_er)) + M_MgO * (M_SiO2 * (
                                                        -dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                                              -dM_FeO_er + dM_MgSiO3_er + dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                                          dM_MgO_er + dM_MgSiO3_er))) + M_Si * (
                                                        M_Fe * (M_FeO * (M_MgO * (
                                                        -2.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_SiO2 * (
                                                                         -dM_MgO_er - dM_MgSiO3_er) + M_m * (
                                                                         dM_MgO_er + dM_MgSiO3_er)) + M_FeSiO3 * (
                                                                M_FeO * (
                                                                -4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_MgO * (
                                                                -2.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 10.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_SiO2 * (
                                                                -dM_MgO_er - dM_MgSiO3_er) + M_m * (
                                                                9.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er)) + M_MgO * (
                                                                M_SiO2 * (
                                                                -2.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 8.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_m * (
                                                                2.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er)) + M_O * (
                                                                M_FeO * (-dM_MgO_er - dM_MgSiO3_er) + M_FeSiO3 * (
                                                                -dM_MgO_er - dM_MgSiO3_er) + M_SiO2 * (
                                                                -dM_MgO_er - dM_MgSiO3_er) + M_m * (
                                                                dM_MgO_er + dM_MgSiO3_er)) + M_SiO2 * M_m * (
                                                                4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er)) + M_FeO * (
                                                        M_MgO * (M_SiO2 * (
                                                        -4.0 * dM_FeO_er - 8.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 8.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_m * (
                                                                 4.0 * dM_FeSiO3_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                        4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er)) + M_FeSiO3 * (M_FeO * (
                                                        M_SiO2 * (-4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_m * (
                                                        4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er)) + M_MgO * (M_SiO2 * (
                                                        -4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 8.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_m * (
                                                                                                          -4.0 * dM_FeO_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                                             4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er)) + M_Mg * (
                                                        M_Fe * (M_FeO * (
                                                        -2.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_MgO * (
                                                                -2.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_SiO2 * (
                                                                -dM_FeO_er - 3.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                                                                3.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er + 3.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er + 6.0 * dM_SiO2_er)) + M_FeO * (
                                                        M_MgO * (
                                                        -2.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_SiO2 * (
                                                        -2 * dM_FeO_er - 4.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                                                        6.0 * dM_FeSiO3_er + 3.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er + 6.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                        M_FeO * (4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er) + M_MgO * (
                                                        2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_SiO2 * (
                                                        -2.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                                                        -6.0 * dM_FeO_er + 3.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er + 6.0 * dM_SiO2_er)) + M_O * (
                                                        M_FeO * (dM_FeO_er + dM_FeSiO3_er) + M_FeSiO3 * (
                                                        dM_FeO_er + dM_FeSiO3_er))) + M_O * (M_FeO * (
                                                        M_SiO2 * (-dM_MgO_er - dM_MgSiO3_er) + M_m * (
                                                        dM_MgO_er + dM_MgSiO3_er)) + M_FeSiO3 * (M_SiO2 * (
                                                        -dM_MgO_er - dM_MgSiO3_er) + M_m * (
                                                                                                 dM_MgO_er + dM_MgSiO3_er))))) + dKMgSiO3_KMgSiO3 * (
                           M_O * (M_Fe * M_MgO * (-4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                           M_FeO * (4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_Mg * (M_Fe * (M_FeSiO3 * (
                           M_FeO * (4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                                                                                                       -4.0 * M_FeO - 4.0 * M_MgO)) + M_MgO * (
                                                                                               -4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                               M_FeO * (
                                                                                               4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)))) + M_Si * (
                           M_Mg * (M_Fe * (
                           M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
                           M_FeO * (M_SiO2 - 4.0 * M_m) + M_FeSiO3 * (
                           9.0 * M_FeO + 9.0 * M_MgO + 4.0 * M_SiO2 - 25.0 * M_m) + M_MgO * (
                           M_SiO2 - 4.0 * M_m) - 9.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (-M_FeO - M_MgO)) + M_MgO * (
                                   -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
                                   M_FeO * (M_MgO * (M_SiO2 - 4.0 * M_m) - 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                   M_FeO * (9.0 * M_MgO + 9.0 * M_SiO2 - 9.0 * M_m) + M_MgO * (
                                   M_SiO2 - 4.0 * M_m) - 9.0 * M_SiO2 * M_m))) + M_MgO * (M_Fe * (
                           -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
                           M_FeO * (M_SiO2 - 4.0 * M_m) + M_FeSiO3 * (
                           9.0 * M_FeO + 4.0 * M_SiO2 - 25.0 * M_m) - 9.0 * M_SiO2 * M_m)) + M_O * (
                                                                                          -9.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                          M_FeO * (
                                                                                          9.0 * M_SiO2 - 9.0 * M_m) - 9.0 * M_SiO2 * M_m)))) + M_c * (
                           M_Mg * (M_Fe * (
                           M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                           M_FeSiO3 * (-M_FeO - M_MgO + M_m) + M_SiO2 * (-M_FeO - M_MgO + M_m)) + M_SiO2 * M_m * (
                           M_FeO + M_MgO)) + M_MgO * (
                                   M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (
                                   M_FeO * M_SiO2 * (-M_MgO + M_m) + M_FeSiO3 * (
                                   M_FeO * (-M_MgO - M_SiO2 + M_m) + M_SiO2 * (-M_MgO + M_m)))) + M_MgO * (M_Fe * (
                           M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                           M_FeSiO3 * (-M_FeO + M_m) + M_SiO2 * (-M_FeO + M_m))) + M_O * (
                                                                                                           M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           -M_SiO2 + M_m) + M_SiO2 * M_m))) + M_Si * (
                           M_Mg * (M_Fe * (M_FeO * (-M_SiO2 + M_m) + M_FeSiO3 * (
                           -4.0 * M_FeO - 4.0 * M_MgO - M_SiO2 + 9.0 * M_m) + M_MgO * (-M_SiO2 + M_m) + M_O * (
                                           -M_FeO - M_FeSiO3 - M_MgO - M_SiO2 + M_m) + 4.0 * M_SiO2 * M_m) + M_FeO * (
                                   M_MgO * (-M_SiO2 + M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                   M_FeO * (-4.0 * M_MgO - 4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                                   -M_SiO2 + M_m) + 4.0 * M_SiO2 * M_m) + M_O * (
                                   M_FeO * (-M_MgO - M_SiO2 + M_m) + M_FeSiO3 * (-M_MgO - M_SiO2 + M_m))) + M_MgO * (
                           M_Fe * (M_FeO * (-M_SiO2 + M_m) + M_FeSiO3 * (-4.0 * M_FeO - M_SiO2 + 9.0 * M_m) + M_O * (
                           -M_FeO - M_FeSiO3 - M_SiO2 + M_m) + 4.0 * M_SiO2 * M_m) + 4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                           M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_O * (
                           M_FeO * (-M_SiO2 + M_m) + M_FeSiO3 * (-M_SiO2 + M_m))))))) / (M_O * (M_Fe * (M_MgO * (
        4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
        M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (M_FeO * (
        M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (M_FeO * (
        -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m))) + M_Mg * (M_Fe * (
        M_FeSiO3 * (
        M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
        M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
        -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (4.0 * M_FeO + 4.0 * M_MgO)) + M_MgO * (
                                                                                                           4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           M_MgO * (
                                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)))) + M_Si * (
                                                                                         M_Fe * (M_MgO * (
                                                                                         M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                         M_FeO * (
                                                                                         -M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                 M_FeO * (M_MgO * (
                                                                                                 -M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                 M_FeO * (
                                                                                                 -M_SiO2 + M_m) + M_MgO * (
                                                                                                 -M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (
                                                                                                 M_MgO * (M_FeO * (
                                                                                                 -M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                                                                                                          -9.0 * M_FeO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                 M_FeO * (
                                                                                                 -9.0 * M_MgO - M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                                                                                                 -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                                                                                                 -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m))) + M_Mg * (
                                                                                         M_Fe * (M_FeSiO3 * (M_FeO * (
                                                                                         -M_SiO2 + M_m) + M_MgO * (
                                                                                                             -M_SiO2 + M_m) + M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                 M_FeO * (
                                                                                                 -M_SiO2 + M_m) + M_MgO * (
                                                                                                 -M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                                                                                                 M_FeO * (
                                                                                                 -M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                                                                                                 -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                                                                                                 -M_SiO2 + 4.0 * M_m) + M_MgSiO3 * (
                                                                                                 -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                                                                                                 M_FeO + M_MgO)) + M_MgO * (
                                                                                         M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                         M_FeO * (
                                                                                         -M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                         M_FeO * (M_MgO * (
                                                                                         -M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                         M_FeO * (
                                                                                         -M_SiO2 + M_m) + M_MgO * (
                                                                                         -M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (
                                                                                         M_FeO * (M_MgO * (
                                                                                         -M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                         M_FeO * (
                                                                                         -9.0 * M_MgO - 9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (
                                                                                         -M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                         M_FeO * (
                                                                                         -9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_FeSiO3 * (
                                                                                         -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m)))) + M_O * (
                                                                                         M_MgO * (
                                                                                         9.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                         M_FeO * (
                                                                                         -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                         M_FeO * (M_MgO * (
                                                                                         -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                         M_FeO * (
                                                                                         -9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (
                                                                                         -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)))) + M_c * (
                                                                                         M_Fe * (M_MgO * (
                                                                                         -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                         M_FeO * (
                                                                                         M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                 M_FeO * (M_MgO * (
                                                                                                 M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                 M_FeO * (
                                                                                                 M_SiO2 - M_m) + M_MgO * (
                                                                                                 M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
                                                                                                 M_MgO * (M_FeSiO3 * (
                                                                                                 M_FeO - M_m) + M_SiO2 * (
                                                                                                          M_FeO - M_m)) + M_MgSiO3 * (
                                                                                                 M_FeO * (
                                                                                                 M_MgO + M_SiO2) + M_FeSiO3 * (
                                                                                                 M_FeO + M_MgO - M_m) + M_MgO * (
                                                                                                 M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Mg * (
                                                                                         M_Fe * (M_FeSiO3 * (M_FeO * (
                                                                                         M_SiO2 - M_m) + M_MgO * (
                                                                                                             M_SiO2 - M_m) - M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                 M_FeO * (
                                                                                                 M_SiO2 - M_m) + M_MgO * (
                                                                                                 M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
                                                                                                 M_FeSiO3 * (
                                                                                                 M_FeO + M_MgO - M_m) + M_MgSiO3 * (
                                                                                                 M_FeO + M_MgO - M_m) + M_SiO2 * (
                                                                                                 M_FeO + M_MgO - M_m)) + M_SiO2 * M_m * (
                                                                                                 -M_FeO - M_MgO)) + M_MgO * (
                                                                                         -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                         M_FeO * (
                                                                                         M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                         M_FeO * (M_MgO * (
                                                                                         M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                         M_FeO * (
                                                                                         M_SiO2 - M_m) + M_MgO * (
                                                                                         M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
                                                                                         M_FeO * M_SiO2 * (
                                                                                         M_MgO - M_m) + M_FeSiO3 * (
                                                                                         M_FeO * (
                                                                                         M_MgO + M_SiO2 - M_m) + M_SiO2 * (
                                                                                         M_MgO - M_m)) + M_MgSiO3 * (
                                                                                         M_FeO * (
                                                                                         M_MgO - M_m) + M_FeSiO3 * (
                                                                                         M_FeO + M_MgO - M_m)))) + M_O * (
                                                                                         M_MgO * (
                                                                                         -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                         M_FeO * (
                                                                                         M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                         M_FeO * (M_MgO * (
                                                                                         M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                         M_FeO * (
                                                                                         M_SiO2 - M_m) + M_MgO * (
                                                                                         M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Si * (
                                                                                         M_Fe * (M_MgO * (M_FeO * (
                                                                                         M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                          4.0 * M_FeO + M_SiO2 - 9.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                 M_FeO * (
                                                                                                 4.0 * M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                 4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                                 4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_O * (
                                                                                                 M_MgO * (
                                                                                                 M_FeO + M_FeSiO3 + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                 M_FeO + M_FeSiO3 + M_SiO2 - M_m))) + M_Mg * (
                                                                                         M_Fe * (M_FeO * (
                                                                                         M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                 4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                                 M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                 4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_O * (
                                                                                                 M_FeO + M_FeSiO3 + M_MgO + M_MgSiO3 + M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeO * (
                                                                                         M_MgO * (
                                                                                         M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                         M_FeO * (
                                                                                         4.0 * M_MgO + 4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                         M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                         M_FeO * (
                                                                                         4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_FeSiO3 * (
                                                                                         4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m)) + M_O * (
                                                                                         M_FeO * (
                                                                                         M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                         M_MgO + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                         M_FeO + M_FeSiO3))) + M_MgO * (
                                                                                         -4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                         M_FeO * (
                                                                                         4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                         M_FeO * (M_MgO * (
                                                                                         4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                         M_FeO * (
                                                                                         4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                         4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_O * (
                                                                                         M_MgO * (M_FeO * (
                                                                                         M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                  M_SiO2 - M_m)) + M_MgSiO3 * (
                                                                                         M_FeO * (
                                                                                         M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                         M_SiO2 - M_m))))))

    def dM_MgO_dTc(self, Moles, dKs, dMi_b):
        dM_Mg_er, dM_Si_er, dM_Fe_er, dM_O_er, dM_c_er, dM_MgO_er, dM_SiO2_er, dM_FeO_er, dM_MgSiO3_er, dM_FeSiO3_er, dM_m_er = self.unwrap_Moles(
            dMi_b)
        M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = Moles
        dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs
        return M_MgO * (M_Fe * (M_O * (M_FeO * M_SiO2 * M_m * (-4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_FeSiO3 * (
        M_FeO * (M_SiO2 * (4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er) + M_m * (
        -4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er)) + M_SiO2 * M_m * (-4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er)) + M_MgSiO3 * (
                                       M_FeO * (M_SiO2 * (
                                       -4.0 * dM_FeO_er - 8.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_m * (
                                                4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er + 4.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                       M_SiO2 * (-4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_m * (
                                       -4.0 * dM_FeO_er - 4.0 * dM_MgO_er + 4.0 * dM_SiO2_er)))) + dKFeO_KFeO * (
                                M_Mg * M_O * (
                                M_FeO * (M_MgSiO3 * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_Si * (M_Mg * (
                                M_FeO * (M_MgSiO3 * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
                                M_FeO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                                M_FeO * (-9.0 * M_MgSiO3 - M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                                -9.0 * M_FeO + 2.0 * M_SiO2 + 10.0 * M_m))) + M_MgSiO3 * M_O * (M_FeO * (
                                -3.0 * M_SiO2 - 6.0 * M_m) + M_FeSiO3 * (6.0 * M_SiO2 - 15.0 * M_m))) + M_c * (M_Mg * (
                                M_FeO * (M_MgSiO3 * (M_SiO2 - M_m) + M_O * (
                                M_FeSiO3 + M_MgSiO3 + M_SiO2) - M_SiO2 * M_m) + M_FeSiO3 * (
                                M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * M_O * (
                                                                                                               M_FeO * M_SiO2 + M_FeSiO3 * M_m) + M_Si * (
                                                                                                               M_Mg * (
                                                                                                               M_FeO * (
                                                                                                               4.0 * M_MgSiO3 + M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                               4.0 * M_FeO - M_SiO2 - 3.0 * M_m) + M_O * (
                                                                                                               M_FeO + M_FeSiO3)) + M_MgSiO3 * (
                                                                                                               M_FeO * (
                                                                                                               2.0 * M_SiO2 + 2.0 * M_m) + M_FeSiO3 * (
                                                                                                               -2.0 * M_SiO2 + 6.0 * M_m)))))) + M_FeSiO3 * dKFeSiO3_KFeSiO3 * (
                        M_O * (M_Fe * M_FeO * M_MgSiO3 * (4.0 * M_SiO2 - 4.0 * M_m) + M_Mg * (
                        4.0 * M_Fe * M_SiO2 * M_m + M_FeO * M_MgSiO3 * (4.0 * M_SiO2 - 4.0 * M_m))) + M_Si * (M_Mg * (
                        M_Fe * (M_O * (2.0 * M_SiO2 + 10.0 * M_m) + M_SiO2 * M_m) + M_FeO * (
                        M_MgSiO3 * (M_SiO2 - M_m) + M_O * (9.0 * M_MgSiO3 + 3.0 * M_SiO2 + 6.0 * M_m))) + M_MgSiO3 * (
                                                                                                              M_Fe * (
                                                                                                              M_FeO * (
                                                                                                              M_SiO2 - M_m) + M_O * (
                                                                                                              9.0 * M_FeO + 6.0 * M_SiO2 - 15.0 * M_m)) + M_FeO * M_O * (
                                                                                                              9.0 * M_SiO2 - 9.0 * M_m))) + M_c * (
                        M_Mg * (-M_Fe * M_SiO2 * M_m + M_FeO * (
                        M_MgSiO3 * (-M_SiO2 + M_m) + M_O * (-M_MgSiO3 - M_SiO2))) + M_MgSiO3 * (
                        M_Fe * (M_FeO * (-M_SiO2 + M_m) + M_O * (-M_FeO + M_m)) + M_FeO * M_O * (
                        -M_SiO2 + M_m)) + M_Si * (M_Mg * (M_Fe * (M_O - M_SiO2 - 3.0 * M_m) + M_FeO * (
                        -4.0 * M_MgSiO3 - 2.0 * M_SiO2 - 2.0 * M_m)) + M_MgSiO3 * (
                                                  M_Fe * (-4.0 * M_FeO - 2.0 * M_SiO2 + 6.0 * M_m) + M_FeO * (
                                                  -4.0 * M_SiO2 + 4.0 * M_m))))) + M_Mg * (M_O * (M_Fe * (M_FeSiO3 * (
        M_SiO2 * (-4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_m * (
        -4.0 * dM_FeO_er - 4.0 * dM_MgO_er + 4.0 * dM_SiO2_er)) + M_MgSiO3 * (M_SiO2 * (
        -4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_m * (
                                                                              -4.0 * dM_FeO_er - 4.0 * dM_MgO_er + 4.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                                          -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er)) + M_FeO * M_SiO2 * M_m * (
                                                                                                  -4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_FeSiO3 * (
                                                                                                  M_FeO * (M_SiO2 * (
                                                                                                  4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er) + M_m * (
                                                                                                           -4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er)) + M_SiO2 * M_m * (
                                                                                                  -4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er)) + M_MgSiO3 * (
                                                                                                  M_FeO * (M_SiO2 * (
                                                                                                  -4.0 * dM_FeO_er - 8.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_m * (
                                                                                                           4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er + 4.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                                                  M_SiO2 * (
                                                                                                  -4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_m * (
                                                                                                  -4.0 * dM_FeO_er - 4.0 * dM_MgO_er + 4.0 * dM_SiO2_er)))) + dKMgO_KMgO * (
                                                                                           M_Fe * M_O * (M_FeO * (
                                                                                           M_MgSiO3 * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                         M_FeO * (
                                                                                                         4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_Si * (
                                                                                           M_Fe * (M_FeO * (M_MgSiO3 * (
                                                                                           M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                   M_FeO * (
                                                                                                   M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
                                                                                                   M_FeO * (
                                                                                                   M_SiO2 - 4.0 * M_m) + M_FeSiO3 * (
                                                                                                   9.0 * M_FeO + 4.0 * M_SiO2 - 25.0 * M_m) + M_MgSiO3 * (
                                                                                                   9.0 * M_FeO + 6.0 * M_SiO2 - 15.0 * M_m) - 9.0 * M_SiO2 * M_m)) + M_O * (
                                                                                           -9.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           9.0 * M_SiO2 - 9.0 * M_m) - 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                           M_FeO * (
                                                                                           6.0 * M_SiO2 - 15.0 * M_m) + M_FeSiO3 * (
                                                                                           6.0 * M_SiO2 - 15.0 * M_m)))) + M_c * (
                                                                                           M_Fe * (M_FeO * (M_MgSiO3 * (
                                                                                           -M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                   M_FeO * (
                                                                                                   -M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                                                                                                   M_FeSiO3 * (
                                                                                                   -M_FeO + M_m) + M_MgSiO3 * (
                                                                                                   -M_FeO + M_m) + M_SiO2 * (
                                                                                                   -M_FeO + M_m))) + M_O * (
                                                                                           M_FeSiO3 * (M_FeO * (
                                                                                           -M_SiO2 + M_m) + M_SiO2 * M_m) + M_m * (
                                                                                           M_FeO * M_SiO2 + M_MgSiO3 * (
                                                                                           M_FeO + M_FeSiO3))) + M_Si * (
                                                                                           M_Fe * (M_FeO * (
                                                                                           -M_SiO2 + M_m) + M_FeSiO3 * (
                                                                                                   -4.0 * M_FeO - M_SiO2 + 9.0 * M_m) + M_MgSiO3 * (
                                                                                                   -4.0 * M_FeO - 2.0 * M_SiO2 + 6.0 * M_m) + M_O * (
                                                                                                   -M_FeO - M_FeSiO3 - M_SiO2 + M_m) + 4.0 * M_SiO2 * M_m) + 4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                           M_FeO * (
                                                                                           -2.0 * M_SiO2 + 6.0 * M_m) + M_FeSiO3 * (
                                                                                           -2.0 * M_SiO2 + 6.0 * M_m)) + M_O * (
                                                                                           M_FeO * (
                                                                                           -M_SiO2 + M_m) + M_FeSiO3 * (
                                                                                           -M_SiO2 + M_m)))))) + M_MgSiO3 * dKMgSiO3_KMgSiO3 * (
                        M_O * (M_Fe * (4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                        M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_Mg * (
                               M_FeSiO3 * (M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                               4.0 * M_Fe + 4.0 * M_FeO))) + M_Si * (M_Fe * (
                        M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                        M_FeO * (-M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                        -9.0 * M_FeO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m)) + M_Mg * (M_Fe * (
                        M_O * (2.0 * M_SiO2 + 10.0 * M_m) + M_SiO2 * M_m) + M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (
                        -M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (M_FeO * (2.0 * M_SiO2 + 10.0 * M_m) + M_FeSiO3 * (
                        -9.0 * M_FeO + 2.0 * M_SiO2 + 10.0 * M_m))) + M_O * (9.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                        M_FeO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m))) + M_c * (M_Fe * (
                        -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
                        M_FeSiO3 * (M_FeO - M_m) + M_SiO2 * (M_FeO - M_m))) + M_Mg * (M_FeSiO3 * (
                        M_FeO * (M_O + M_SiO2 - M_m) - M_SiO2 * M_m) + M_SiO2 * M_m * (-M_Fe - M_FeO)) + M_O * (
                                                                                             -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                             M_FeO * (
                                                                                             M_SiO2 - M_m) - M_SiO2 * M_m)) + M_Si * (
                                                                                             M_Fe * (M_FeO * (
                                                                                             M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                     4.0 * M_FeO + M_SiO2 - 9.0 * M_m) + M_O * (
                                                                                                     M_FeO + M_FeSiO3 + M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) - 4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                             M_FeO * (
                                                                                             4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_Mg * (
                                                                                             M_Fe * (
                                                                                             M_O - M_SiO2 - 3.0 * M_m) + M_FeO * (
                                                                                             -M_SiO2 - 3.0 * M_m) + M_FeSiO3 * (
                                                                                             4.0 * M_FeO - M_SiO2 - 3.0 * M_m) + M_O * (
                                                                                             M_FeO + M_FeSiO3)) + M_O * (
                                                                                             M_FeO * (
                                                                                             M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                             M_SiO2 - M_m))))) + M_Si * (
                        M_Fe * (M_FeO * M_SiO2 * M_m * (-dM_MgO_er - dM_MgSiO3_er) + M_FeSiO3 * (M_FeO * (
                        M_SiO2 * (dM_MgO_er + dM_MgSiO3_er) + M_m * (-dM_MgO_er - dM_MgSiO3_er)) + M_SiO2 * M_m * (
                                                                                                 -dM_MgO_er - dM_MgSiO3_er)) + M_MgSiO3 * (
                                M_FeO * (M_SiO2 * (-dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er + 1.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                M_SiO2 * (-1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
                                -1.0 * dM_FeO_er - 1.0 * dM_MgO_er + 1.0 * dM_SiO2_er))) + M_O * (M_FeO * (
                        M_SiO2 * (dM_MgO_er + dM_MgSiO3_er) + M_m * (
                        -4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er)) + M_FeSiO3 * (M_FeO * (
                        9.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                              4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er) + M_m * (
                                                                              -25.0 * dM_MgO_er - 25.0 * dM_MgSiO3_er)) + M_MgSiO3 * (
                                                                                                  M_FeO * (
                                                                                                  -6.0 * dM_FeO_er - 15.0 * dM_FeSiO3_er + 3.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                                                  -6.0 * dM_FeO_er - 15.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 15.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                                  -6.0 * dM_FeO_er - 15.0 * dM_FeSiO3_er - 9.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er) + M_m * (
                                                                                                  6.0 * dM_FeO_er + 15.0 * dM_FeSiO3_er - 9.0 * dM_MgO_er + 9.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                                  -9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er))) + M_Mg * (
                        M_Fe * (M_FeSiO3 * (
                        M_SiO2 * (-1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
                        -1.0 * dM_FeO_er - 1.0 * dM_MgO_er + 1.0 * dM_SiO2_er)) + M_MgSiO3 * (
                                M_SiO2 * (-1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
                                -1.0 * dM_FeO_er - 1.0 * dM_MgO_er + 1.0 * dM_SiO2_er)) + M_O * (M_FeSiO3 * (
                        -6.0 * dM_FeO_er - 15.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 15.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er) + M_MgSiO3 * (
                                                                                                 -6.0 * dM_FeO_er - 15.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 15.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                                 -2 * dM_FeO_er - 5.0 * dM_FeSiO3_er - 2 * dM_MgO_er - 5.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_m * (
                                                                                                 -4.0 * dM_FeO_er - 10.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 10.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er)) + M_FeO * M_SiO2 * M_m * (
                        -1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_FeSiO3 * (M_FeO * (
                        M_SiO2 * (1.0 * dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_m * (
                        -1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er)) + M_SiO2 * M_m * (
                                                                             -1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er)) + M_MgSiO3 * (
                        M_FeO * (
                        M_SiO2 * (-dM_FeO_er - 2.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
                        1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er + 1.0 * dM_SiO2_er)) + M_FeSiO3 * (
                        M_SiO2 * (-1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
                        -1.0 * dM_FeO_er - 1.0 * dM_MgO_er + 1.0 * dM_SiO2_er))) + M_O * (M_FeO * (M_SiO2 * (
                        -3.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 2 * dM_MgO_er - 5.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_m * (
                                                                                                   -6.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 10.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                                          M_FeO * (
                                                                                          -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er) + M_SiO2 * (
                                                                                          -3.0 * dM_FeSiO3_er - 2 * dM_MgO_er - 5.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_m * (
                                                                                          6.0 * dM_FeO_er - 4.0 * dM_MgO_er - 10.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er)) + M_MgSiO3 * (
                                                                                          M_FeO * (
                                                                                          -15.0 * dM_FeO_er - 24.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 15.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                                          -6.0 * dM_FeO_er - 15.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 15.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er)))) + M_O * (
                        M_FeO * M_SiO2 * M_m * (-9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er) + M_FeSiO3 * (M_FeO * (
                        M_SiO2 * (9.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er) + M_m * (
                        -9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er)) + M_SiO2 * M_m * (
                                                                                                     -9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er)) + M_MgSiO3 * (
                        M_FeO * (M_SiO2 * (
                        -9.0 * dM_FeO_er - 18.0 * dM_FeSiO3_er - 9.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er) + M_m * (
                                 9.0 * dM_FeSiO3_er - 9.0 * dM_MgO_er + 9.0 * dM_SiO2_er)) + M_FeSiO3 * (
                        M_SiO2 * (-9.0 * dM_FeSiO3_er - 9.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er) + M_m * (
                        -9.0 * dM_FeO_er - 9.0 * dM_MgO_er + 9.0 * dM_SiO2_er)))) + dKSiO2_KSiO2 * (M_O * (
                        M_Fe * M_MgSiO3 * (
                        M_FeO * (2.0 * M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (-4.0 * M_SiO2 + 10.0 * M_m)) + M_Mg * (M_Fe * (
                        M_FeSiO3 * (-4.0 * M_SiO2 + 10.0 * M_m) + M_MgSiO3 * (
                        -4.0 * M_SiO2 + 10.0 * M_m) + 6.0 * M_SiO2 * M_m) + 6.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                                M_FeO * (
                                                                                                                -6.0 * M_SiO2 + 6.0 * M_m) + 6.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                                M_FeO * (
                                                                                                                -4.0 * M_SiO2 + 10.0 * M_m) + M_FeSiO3 * (
                                                                                                                -4.0 * M_SiO2 + 10.0 * M_m)))) + M_c * (
                                                                                                    M_Mg * (M_Fe * (
                                                                                                    M_FeSiO3 * (
                                                                                                    M_SiO2 - 3.0 * M_m) + M_MgSiO3 * (
                                                                                                    M_SiO2 - 3.0 * M_m) + M_O * (
                                                                                                    M_FeSiO3 + M_MgSiO3 + M_SiO2) - 2.0 * M_SiO2 * M_m) - 2.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                            M_FeO * (
                                                                                                            2.0 * M_SiO2 - 2.0 * M_m) - 2.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                            M_FeO * (
                                                                                                            M_SiO2 - 3.0 * M_m) + M_FeSiO3 * (
                                                                                                            M_SiO2 - 3.0 * M_m)) + M_O * (
                                                                                                            M_MgSiO3 * (
                                                                                                            M_FeO + M_FeSiO3) + M_SiO2 * (
                                                                                                            M_FeO + M_FeSiO3))) + M_MgSiO3 * (
                                                                                                    M_Fe * (M_FeO * (
                                                                                                    -M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                            M_SiO2 - 3.0 * M_m) + M_O * (
                                                                                                            M_FeO + M_FeSiO3 + M_SiO2 - M_m)) + M_O * (
                                                                                                    M_FeO * (
                                                                                                    M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                    M_SiO2 - M_m)))))) + M_c * (
                        M_Fe * (M_FeO * M_SiO2 * M_m * (dM_MgO_er + dM_MgSiO3_er) + M_FeSiO3 * (M_FeO * (
                        M_SiO2 * (-dM_MgO_er - dM_MgSiO3_er) + M_m * (dM_MgO_er + dM_MgSiO3_er)) + M_SiO2 * M_m * (
                                                                                                dM_MgO_er + dM_MgSiO3_er)) + M_MgSiO3 * (
                                M_FeO * (M_SiO2 * (
                                dM_FeO_er + 2.0 * dM_FeSiO3_er + 1.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er) + M_m * (
                                         -dM_FeSiO3_er + dM_MgO_er - dM_SiO2_er)) + M_FeSiO3 * (
                                M_SiO2 * (dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                                dM_FeO_er + dM_MgO_er - 1.0 * dM_SiO2_er))) + M_O * (M_FeSiO3 * (
                        M_FeO * (-dM_MgO_er - dM_MgSiO3_er) + M_m * (dM_MgO_er + dM_MgSiO3_er)) + M_MgSiO3 * (M_FeO * (
                        dM_FeSiO3_er - dM_MgO_er + dM_SiO2_er) + M_FeSiO3 * (
                                                                                                              dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_SiO2 * (
                                                                                                              dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                                                                                                              -dM_FeSiO3_er + dM_MgO_er - dM_SiO2_er)) + M_SiO2 * (
                                                                                     M_FeO * (
                                                                                     -dM_MgO_er - dM_MgSiO3_er) + M_m * (
                                                                                     dM_MgO_er + dM_MgSiO3_er)))) + M_Mg * (
                        M_Fe * (M_FeSiO3 * (M_SiO2 * (dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                        dM_FeO_er + dM_MgO_er - 1.0 * dM_SiO2_er)) + M_MgSiO3 * (
                                M_SiO2 * (dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                                dM_FeO_er + dM_MgO_er - 1.0 * dM_SiO2_er)) + M_O * (
                                M_FeSiO3 * (dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_MgSiO3 * (
                                dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_SiO2 * (
                                dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_SiO2 * M_m * (
                                dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er)) + M_FeO * M_SiO2 * M_m * (
                        dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_FeSiO3 * (M_FeO * (
                        M_SiO2 * (-dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_m * (
                        dM_MgO_er + 1.0 * dM_MgSiO3_er)) + M_SiO2 * M_m * (
                                                                      dM_MgO_er + 1.0 * dM_MgSiO3_er)) + M_MgSiO3 * (
                        M_FeO * (M_SiO2 * (dM_FeO_er + 2 * dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                        -1.0 * dM_FeSiO3_er + dM_MgO_er - 1.0 * dM_SiO2_er)) + M_FeSiO3 * (
                        M_SiO2 * (dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                        dM_FeO_er + dM_MgO_er - 1.0 * dM_SiO2_er))) + M_O * (
                        M_FeO * M_SiO2 * (dM_FeO_er + 2 * dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_FeSiO3 * (
                        M_FeO * (dM_FeO_er + dM_FeSiO3_er) + M_SiO2 * (
                        dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_MgSiO3 * (
                        M_FeO * (dM_FeO_er + 2 * dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_FeSiO3 * (
                        dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)))) + M_O * (
                        M_FeO * M_SiO2 * M_m * (dM_MgO_er + dM_MgSiO3_er) + M_FeSiO3 * (M_FeO * (
                        M_SiO2 * (-dM_MgO_er - dM_MgSiO3_er) + M_m * (dM_MgO_er + dM_MgSiO3_er)) + M_SiO2 * M_m * (
                                                                                        dM_MgO_er + dM_MgSiO3_er)) + M_MgSiO3 * (
                        M_FeO * (M_SiO2 * (dM_FeO_er + 2 * dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                        -dM_FeSiO3_er + dM_MgO_er - dM_SiO2_er)) + M_FeSiO3 * (
                        M_SiO2 * (dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                        dM_FeO_er + dM_MgO_er - dM_SiO2_er)))) + M_Si * (M_Fe * (
                        M_FeO * (M_SiO2 * (-dM_MgO_er - dM_MgSiO3_er) + M_m * (dM_MgO_er + dM_MgSiO3_er)) + M_FeSiO3 * (
                        M_FeO * (-4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_SiO2 * (-dM_MgO_er - dM_MgSiO3_er) + M_m * (
                        9.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er)) + M_MgSiO3 * (M_FeO * (
                        2.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                             2.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_SiO2 * (
                                                                             2.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
                                                                             -2.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er - 4.0 * dM_SiO2_er)) + M_O * (
                        M_FeO * (-dM_MgO_er - dM_MgSiO3_er) + M_FeSiO3 * (-dM_MgO_er - dM_MgSiO3_er) + M_SiO2 * (
                        -dM_MgO_er - dM_MgSiO3_er) + M_m * (dM_MgO_er + dM_MgSiO3_er)) + M_SiO2 * M_m * (
                        4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er)) + M_FeO * M_SiO2 * M_m * (
                                                                         4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er) + M_FeSiO3 * (
                                                                         M_FeO * (M_SiO2 * (
                                                                         -4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_m * (
                                                                                  4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er)) + M_SiO2 * M_m * (
                                                                         4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er)) + M_Mg * (
                                                                         M_Fe * (M_FeSiO3 * (
                                                                         2.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_MgSiO3 * (
                                                                                 2.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                 dM_FeO_er + 3.0 * dM_FeSiO3_er + dM_MgO_er + 3.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_m * (
                                                                                 dM_FeO_er + 3.0 * dM_FeSiO3_er + dM_MgO_er + 3.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er)) + M_FeO * (
                                                                         M_SiO2 * (
                                                                         2 * dM_FeO_er + 4.0 * dM_FeSiO3_er + dM_MgO_er + 3.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_m * (
                                                                         2.0 * dM_FeSiO3_er + dM_MgO_er + 3.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                         M_FeO * (
                                                                         4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er) + M_SiO2 * (
                                                                         2.0 * dM_FeSiO3_er + dM_MgO_er + 3.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_m * (
                                                                         -2.0 * dM_FeO_er + dM_MgO_er + 3.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er)) + M_MgSiO3 * (
                                                                         M_FeO * (
                                                                         6.0 * dM_FeO_er + 10.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                         2.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er)) + M_O * (
                                                                         M_FeO * (
                                                                         dM_FeO_er + dM_FeSiO3_er) + M_FeSiO3 * (
                                                                         dM_FeO_er + dM_FeSiO3_er))) + M_MgSiO3 * (
                                                                         M_FeO * (M_SiO2 * (
                                                                         4.0 * dM_FeO_er + 8.0 * dM_FeSiO3_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
                                                                                  -4.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er - 4.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                         M_SiO2 * (
                                                                         4.0 * dM_FeSiO3_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
                                                                         4.0 * dM_FeO_er + 4.0 * dM_MgO_er - 4.0 * dM_SiO2_er))) + M_O * (
                                                                         M_FeO * (
                                                                         M_SiO2 * (-dM_MgO_er - dM_MgSiO3_er) + M_m * (
                                                                         dM_MgO_er + dM_MgSiO3_er)) + M_FeSiO3 * (
                                                                         M_SiO2 * (-dM_MgO_er - dM_MgSiO3_er) + M_m * (
                                                                         dM_MgO_er + dM_MgSiO3_er)))))) / (M_O * (
        M_Fe * (M_MgO * (4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
        M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                M_FeO * (M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m))) + M_Mg * (M_Fe * (M_FeSiO3 * (
        M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                     M_FeO * (
                                                                                     -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                                                                                     -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                                                                                     4.0 * M_FeO + 4.0 * M_MgO)) + M_MgO * (
                                                                             4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                             M_FeO * (
                                                                             -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                             M_FeO * (M_MgO * (
                                                                             -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                             M_FeO * (
                                                                             -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                                                                             -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)))) + M_Si * (
                                                                                                           M_Fe * (
                                                                                                           M_MgO * (
                                                                                                           M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           -M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           M_MgO * (
                                                                                                           -M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           -M_SiO2 + M_m) + M_MgO * (
                                                                                                           -M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (
                                                                                                           M_MgO * (
                                                                                                           M_FeO * (
                                                                                                           -M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                                                                                                           -9.0 * M_FeO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           -9.0 * M_MgO - M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                                                                                                           -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                                                                                                           -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m))) + M_Mg * (
                                                                                                           M_Fe * (
                                                                                                           M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           -M_SiO2 + M_m) + M_MgO * (
                                                                                                           -M_SiO2 + M_m) + M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           -M_SiO2 + M_m) + M_MgO * (
                                                                                                           -M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                                                                                                           M_FeO * (
                                                                                                           -M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                                                                                                           -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                                                                                                           -M_SiO2 + 4.0 * M_m) + M_MgSiO3 * (
                                                                                                           -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                                                                                                           M_FeO + M_MgO)) + M_MgO * (
                                                                                                           M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           -M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           M_MgO * (
                                                                                                           -M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           -M_SiO2 + M_m) + M_MgO * (
                                                                                                           -M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (
                                                                                                           M_FeO * (
                                                                                                           M_MgO * (
                                                                                                           -M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           -9.0 * M_MgO - 9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (
                                                                                                           -M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           -9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_FeSiO3 * (
                                                                                                           -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m)))) + M_O * (
                                                                                                           M_MgO * (
                                                                                                           9.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           M_MgO * (
                                                                                                           -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           -9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (
                                                                                                           -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)))) + M_c * (
                                                                                                           M_Fe * (
                                                                                                           M_MgO * (
                                                                                                           -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           M_MgO * (
                                                                                                           M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           M_SiO2 - M_m) + M_MgO * (
                                                                                                           M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
                                                                                                           M_MgO * (
                                                                                                           M_FeSiO3 * (
                                                                                                           M_FeO - M_m) + M_SiO2 * (
                                                                                                           M_FeO - M_m)) + M_MgSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           M_MgO + M_SiO2) + M_FeSiO3 * (
                                                                                                           M_FeO + M_MgO - M_m) + M_MgO * (
                                                                                                           M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Mg * (
                                                                                                           M_Fe * (
                                                                                                           M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           M_SiO2 - M_m) + M_MgO * (
                                                                                                           M_SiO2 - M_m) - M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           M_SiO2 - M_m) + M_MgO * (
                                                                                                           M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
                                                                                                           M_FeSiO3 * (
                                                                                                           M_FeO + M_MgO - M_m) + M_MgSiO3 * (
                                                                                                           M_FeO + M_MgO - M_m) + M_SiO2 * (
                                                                                                           M_FeO + M_MgO - M_m)) + M_SiO2 * M_m * (
                                                                                                           -M_FeO - M_MgO)) + M_MgO * (
                                                                                                           -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           M_MgO * (
                                                                                                           M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           M_SiO2 - M_m) + M_MgO * (
                                                                                                           M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
                                                                                                           M_FeO * M_SiO2 * (
                                                                                                           M_MgO - M_m) + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           M_MgO + M_SiO2 - M_m) + M_SiO2 * (
                                                                                                           M_MgO - M_m)) + M_MgSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           M_MgO - M_m) + M_FeSiO3 * (
                                                                                                           M_FeO + M_MgO - M_m)))) + M_O * (
                                                                                                           M_MgO * (
                                                                                                           -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           M_MgO * (
                                                                                                           M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           M_SiO2 - M_m) + M_MgO * (
                                                                                                           M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Si * (
                                                                                                           M_Fe * (
                                                                                                           M_MgO * (
                                                                                                           M_FeO * (
                                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                           4.0 * M_FeO + M_SiO2 - 9.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           4.0 * M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                           4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_O * (
                                                                                                           M_MgO * (
                                                                                                           M_FeO + M_FeSiO3 + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                           M_FeO + M_FeSiO3 + M_SiO2 - M_m))) + M_Mg * (
                                                                                                           M_Fe * (
                                                                                                           M_FeO * (
                                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                           4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                                           M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                           4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_O * (
                                                                                                           M_FeO + M_FeSiO3 + M_MgO + M_MgSiO3 + M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeO * (
                                                                                                           M_MgO * (
                                                                                                           M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           4.0 * M_MgO + 4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                                           M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_FeSiO3 * (
                                                                                                           4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m)) + M_O * (
                                                                                                           M_FeO * (
                                                                                                           M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                           M_MgO + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                           M_FeO + M_FeSiO3))) + M_MgO * (
                                                                                                           -4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           M_MgO * (
                                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_O * (
                                                                                                           M_MgO * (
                                                                                                           M_FeO * (
                                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                           M_SiO2 - M_m)) + M_MgSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                           M_SiO2 - M_m))))))

    def dM_SiO2_dTc(self, Moles, dKs, dMi_b):
        dM_Mg_er, dM_Si_er, dM_Fe_er, dM_O_er, dM_c_er, dM_MgO_er, dM_SiO2_er, dM_FeO_er, dM_MgSiO3_er, dM_FeSiO3_er, dM_m_er = self.unwrap_Moles(
            dMi_b)
        M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = Moles
        dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs
        return M_SiO2 * (M_Fe * (M_O * (M_MgO * (
        M_FeO * M_m * (-4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_FeSiO3 * (
        M_FeO * (-4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_m * (
        4.0 * dM_FeO_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er))) + M_MgSiO3 * (M_FeO * (
        M_MgO * (-4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_m * (
        -4.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er - 4.0 * dM_SiO2_er)) + M_FeSiO3 * (M_FeO * (
        -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                 -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_m * (
                                                                                 4.0 * dM_FeO_er + 4.0 * dM_MgO_er - 4.0 * dM_SiO2_er)))) + dKFeO_KFeO * (
                                 M_Mg * M_O * M_m * (-4.0 * M_FeO * M_MgSiO3 + 4.0 * M_FeSiO3 * M_MgO) + M_Si * (
                                 M_Mg * (M_O * (M_FeO * (-6.0 * M_MgSiO3 + 6.0 * M_m) + M_FeSiO3 * (
                                 -9.0 * M_FeO - 3.0 * M_MgO + 15.0 * M_m)) + M_m * (
                                         -M_FeO * M_MgSiO3 + M_FeSiO3 * M_MgO)) + M_O * (
                                 M_MgO * (6.0 * M_FeO * M_m + M_FeSiO3 * (-9.0 * M_FeO + 15.0 * M_m)) + M_MgSiO3 * (
                                 M_FeO * (-9.0 * M_MgO + 6.0 * M_m) + M_FeSiO3 * (
                                 -9.0 * M_FeO - 9.0 * M_MgO + 15.0 * M_m)))) + M_c * (M_Mg * (
                                 M_FeO * M_MgSiO3 * M_m + M_FeSiO3 * (
                                 -M_MgO * M_m + M_O * (M_FeO + M_MgO - M_m))) + M_O * (M_FeSiO3 * M_MgO * (
                                 M_FeO - M_m) + M_MgSiO3 * (M_FeO * M_MgO + M_FeSiO3 * (
                                 M_FeO + M_MgO - M_m))) + M_Si * (M_Mg * (
                                 M_FeO * (2.0 * M_MgSiO3 - 2.0 * M_m) + M_FeSiO3 * (
                                 4.0 * M_FeO + 2.0 * M_MgO - 6.0 * M_m) + M_O * (M_FeO + M_FeSiO3)) + M_MgO * (
                                                                  -2.0 * M_FeO * M_m + M_FeSiO3 * (
                                                                  4.0 * M_FeO - 6.0 * M_m)) + M_MgSiO3 * (
                                                                  M_FeO * (4.0 * M_MgO - 2.0 * M_m) + M_FeSiO3 * (
                                                                  4.0 * M_FeO + 4.0 * M_MgO - 6.0 * M_m)) + M_O * (
                                                                  M_MgO * (M_FeO + M_FeSiO3) + M_MgSiO3 * (
                                                                  M_FeO + M_FeSiO3)))))) + M_FeSiO3 * dKFeSiO3_KFeSiO3 * (
                         M_O * M_m * (M_Fe * M_FeO * (4.0 * M_MgO + 4.0 * M_MgSiO3) + M_Mg * (
                         M_Fe * (4.0 * M_FeO + 4.0 * M_MgO) + M_FeO * (4.0 * M_MgO + 4.0 * M_MgSiO3))) + M_Si * (
                         M_Fe * (M_FeO * M_m * (M_MgO + M_MgSiO3) + M_O * (
                         M_MgO * (-3.0 * M_FeO + 15.0 * M_m) + M_MgSiO3 * (
                         -3.0 * M_FeO - 9.0 * M_MgO + 15.0 * M_m))) + M_FeO * M_O * M_m * (
                         9.0 * M_MgO + 9.0 * M_MgSiO3) + M_Mg * (
                         M_Fe * (M_O * (-3.0 * M_FeO - 3.0 * M_MgO + 15.0 * M_m) + M_m * (M_FeO + M_MgO)) + M_FeO * (
                         M_O * (-3.0 * M_MgO + 6.0 * M_MgSiO3 + 9.0 * M_m) + M_m * (M_MgO + M_MgSiO3)))) + M_c * (
                         M_Fe * (M_FeO * M_m * (-M_MgO - M_MgSiO3) + M_O * (
                         M_MgO * (M_FeO - M_m) + M_MgSiO3 * (M_FeO + M_MgO - M_m))) + M_FeO * M_O * M_m * (
                         -M_MgO - M_MgSiO3) + M_Mg * (
                         M_Fe * (M_O * (M_FeO + M_MgO - M_m) + M_m * (-M_FeO - M_MgO)) + M_FeO * (
                         M_O * (M_MgO - M_m) + M_m * (-M_MgO - M_MgSiO3))) + M_Si * (M_Fe * (
                         M_MgO * (2.0 * M_FeO - 6.0 * M_m) + M_MgSiO3 * (
                         2.0 * M_FeO + 4.0 * M_MgO - 6.0 * M_m) + M_O * (M_MgO + M_MgSiO3)) + M_FeO * M_m * (
                                                                                     -4.0 * M_MgO - 4.0 * M_MgSiO3) + M_Mg * (
                                                                                     M_Fe * (
                                                                                     2.0 * M_FeO + 2.0 * M_MgO + M_O - 6.0 * M_m) + M_FeO * (
                                                                                     2.0 * M_MgO - 2.0 * M_MgSiO3 - 4.0 * M_m))))) + M_Mg * (
                         M_O * (M_Fe * (M_FeSiO3 * (M_FeO * (
                         -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_MgO * (
                                                    -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_m * (
                                                    4.0 * dM_FeO_er + 4.0 * dM_MgO_er - 4.0 * dM_SiO2_er)) + M_MgSiO3 * (
                                        M_FeO * (
                                        -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_MgO * (
                                        -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_m * (
                                        4.0 * dM_FeO_er + 4.0 * dM_MgO_er - 4.0 * dM_SiO2_er)) + M_m * (M_FeO * (
                         -4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_MgO * (
                                                                                                        -4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er))) + M_MgO * (
                                M_FeO * M_m * (
                                -4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_FeSiO3 * (M_FeO * (
                                -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_m * (
                                                                                                           4.0 * dM_FeO_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er))) + M_MgSiO3 * (
                                M_FeO * (M_MgO * (
                                -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_m * (
                                         -4.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er - 4.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                M_FeO * (
                                -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_MgO * (
                                -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_m * (
                                4.0 * dM_FeO_er + 4.0 * dM_MgO_er - 4.0 * dM_SiO2_er)))) + dKMgO_KMgO * (
                         M_Fe * M_O * M_m * (4.0 * M_FeO * M_MgSiO3 - 4.0 * M_FeSiO3 * M_MgO) + M_Si * (M_Fe * (M_O * (
                         M_MgO * (-6.0 * M_FeSiO3 + 6.0 * M_m) + M_MgSiO3 * (
                         -3.0 * M_FeO - 9.0 * M_MgO + 15.0 * M_m)) + M_m * (
                                                                                                                M_FeO * M_MgSiO3 - M_FeSiO3 * M_MgO)) + M_O * (
                                                                                                        M_MgO * (
                                                                                                        6.0 * M_FeO * M_m + M_FeSiO3 * (
                                                                                                        -9.0 * M_FeO + 6.0 * M_m)) + M_MgSiO3 * (
                                                                                                        M_FeO * (
                                                                                                        -9.0 * M_MgO + 15.0 * M_m) + M_FeSiO3 * (
                                                                                                        -9.0 * M_FeO - 9.0 * M_MgO + 15.0 * M_m)))) + M_c * (
                         M_Fe * (
                         M_FeSiO3 * M_MgO * M_m + M_MgSiO3 * (-M_FeO * M_m + M_O * (M_FeO + M_MgO - M_m))) + M_O * (
                         M_FeO * M_FeSiO3 * M_MgO + M_MgSiO3 * (
                         M_FeO * (M_MgO - M_m) + M_FeSiO3 * (M_FeO + M_MgO - M_m))) + M_Si * (M_Fe * (
                         M_MgO * (2.0 * M_FeSiO3 - 2.0 * M_m) + M_MgSiO3 * (
                         2.0 * M_FeO + 4.0 * M_MgO - 6.0 * M_m) + M_O * (M_MgO + M_MgSiO3)) + M_MgO * (
                                                                                              -2.0 * M_FeO * M_m + M_FeSiO3 * (
                                                                                              4.0 * M_FeO - 2.0 * M_m)) + M_MgSiO3 * (
                                                                                              M_FeO * (
                                                                                              4.0 * M_MgO - 6.0 * M_m) + M_FeSiO3 * (
                                                                                              4.0 * M_FeO + 4.0 * M_MgO - 6.0 * M_m)) + M_O * (
                                                                                              M_MgO * (
                                                                                              M_FeO + M_FeSiO3) + M_MgSiO3 * (
                                                                                              M_FeO + M_FeSiO3)))))) + M_MgSiO3 * dKMgSiO3_KMgSiO3 * (
                         M_O * M_m * (M_Fe * M_MgO * (4.0 * M_FeO + 4.0 * M_FeSiO3) + M_Mg * (
                         M_Fe * (4.0 * M_FeO + 4.0 * M_MgO) + M_MgO * (4.0 * M_FeO + 4.0 * M_FeSiO3))) + M_Si * (
                         M_Mg * (M_Fe * (
                         M_O * (-3.0 * M_FeO - 3.0 * M_MgO + 15.0 * M_m) + M_m * (M_FeO + M_MgO)) + M_MgO * M_m * (
                                 M_FeO + M_FeSiO3) + M_O * (M_FeO * (-3.0 * M_MgO + 15.0 * M_m) + M_FeSiO3 * (
                         -9.0 * M_FeO - 3.0 * M_MgO + 15.0 * M_m))) + M_MgO * (M_Fe * (
                         M_O * (-3.0 * M_FeO + 6.0 * M_FeSiO3 + 9.0 * M_m) + M_m * (M_FeO + M_FeSiO3)) + M_O * M_m * (
                                                                               9.0 * M_FeO + 9.0 * M_FeSiO3))) + M_c * (
                         M_Mg * (M_Fe * (M_O * (M_FeO + M_MgO - M_m) + M_m * (-M_FeO - M_MgO)) + M_MgO * M_m * (
                         -M_FeO - M_FeSiO3) + M_O * (
                                 M_FeO * (M_MgO - M_m) + M_FeSiO3 * (M_FeO + M_MgO - M_m))) + M_MgO * (
                         M_Fe * (M_O * (M_FeO - M_m) + M_m * (-M_FeO - M_FeSiO3)) + M_O * M_m * (
                         -M_FeO - M_FeSiO3)) + M_Si * (M_Mg * (
                         M_Fe * (2.0 * M_FeO + 2.0 * M_MgO + M_O - 6.0 * M_m) + M_FeO * (
                         2.0 * M_MgO - 6.0 * M_m) + M_FeSiO3 * (4.0 * M_FeO + 2.0 * M_MgO - 6.0 * M_m) + M_O * (
                         M_FeO + M_FeSiO3)) + M_MgO * (M_Fe * (2.0 * M_FeO - 2.0 * M_FeSiO3 - 4.0 * M_m) + M_m * (
                         -4.0 * M_FeO - 4.0 * M_FeSiO3))))) + M_Si * (M_Fe * (M_MgO * (
        M_FeO * M_m * (-1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_FeSiO3 * (
        M_FeO * (-1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_m * (
        1.0 * dM_FeO_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er))) + M_MgSiO3 * (M_FeO * (
        M_MgO * (-1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_m * (
        -1.0 * dM_FeSiO3_er + 1.0 * dM_MgO_er - 1.0 * dM_SiO2_er)) + M_FeSiO3 * (M_FeO * (
        -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                 -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_m * (
                                                                                 1.0 * dM_FeO_er + 1.0 * dM_MgO_er - 1.0 * dM_SiO2_er))) + M_O * (
                                                                              M_MgO * (M_FeO * (
                                                                              2 * dM_FeO_er + 5.0 * dM_FeSiO3_er + 2 * dM_MgO_er + 5.0 * dM_MgSiO3_er + 3.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                                       -4.0 * dM_FeO_er - 10.0 * dM_FeSiO3_er - 10.0 * dM_MgO_er - 16.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er) + M_m * (
                                                                                       -6.0 * dM_FeO_er - 15.0 * dM_FeSiO3_er - 9.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er)) + M_MgSiO3 * (
                                                                              M_FeO * (
                                                                              2 * dM_FeO_er + 5.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er + 3.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                              -4.0 * dM_FeO_er - 10.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 10.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er) + M_MgO * (
                                                                              -9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er) + M_m * (
                                                                              -6.0 * dM_FeO_er - 15.0 * dM_FeSiO3_er + 9.0 * dM_MgO_er - 9.0 * dM_SiO2_er)))) + M_Mg * (
                                                                      M_Fe * (M_FeSiO3 * (M_FeO * (
                                                                      -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                          -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_m * (
                                                                                          1.0 * dM_FeO_er + 1.0 * dM_MgO_er - 1.0 * dM_SiO2_er)) + M_MgSiO3 * (
                                                                              M_FeO * (
                                                                              -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_MgO * (
                                                                              -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_m * (
                                                                              1.0 * dM_FeO_er + 1.0 * dM_MgO_er - 1.0 * dM_SiO2_er)) + M_O * (
                                                                              M_FeO * (
                                                                              2 * dM_FeO_er + 5.0 * dM_FeSiO3_er + 2 * dM_MgO_er + 5.0 * dM_MgSiO3_er + 3.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                              -4.0 * dM_FeO_er - 10.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 10.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er) + M_MgO * (
                                                                              2 * dM_FeO_er + 5.0 * dM_FeSiO3_er + 2 * dM_MgO_er + 5.0 * dM_MgSiO3_er + 3.0 * dM_SiO2_er) + M_MgSiO3 * (
                                                                              -4.0 * dM_FeO_er - 10.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 10.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er) + M_m * (
                                                                              -6.0 * dM_FeO_er - 15.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 15.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er)) + M_m * (
                                                                              M_FeO * (
                                                                              -1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_MgO * (
                                                                              -1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er))) + M_MgO * (
                                                                      M_FeO * M_m * (
                                                                      -1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                      M_FeO * (
                                                                      -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_m * (
                                                                      1.0 * dM_FeO_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er))) + M_MgSiO3 * (
                                                                      M_FeO * (M_MgO * (
                                                                      -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_m * (
                                                                               -1.0 * dM_FeSiO3_er + 1.0 * dM_MgO_er - 1.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                      M_FeO * (
                                                                      -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_MgO * (
                                                                      -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_m * (
                                                                      1.0 * dM_FeO_er + 1.0 * dM_MgO_er - 1.0 * dM_SiO2_er))) + M_O * (
                                                                      M_FeO * (M_MgO * (
                                                                      2 * dM_FeO_er + 5.0 * dM_FeSiO3_er + 2 * dM_MgO_er + 5.0 * dM_MgSiO3_er + 3.0 * dM_SiO2_er) + M_m * (
                                                                               -9.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 15.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                      M_FeO * (
                                                                      -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er) + M_MgO * (
                                                                      -1.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er + 2 * dM_MgO_er + 5.0 * dM_MgSiO3_er + 3.0 * dM_SiO2_er) + M_m * (
                                                                      9.0 * dM_FeO_er - 6.0 * dM_MgO_er - 15.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er)) + M_MgSiO3 * (
                                                                      M_FeO * (
                                                                      -10.0 * dM_FeO_er - 16.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 10.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                      -4.0 * dM_FeO_er - 10.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 10.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er)))) + M_O * (
                                                                      M_MgO * (M_FeO * M_m * (
                                                                      -9.0 * dM_FeSiO3_er - 9.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                               M_FeO * (
                                                                               -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er - 9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er) + M_m * (
                                                                               9.0 * dM_FeO_er - 9.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er))) + M_MgSiO3 * (
                                                                      M_FeO * (M_MgO * (
                                                                      -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er - 9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er) + M_m * (
                                                                               -9.0 * dM_FeSiO3_er + 9.0 * dM_MgO_er - 9.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                      M_FeO * (
                                                                      -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er - 9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er) + M_MgO * (
                                                                      -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er - 9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er) + M_m * (
                                                                      9.0 * dM_FeO_er + 9.0 * dM_MgO_er - 9.0 * dM_SiO2_er)))) + dKSiO2_KSiO2 * (
                                                                      M_O * (M_Fe * (M_MgO * (
                                                                      -4.0 * M_FeO * M_m + M_FeSiO3 * (
                                                                      6.0 * M_FeO - 10.0 * M_m)) + M_MgSiO3 * (M_FeO * (
                                                                      6.0 * M_MgO - 4.0 * M_m) + M_FeSiO3 * (
                                                                                                               6.0 * M_FeO + 6.0 * M_MgO - 10.0 * M_m))) + M_Mg * (
                                                                             M_Fe * (M_FeSiO3 * (
                                                                             6.0 * M_FeO + 6.0 * M_MgO - 10.0 * M_m) + M_MgSiO3 * (
                                                                                     6.0 * M_FeO + 6.0 * M_MgO - 10.0 * M_m) + M_m * (
                                                                                     -4.0 * M_FeO - 4.0 * M_MgO)) + M_MgO * (
                                                                             -4.0 * M_FeO * M_m + M_FeSiO3 * (
                                                                             6.0 * M_FeO - 4.0 * M_m)) + M_MgSiO3 * (
                                                                             M_FeO * (
                                                                             6.0 * M_MgO - 10.0 * M_m) + M_FeSiO3 * (
                                                                             6.0 * M_FeO + 6.0 * M_MgO - 10.0 * M_m)))) + M_c * (
                                                                      M_Fe * (M_MgO * (M_FeO * M_m + M_FeSiO3 * (
                                                                      -2.0 * M_FeO + 3.0 * M_m)) + M_MgSiO3 * (M_FeO * (
                                                                      -2.0 * M_MgO + M_m) + M_FeSiO3 * (
                                                                                                               -2.0 * M_FeO - 2.0 * M_MgO + 3.0 * M_m)) + M_O * (
                                                                              M_MgO * (-M_FeO + M_m) + M_MgSiO3 * (
                                                                              -M_FeO + M_m))) + M_Mg * (M_Fe * (
                                                                      M_FeSiO3 * (
                                                                      -2.0 * M_FeO - 2.0 * M_MgO + 3.0 * M_m) + M_MgSiO3 * (
                                                                      -2.0 * M_FeO - 2.0 * M_MgO + 3.0 * M_m) + M_O * (
                                                                      -M_FeO - M_MgO + M_m) + M_m * (
                                                                      M_FeO + M_MgO)) + M_MgO * (
                                                                                                        M_FeO * M_m + M_FeSiO3 * (
                                                                                                        -2.0 * M_FeO + M_m)) + M_MgSiO3 * (
                                                                                                        M_FeO * (
                                                                                                        -2.0 * M_MgO + 3.0 * M_m) + M_FeSiO3 * (
                                                                                                        -2.0 * M_FeO - 2.0 * M_MgO + 3.0 * M_m)) + M_O * (
                                                                                                        M_FeO * (
                                                                                                        -M_MgO + M_m) + M_FeSiO3 * (
                                                                                                        -M_MgO + M_m))) + M_O * M_m * (
                                                                      M_MgO * (M_FeO + M_FeSiO3) + M_MgSiO3 * (
                                                                      M_FeO + M_FeSiO3))))) + M_c * (M_Fe * (M_MgO * (
        M_FeO * M_m * (dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_FeSiO3 * (
        M_FeO * (dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_m * (
        -dM_FeO_er + 1.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er))) + M_MgSiO3 * (M_FeO * (
        M_MgO * (dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_m * (
        dM_FeSiO3_er - dM_MgO_er + dM_SiO2_er)) + M_FeSiO3 * (M_FeO * (
        dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_MgO * (
                                                              dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_m * (
                                                              -dM_FeO_er - dM_MgO_er + 1.0 * dM_SiO2_er))) + M_O * (
                                                                                                             M_MgO * (
                                                                                                             M_FeO * (
                                                                                                             -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                                                             dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_MgSiO3 * (
                                                                                                             M_FeO * (
                                                                                                             -dM_FeSiO3_er + dM_MgO_er - dM_SiO2_er) + M_MgO * (
                                                                                                             dM_MgO_er + dM_MgSiO3_er) + M_m * (
                                                                                                             dM_FeSiO3_er - dM_MgO_er + dM_SiO2_er)))) + M_Mg * (
                                                                                                     M_Fe * (
                                                                                                     M_FeSiO3 * (
                                                                                                     M_FeO * (
                                                                                                     dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                                     dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_m * (
                                                                                                     -dM_FeO_er - dM_MgO_er + 1.0 * dM_SiO2_er)) + M_MgSiO3 * (
                                                                                                     M_FeO * (
                                                                                                     dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                                     dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_m * (
                                                                                                     -dM_FeO_er - dM_MgO_er + 1.0 * dM_SiO2_er)) + M_O * (
                                                                                                     M_FeO * (
                                                                                                     -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_MgO * (
                                                                                                     -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                                                     dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_m * (
                                                                                                     M_FeO * (
                                                                                                     dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_MgO * (
                                                                                                     dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er))) + M_MgO * (
                                                                                                     M_FeO * M_m * (
                                                                                                     dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_FeSiO3 * (
                                                                                                     M_FeO * (
                                                                                                     dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_m * (
                                                                                                     -dM_FeO_er + dM_MgSiO3_er + dM_SiO2_er))) + M_MgSiO3 * (
                                                                                                     M_FeO * (M_MgO * (
                                                                                                     dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_m * (
                                                                                                              1.0 * dM_FeSiO3_er - dM_MgO_er + 1.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                                                     M_FeO * (
                                                                                                     dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                                     dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_m * (
                                                                                                     -dM_FeO_er - dM_MgO_er + 1.0 * dM_SiO2_er))) + M_O * (
                                                                                                     M_FeO * (M_MgO * (
                                                                                                     -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                                                              dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_FeSiO3 * (
                                                                                                     M_FeO * (
                                                                                                     dM_FeO_er + dM_FeSiO3_er) + M_MgO * (
                                                                                                     dM_FeO_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                                                     -dM_FeO_er + dM_MgSiO3_er + dM_SiO2_er)))) + M_O * (
                                                                                                     M_MgO * (
                                                                                                     M_FeO * M_m * (
                                                                                                     dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_FeSiO3 * (
                                                                                                     M_FeO * (
                                                                                                     dM_FeO_er + dM_FeSiO3_er + dM_MgO_er + dM_MgSiO3_er) + M_m * (
                                                                                                     -dM_FeO_er + dM_MgSiO3_er + dM_SiO2_er))) + M_MgSiO3 * (
                                                                                                     M_FeO * (M_MgO * (
                                                                                                     dM_FeO_er + dM_FeSiO3_er + dM_MgO_er + dM_MgSiO3_er) + M_m * (
                                                                                                              dM_FeSiO3_er - dM_MgO_er + dM_SiO2_er)) + M_FeSiO3 * (
                                                                                                     M_FeO * (
                                                                                                     dM_FeO_er + dM_FeSiO3_er + dM_MgO_er + dM_MgSiO3_er) + M_MgO * (
                                                                                                     dM_FeO_er + dM_FeSiO3_er + dM_MgO_er + dM_MgSiO3_er) + M_m * (
                                                                                                     -dM_FeO_er - dM_MgO_er + dM_SiO2_er)))) + M_Si * (
                                                                                                     M_Fe * (M_MgO * (
                                                                                                     M_FeO * (
                                                                                                     -dM_FeO_er - 3.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                                                     dM_FeO_er + 3.0 * dM_FeSiO3_er + 3.0 * dM_MgO_er + 5.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_m * (
                                                                                                     2.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er)) + M_MgSiO3 * (
                                                                                                             M_FeO * (
                                                                                                             -dM_FeO_er - 3.0 * dM_FeSiO3_er + 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                                                             dM_FeO_er + 3.0 * dM_FeSiO3_er + dM_MgO_er + 3.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_MgO * (
                                                                                                             4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er) + M_m * (
                                                                                                             2.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er + 4.0 * dM_SiO2_er)) + M_O * (
                                                                                                             M_MgO * (
                                                                                                             dM_MgO_er + dM_MgSiO3_er) + M_MgSiO3 * (
                                                                                                             dM_MgO_er + dM_MgSiO3_er))) + M_Mg * (
                                                                                                     M_Fe * (M_FeO * (
                                                                                                     -dM_FeO_er - 3.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                                                             dM_FeO_er + 3.0 * dM_FeSiO3_er + dM_MgO_er + 3.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_MgO * (
                                                                                                             -dM_FeO_er - 3.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_MgSiO3 * (
                                                                                                             dM_FeO_er + 3.0 * dM_FeSiO3_er + dM_MgO_er + 3.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_m * (
                                                                                                             2.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er)) + M_FeO * (
                                                                                                     M_MgO * (
                                                                                                     -dM_FeO_er - 3.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                                                                                                     4.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                                                     M_FeO * (
                                                                                                     4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er) + M_MgO * (
                                                                                                     1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                                                                                                     -4.0 * dM_FeO_er + 2.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er)) + M_MgSiO3 * (
                                                                                                     M_FeO * (
                                                                                                     3.0 * dM_FeO_er + 5.0 * dM_FeSiO3_er + dM_MgO_er + 3.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                                                     dM_FeO_er + 3.0 * dM_FeSiO3_er + dM_MgO_er + 3.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er)) + M_O * (
                                                                                                     M_FeO * (
                                                                                                     dM_FeO_er + dM_FeSiO3_er) + M_FeSiO3 * (
                                                                                                     dM_FeO_er + dM_FeSiO3_er))) + M_MgO * (
                                                                                                     M_FeO * M_m * (
                                                                                                     4.0 * dM_FeSiO3_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                                                     M_FeO * (
                                                                                                     4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er) + M_m * (
                                                                                                     -4.0 * dM_FeO_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er))) + M_MgSiO3 * (
                                                                                                     M_FeO * (M_MgO * (
                                                                                                     4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er) + M_m * (
                                                                                                              4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er + 4.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                                                     M_FeO * (
                                                                                                     4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                                     4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er) + M_m * (
                                                                                                     -4.0 * dM_FeO_er - 4.0 * dM_MgO_er + 4.0 * dM_SiO2_er))) + M_O * (
                                                                                                     M_MgO * (M_FeO * (
                                                                                                     dM_FeO_er + dM_FeSiO3_er + dM_MgO_er + dM_MgSiO3_er) + M_FeSiO3 * (
                                                                                                              dM_FeO_er + dM_FeSiO3_er + dM_MgO_er + dM_MgSiO3_er)) + M_MgSiO3 * (
                                                                                                     M_FeO * (
                                                                                                     dM_FeO_er + dM_FeSiO3_er + dM_MgO_er + dM_MgSiO3_er) + M_FeSiO3 * (
                                                                                                     dM_FeO_er + dM_FeSiO3_er + dM_MgO_er + dM_MgSiO3_er)))))) / (
               M_O * (M_Fe * (M_MgO * (4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                              M_FeO * (M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                              M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                              -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m))) + M_Mg * (M_Fe * (M_FeSiO3 * (
               M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
               -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO)) + M_MgO * (
                                                                                           4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (M_MgO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)))) + M_Si * (
               M_Fe * (
               M_MgO * (M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (M_MgO * (
               M_FeO * (-M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
               -9.0 * M_FeO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               -9.0 * M_MgO - M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                                                                                             -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                                                                                             -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m))) + M_Mg * (
               M_Fe * (M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                       M_FeO * (-M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                       -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                       -M_SiO2 + 4.0 * M_m) + M_MgSiO3 * (
                       -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                       M_FeO + M_MgO)) + M_MgO * (
               M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (
               M_FeO * (M_MgO * (-M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-9.0 * M_MgO - 9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (
               -M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (-9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_FeSiO3 * (
               -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m)))) + M_O * (M_MgO * (
               9.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)) + M_MgSiO3 * (M_FeO * (
               M_MgO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (M_FeO * (
               -9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)))) + M_c * (
               M_Fe * (
               M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
               M_MgO * (M_FeSiO3 * (M_FeO - M_m) + M_SiO2 * (M_FeO - M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO + M_SiO2) + M_FeSiO3 * (M_FeO + M_MgO - M_m) + M_MgO * (
               M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Mg * (M_Fe * (
               M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
               M_FeSiO3 * (M_FeO + M_MgO - M_m) + M_MgSiO3 * (M_FeO + M_MgO - M_m) + M_SiO2 * (
               M_FeO + M_MgO - M_m)) + M_SiO2 * M_m * (-M_FeO - M_MgO)) + M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                         M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                         M_FeO * (M_SiO2 - M_m) + M_MgO * (
                                                         M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
                                                         M_FeO * M_SiO2 * (M_MgO - M_m) + M_FeSiO3 * (
                                                         M_FeO * (M_MgO + M_SiO2 - M_m) + M_SiO2 * (
                                                         M_MgO - M_m)) + M_MgSiO3 * (
                                                         M_FeO * (M_MgO - M_m) + M_FeSiO3 * (
                                                         M_FeO + M_MgO - M_m)))) + M_O * (
               M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Si * (M_Fe * (M_MgO * (
               M_FeO * (M_SiO2 - M_m) + M_FeSiO3 * (
               4.0 * M_FeO + M_SiO2 - 9.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               4.0 * M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                     4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_O * (
                                                                                                   M_MgO * (
                                                                                                   M_FeO + M_FeSiO3 + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                   M_FeO + M_FeSiO3 + M_SiO2 - M_m))) + M_Mg * (
                                                                                           M_Fe * (M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                                   M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_O * (
                                                                                                   M_FeO + M_FeSiO3 + M_MgO + M_MgSiO3 + M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeO * (
                                                                                           M_MgO * (
                                                                                           M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_MgO + 4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                           M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_FeSiO3 * (
                                                                                           4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m)) + M_O * (
                                                                                           M_FeO * (
                                                                                           M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                           M_MgO + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                           M_FeO + M_FeSiO3))) + M_MgO * (
                                                                                           -4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (M_MgO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_O * (
                                                                                           M_MgO * (M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                    M_SiO2 - M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                           M_SiO2 - M_m))))))

    def dM_FeSiO3_dTc(self, Moles, dKs, dMi_b):
        dM_Mg_er, dM_Si_er, dM_Fe_er, dM_O_er, dM_c_er, dM_MgO_er, dM_SiO2_er, dM_FeO_er, dM_MgSiO3_er, dM_FeSiO3_er, dM_m_er = self.unwrap_Moles(
            dMi_b)
        M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = Moles
        dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs
        return M_FeSiO3 * (M_Fe * (M_O * (M_MgO * (M_FeO * (M_SiO2 * (
        4.0 * dM_FeO_er + 8.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 8.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
                                                            -4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                   -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er)) + M_MgSiO3 * (M_FeO * (
        M_SiO2 * (4.0 * dM_FeO_er + 8.0 * dM_FeSiO3_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
        -4.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er - 4.0 * dM_SiO2_er)) + M_MgO * (M_SiO2 * (
        4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er) + M_m * (-4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                                         -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er))) + dKFeO_KFeO * (
                                   M_Mg * M_O * (-4.0 * M_MgO * M_SiO2 * M_m + M_MgSiO3 * (
                                   M_FeO * (4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                   4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_Si * (M_Mg * (
                                   -M_MgO * M_SiO2 * M_m + M_MgSiO3 * (
                                   M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
                                   M_FeO * (3.0 * M_SiO2 + 6.0 * M_m) + M_MgO * (M_SiO2 - 4.0 * M_m) + M_MgSiO3 * (
                                   9.0 * M_FeO + 9.0 * M_MgO + 4.0 * M_SiO2 - 25.0 * M_m) - 9.0 * M_SiO2 * M_m)) + M_O * (
                                                                                              M_MgO * (M_FeO * (
                                                                                              3.0 * M_SiO2 + 6.0 * M_m) - 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                              M_FeO * (
                                                                                              3.0 * M_SiO2 + 6.0 * M_m) + M_MgO * (
                                                                                              9.0 * M_SiO2 - 9.0 * M_m) - 9.0 * M_SiO2 * M_m))) + M_c * (
                                   M_Mg * (M_MgO * M_SiO2 * M_m + M_MgSiO3 * (
                                   M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                                           M_MgSiO3 * (-M_FeO - M_MgO + M_m) + M_SiO2 * (
                                           -M_FeO - M_MgO + M_m))) + M_O * (
                                   M_MgO * M_SiO2 * (-M_FeO + M_m) + M_MgSiO3 * (
                                   M_MgO * (-M_SiO2 + M_m) + M_SiO2 * (-M_FeO + M_m))) + M_Si * (M_Mg * (
                                   M_FeO * (-2.0 * M_SiO2 - 2.0 * M_m) + M_MgO * (-M_SiO2 + M_m) + M_MgSiO3 * (
                                   -4.0 * M_FeO - 4.0 * M_MgO - M_SiO2 + 9.0 * M_m) + M_O * (
                                   -M_MgO - M_MgSiO3 - M_SiO2 + M_m) + 4.0 * M_SiO2 * M_m) + M_MgO * (M_FeO * (
                                   -2.0 * M_SiO2 - 2.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
                                   -2.0 * M_SiO2 - 2.0 * M_m) + M_MgO * (
                                                                                                  -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_O * (
                                                                                                 M_MgO * (
                                                                                                 -M_SiO2 + M_m) + M_MgSiO3 * (
                                                                                                 -M_SiO2 + M_m)))))) + M_Mg * (
                           M_O * (M_Fe * (M_FeO * (M_SiO2 * (
                           4.0 * dM_FeO_er + 8.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 8.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
                                                   -4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)) + M_MgO * (
                                          M_SiO2 * (
                                          4.0 * dM_FeO_er + 8.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 8.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
                                          -4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                          -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er)) + M_MgO * (
                                  M_FeO * (M_SiO2 * (
                                  4.0 * dM_FeO_er + 8.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 8.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
                                           -4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                  -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er)) + M_MgSiO3 * (M_FeO * (M_SiO2 * (
                           4.0 * dM_FeO_er + 8.0 * dM_FeSiO3_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
                                                                                                 -4.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er - 4.0 * dM_SiO2_er)) + M_MgO * (
                                                                                        M_SiO2 * (
                                                                                        4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er) + M_m * (
                                                                                        -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                        -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er))) + dKMgO_KMgO * (
                           M_Fe * M_O * (4.0 * M_MgO * M_SiO2 * M_m + M_MgSiO3 * (
                           M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_Si * (M_Fe * (
                           M_MgO * M_SiO2 * M_m + M_MgSiO3 * (
                           M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                           M_MgO * (2.0 * M_SiO2 + 10.0 * M_m) + M_MgSiO3 * (
                           -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m))) + M_FeO * M_O * (M_MgO * (
                           3.0 * M_SiO2 + 6.0 * M_m) + M_MgSiO3 * (-6.0 * M_SiO2 + 15.0 * M_m))) + M_c * (M_Fe * (
                           -M_MgO * M_SiO2 * M_m + M_MgSiO3 * (M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) + M_O * (
                           M_FeO + M_MgO - M_m) - M_SiO2 * M_m)) + M_FeO * M_O * (
                                                                                                          -M_MgO * M_SiO2 - M_MgSiO3 * M_m) + M_Si * (
                                                                                                          M_Fe * (
                                                                                                          M_MgO * (
                                                                                                          -M_SiO2 - 3.0 * M_m) + M_MgSiO3 * (
                                                                                                          4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_O * (
                                                                                                          M_MgO + M_MgSiO3)) + M_FeO * (
                                                                                                          M_MgO * (
                                                                                                          -2.0 * M_SiO2 - 2.0 * M_m) + M_MgSiO3 * (
                                                                                                          2.0 * M_SiO2 - 6.0 * M_m)))))) + M_MgSiO3 * dKMgSiO3_KMgSiO3 * (
                           M_O * (M_Fe * M_FeO * M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_Mg * (M_Fe * (
                           M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeO * M_MgO * (
                                                                                               -4.0 * M_SiO2 + 4.0 * M_m))) + M_Si * (
                           M_Mg * (M_Fe * (M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_O * (
                           -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_SiO2 * M_m) + M_FeO * (
                                   M_MgO * (-M_SiO2 + M_m) + M_O * (
                                   -9.0 * M_MgO - 6.0 * M_SiO2 + 15.0 * M_m))) + M_MgO * (M_Fe * (
                           M_FeO * (-M_SiO2 + M_m) + M_O * (-9.0 * M_FeO - 6.0 * M_SiO2 + 15.0 * M_m)) + M_FeO * M_O * (
                                                                                          -9.0 * M_SiO2 + 9.0 * M_m))) + M_c * (
                           M_Mg * (M_Fe * (M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) + M_O * (
                           M_FeO + M_MgO - M_m) - M_SiO2 * M_m) + M_FeO * (
                                   M_MgO * (M_SiO2 - M_m) + M_O * (M_MgO - M_m))) + M_MgO * (
                           M_Fe * (M_FeO * (M_SiO2 - M_m) + M_O * (M_FeO - M_m)) + M_FeO * M_O * (
                           M_SiO2 - M_m)) + M_Si * (M_Mg * (
                           M_Fe * (4.0 * M_FeO + 4.0 * M_MgO + M_O + M_SiO2 - 9.0 * M_m) + M_FeO * (
                           4.0 * M_MgO + 2.0 * M_SiO2 - 6.0 * M_m)) + M_MgO * (
                                                    M_Fe * (4.0 * M_FeO + 2.0 * M_SiO2 - 6.0 * M_m) + M_FeO * (
                                                    4.0 * M_SiO2 - 4.0 * M_m))))) + M_Si * (M_Fe * (M_MgO * (M_FeO * (
        M_SiO2 * (dM_FeO_er + 2 * dM_FeSiO3_er + dM_MgO_er + 2 * dM_MgSiO3_er + dM_SiO2_er) + M_m * (
        -1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                                             -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er)) + M_MgSiO3 * (
                                                                                                    M_FeO * (M_SiO2 * (
                                                                                                    dM_FeO_er + 2 * dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                                                                                                             -1.0 * dM_FeSiO3_er + 1.0 * dM_MgO_er - 1.0 * dM_SiO2_er)) + M_MgO * (
                                                                                                    M_SiO2 * (
                                                                                                    1.0 * dM_FeO_er + 1.0 * dM_FeSiO3_er) + M_m * (
                                                                                                    -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                                    -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er)) + M_O * (
                                                                                                    M_MgO * (M_FeO * (
                                                                                                    6.0 * dM_FeO_er + 15.0 * dM_FeSiO3_er + 6.0 * dM_MgO_er + 15.0 * dM_MgSiO3_er + 9.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                                             4.0 * dM_FeO_er + 10.0 * dM_FeSiO3_er + 6.0 * dM_MgO_er + 12.0 * dM_MgSiO3_er + 6.0 * dM_SiO2_er) + M_m * (
                                                                                                             -10.0 * dM_FeO_er - 25.0 * dM_FeSiO3_er - 15.0 * dM_MgSiO3_er - 15.0 * dM_SiO2_er)) + M_MgSiO3 * (
                                                                                                    M_FeO * (
                                                                                                    6.0 * dM_FeO_er + 15.0 * dM_FeSiO3_er - 3.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 9.0 * dM_SiO2_er) + M_MgO * (
                                                                                                    -9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                    4.0 * dM_FeO_er + 10.0 * dM_FeSiO3_er + 6.0 * dM_MgSiO3_er + 6.0 * dM_SiO2_er) + M_m * (
                                                                                                    -10.0 * dM_FeO_er - 25.0 * dM_FeSiO3_er + 15.0 * dM_MgO_er - 15.0 * dM_SiO2_er)))) + M_Mg * (
                                                                                            M_Fe * (M_FeO * (M_SiO2 * (
                                                                                            dM_FeO_er + 2 * dM_FeSiO3_er + dM_MgO_er + 2 * dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                                                                                                             -1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er)) + M_MgO * (
                                                                                                    M_SiO2 * (
                                                                                                    dM_FeO_er + 2 * dM_FeSiO3_er + dM_MgO_er + 2 * dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                                                                                                    -1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er)) + M_O * (
                                                                                                    M_FeO * (
                                                                                                    6.0 * dM_FeO_er + 15.0 * dM_FeSiO3_er + 6.0 * dM_MgO_er + 15.0 * dM_MgSiO3_er + 9.0 * dM_SiO2_er) + M_MgO * (
                                                                                                    6.0 * dM_FeO_er + 15.0 * dM_FeSiO3_er + 6.0 * dM_MgO_er + 15.0 * dM_MgSiO3_er + 9.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                                    4.0 * dM_FeO_er + 10.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 10.0 * dM_MgSiO3_er + 6.0 * dM_SiO2_er) + M_m * (
                                                                                                    -10.0 * dM_FeO_er - 25.0 * dM_FeSiO3_er - 10.0 * dM_MgO_er - 25.0 * dM_MgSiO3_er - 15.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                                    -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er)) + M_MgO * (
                                                                                            M_FeO * (M_SiO2 * (
                                                                                            dM_FeO_er + 2 * dM_FeSiO3_er + dM_MgO_er + 2 * dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                                                                                                     -1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                            -dM_FeO_er - dM_FeSiO3_er)) + M_MgSiO3 * (
                                                                                            M_FeO * (M_SiO2 * (
                                                                                            dM_FeO_er + 2.0 * dM_FeSiO3_er + 1.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er) + M_m * (
                                                                                                     -1.0 * dM_FeSiO3_er + 1.0 * dM_MgO_er - 1.0 * dM_SiO2_er)) + M_MgO * (
                                                                                            M_SiO2 * (
                                                                                            dM_FeO_er + dM_FeSiO3_er) + M_m * (
                                                                                            -dM_FeO_er - dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                            -dM_FeO_er - dM_FeSiO3_er)) + M_O * (
                                                                                            M_FeO * (M_MgO * (
                                                                                            6.0 * dM_FeO_er + 15.0 * dM_FeSiO3_er + 6.0 * dM_MgO_er + 15.0 * dM_MgSiO3_er + 9.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                                     9.0 * dM_FeO_er + 18.0 * dM_FeSiO3_er + 6.0 * dM_MgO_er + 15.0 * dM_MgSiO3_er + 9.0 * dM_SiO2_er) + M_m * (
                                                                                                     -9.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 15.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er)) + M_MgO * (
                                                                                            M_SiO2 * (
                                                                                            dM_FeO_er + dM_FeSiO3_er) + M_m * (
                                                                                            -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er)) + M_MgSiO3 * (
                                                                                            M_FeO * (
                                                                                            15.0 * dM_FeO_er + 24.0 * dM_FeSiO3_er + 6.0 * dM_MgO_er + 15.0 * dM_MgSiO3_er + 9.0 * dM_SiO2_er) + M_MgO * (
                                                                                            9.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er) + M_SiO2 * (
                                                                                            4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er) + M_m * (
                                                                                            -25.0 * dM_FeO_er - 25.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                            -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er))) + M_O * (
                                                                                            M_MgO * (M_FeO * (M_SiO2 * (
                                                                                            9.0 * dM_FeO_er + 18.0 * dM_FeSiO3_er + 9.0 * dM_MgO_er + 18.0 * dM_MgSiO3_er + 9.0 * dM_SiO2_er) + M_m * (
                                                                                                              -9.0 * dM_FeSiO3_er - 9.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                                     -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er)) + M_MgSiO3 * (
                                                                                            M_FeO * (M_SiO2 * (
                                                                                            9.0 * dM_FeO_er + 18.0 * dM_FeSiO3_er + 9.0 * dM_MgSiO3_er + 9.0 * dM_SiO2_er) + M_m * (
                                                                                                     -9.0 * dM_FeSiO3_er + 9.0 * dM_MgO_er - 9.0 * dM_SiO2_er)) + M_MgO * (
                                                                                            M_SiO2 * (
                                                                                            9.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er) + M_m * (
                                                                                            -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                            -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er))) + dKSiO2_KSiO2 * (
                                                                                            M_O * (M_Fe * (M_MgO * (
                                                                                            M_FeO * (
                                                                                            -2.0 * M_SiO2 - 4.0 * M_m) + 6.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           -2.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                                           -6.0 * M_SiO2 + 6.0 * M_m) + 6.0 * M_SiO2 * M_m)) + M_Mg * (
                                                                                                   M_Fe * (M_FeO * (
                                                                                                   -2.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                                           -2.0 * M_SiO2 - 4.0 * M_m) + 6.0 * M_SiO2 * M_m) + M_FeO * (
                                                                                                   M_MgO * (
                                                                                                   -2.0 * M_SiO2 - 4.0 * M_m) + M_MgSiO3 * (
                                                                                                   4.0 * M_SiO2 - 10.0 * M_m)))) + M_c * (
                                                                                            M_Fe * (M_MgO * (M_FeO * (
                                                                                            M_SiO2 + M_m) - 2.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                    M_FeO * (
                                                                                                    M_SiO2 + M_m) + M_MgO * (
                                                                                                    2.0 * M_SiO2 - 2.0 * M_m) - 2.0 * M_SiO2 * M_m) + M_O * (
                                                                                                    M_MgO * (
                                                                                                    -M_FeO + M_m) + M_MgSiO3 * (
                                                                                                    -M_FeO + M_m))) + M_FeO * M_O * (
                                                                                            M_MgO * (
                                                                                            -M_SiO2 + M_m) + M_MgSiO3 * (
                                                                                            -M_SiO2 + M_m)) + M_Mg * (
                                                                                            M_Fe * (M_FeO * (
                                                                                            M_SiO2 + M_m) + M_MgO * (
                                                                                                    M_SiO2 + M_m) + M_O * (
                                                                                                    -M_FeO - M_MgO + M_m) - 2.0 * M_SiO2 * M_m) + M_FeO * (
                                                                                            M_MgO * (
                                                                                            M_SiO2 + M_m) + M_MgSiO3 * (
                                                                                            -M_SiO2 + 3.0 * M_m) + M_O * (
                                                                                            -M_MgO - M_MgSiO3 - M_SiO2 + M_m)))))) + M_c * (
                           M_Fe * (M_MgO * (M_FeO * (M_SiO2 * (
                           -dM_FeO_er - 2.0 * dM_FeSiO3_er - dM_MgO_er - 2.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
                                                     dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_SiO2 * M_m * (
                                            dM_FeO_er + 1.0 * dM_FeSiO3_er)) + M_MgSiO3 * (M_FeO * (
                           M_SiO2 * (-dM_FeO_er - 2.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
                           dM_FeSiO3_er - dM_MgO_er + dM_SiO2_er)) + M_MgO * (M_SiO2 * (
                           -dM_FeO_er - 1.0 * dM_FeSiO3_er) + M_m * (dM_FeO_er + 1.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                           dM_FeO_er + 1.0 * dM_FeSiO3_er)) + M_O * (
                                   M_MgO * (M_FeO * (-dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                   dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_MgSiO3 * (
                                   M_FeO * (-dM_FeSiO3_er + dM_MgO_er - dM_SiO2_er) + M_MgO * (
                                   dM_MgO_er + dM_MgSiO3_er) + M_m * (
                                   dM_FeSiO3_er - dM_MgO_er + dM_SiO2_er)))) + M_Mg * (M_Fe * (M_FeO * (M_SiO2 * (
                           -dM_FeO_er - 2.0 * dM_FeSiO3_er - dM_MgO_er - 2.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
                                                                                                        dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_MgO * (
                                                                                               M_SiO2 * (
                                                                                               -dM_FeO_er - 2.0 * dM_FeSiO3_er - dM_MgO_er - 2.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
                                                                                               dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_O * (
                                                                                               M_FeO * (
                                                                                               -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_MgO * (
                                                                                               -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                                               dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                               dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er)) + M_MgO * (
                                                                                       M_FeO * (M_SiO2 * (
                                                                                       -dM_FeO_er - 2.0 * dM_FeSiO3_er - dM_MgO_er - 2.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
                                                                                                dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                       dM_FeO_er + dM_FeSiO3_er)) + M_MgSiO3 * (
                                                                                       M_FeO * (M_SiO2 * (
                                                                                       -dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                                                1.0 * dM_FeSiO3_er - dM_MgO_er + 1.0 * dM_SiO2_er)) + M_MgO * (
                                                                                       M_SiO2 * (
                                                                                       -dM_FeO_er - dM_FeSiO3_er) + M_m * (
                                                                                       dM_FeO_er + dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                       dM_FeO_er + dM_FeSiO3_er)) + M_O * (
                                                                                       M_FeO * (M_MgO * (
                                                                                       -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_SiO2 * (
                                                                                                -dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                                                dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_MgSiO3 * (
                                                                                       M_FeO * (
                                                                                       -dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_MgO * (
                                                                                       -dM_FeO_er - dM_FeSiO3_er) + M_m * (
                                                                                       dM_FeO_er + dM_FeSiO3_er)) + M_SiO2 * (
                                                                                       M_MgO * (
                                                                                       -dM_FeO_er - dM_FeSiO3_er) + M_m * (
                                                                                       dM_FeO_er + dM_FeSiO3_er)))) + M_O * (
                           M_MgO * (M_FeO * (M_SiO2 * (
                           -dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                             dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_SiO2 * M_m * (
                                    dM_FeO_er + dM_FeSiO3_er)) + M_MgSiO3 * (M_FeO * (
                           M_SiO2 * (-dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                           dM_FeSiO3_er - dM_MgO_er + dM_SiO2_er)) + M_MgO * (M_SiO2 * (
                           -dM_FeO_er - dM_FeSiO3_er) + M_m * (dM_FeO_er + dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                             dM_FeO_er + dM_FeSiO3_er))) + M_Si * (
                           M_Fe * (M_MgO * (M_FeO * (
                           -2.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_SiO2 * (
                                            -dM_FeO_er - 3.0 * dM_FeSiO3_er - 2 * dM_MgO_er - 4.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                                            3.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er + 6.0 * dM_MgSiO3_er + 6.0 * dM_SiO2_er)) + M_MgSiO3 * (
                                   M_FeO * (
                                   -2.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_MgO * (
                                   4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er) + M_SiO2 * (
                                   -dM_FeO_er - 3.0 * dM_FeSiO3_er - 2.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                                   3.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er + 6.0 * dM_SiO2_er)) + M_O * (
                                   M_MgO * (dM_MgO_er + dM_MgSiO3_er) + M_MgSiO3 * (
                                   dM_MgO_er + dM_MgSiO3_er))) + M_Mg * (M_Fe * (M_FeO * (
                           -2.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_MgO * (
                                                                                 -2.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                 -dM_FeO_er - 3.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                                                                                 3.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er + 3.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er + 6.0 * dM_SiO2_er)) + M_FeO * (
                                                                         M_MgO * (
                                                                         -2.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_SiO2 * (
                                                                         -4.0 * dM_FeO_er - 8.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_m * (
                                                                         4.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er)) + M_MgO * (
                                                                         M_SiO2 * (-dM_FeO_er - dM_FeSiO3_er) + M_m * (
                                                                         dM_FeO_er + dM_FeSiO3_er)) + M_MgSiO3 * (
                                                                         M_FeO * (
                                                                         -6.0 * dM_FeO_er - 10.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_MgO * (
                                                                         -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er) + M_SiO2 * (
                                                                         -dM_FeO_er - dM_FeSiO3_er) + M_m * (
                                                                         9.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er)) + M_O * (
                                                                         M_MgO * (
                                                                         -dM_FeO_er - dM_FeSiO3_er) + M_MgSiO3 * (
                                                                         -dM_FeO_er - dM_FeSiO3_er) + M_SiO2 * (
                                                                         -dM_FeO_er - dM_FeSiO3_er) + M_m * (
                                                                         dM_FeO_er + dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                         4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er)) + M_MgO * (
                           M_FeO * (M_SiO2 * (
                           -4.0 * dM_FeO_er - 8.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 8.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_m * (
                                    4.0 * dM_FeSiO3_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                           4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er)) + M_MgSiO3 * (M_FeO * (M_SiO2 * (
                           -4.0 * dM_FeO_er - 8.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_m * (
                                                                                         4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er + 4.0 * dM_SiO2_er)) + M_MgO * (
                                                                                M_SiO2 * (
                                                                                -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er) + M_m * (
                                                                                4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er)) + M_O * (
                           M_MgO * (
                           M_SiO2 * (-dM_FeO_er - dM_FeSiO3_er) + M_m * (dM_FeO_er + dM_FeSiO3_er)) + M_MgSiO3 * (
                           M_SiO2 * (-dM_FeO_er - dM_FeSiO3_er) + M_m * (
                           dM_FeO_er + dM_FeSiO3_er))))) + dKFeSiO3_KFeSiO3 * (M_O * (M_Fe * M_FeO * (
        -4.0 * M_MgO * M_SiO2 * M_m + M_MgSiO3 * (M_MgO * (4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_Mg * (
                                                                                      M_Fe * (M_MgSiO3 * (M_FeO * (
                                                                                      4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                                          4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                                                                                              -4.0 * M_FeO - 4.0 * M_MgO)) + M_FeO * (
                                                                                      -4.0 * M_MgO * M_SiO2 * M_m + M_MgSiO3 * (
                                                                                      M_MgO * (
                                                                                      4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)))) + M_Si * (
                                                                               M_Fe * (M_FeO * (
                                                                               -M_MgO * M_SiO2 * M_m + M_MgSiO3 * (
                                                                               M_MgO * (
                                                                               M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
                                                                                       M_MgO * (M_FeO * (
                                                                                       M_SiO2 - 4.0 * M_m) - 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                       M_FeO * (
                                                                                       9.0 * M_MgO + M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                       9.0 * M_SiO2 - 9.0 * M_m) - 9.0 * M_SiO2 * M_m))) + M_FeO * M_O * (
                                                                               -9.0 * M_MgO * M_SiO2 * M_m + M_MgSiO3 * (
                                                                               M_MgO * (
                                                                               9.0 * M_SiO2 - 9.0 * M_m) - 9.0 * M_SiO2 * M_m)) + M_Mg * (
                                                                               M_Fe * (M_MgSiO3 * (
                                                                               M_FeO * (M_SiO2 - M_m) + M_MgO * (
                                                                               M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
                                                                                       M_FeO * (
                                                                                       M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                       M_SiO2 - 4.0 * M_m) + M_MgSiO3 * (
                                                                                       9.0 * M_FeO + 9.0 * M_MgO + 4.0 * M_SiO2 - 25.0 * M_m) - 9.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                                                                                       -M_FeO - M_MgO)) + M_FeO * (
                                                                               -M_MgO * M_SiO2 * M_m + M_MgSiO3 * (
                                                                               M_MgO * (
                                                                               M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
                                                                               M_MgO * (
                                                                               M_SiO2 - 4.0 * M_m) + M_MgSiO3 * (
                                                                               9.0 * M_MgO + 4.0 * M_SiO2 - 25.0 * M_m) - 9.0 * M_SiO2 * M_m)))) + M_c * (
                                                                               M_Fe * (M_FeO * (
                                                                               M_MgO * M_SiO2 * M_m + M_MgSiO3 * (
                                                                               M_MgO * (
                                                                               -M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (
                                                                                       M_MgO * M_SiO2 * (
                                                                                       -M_FeO + M_m) + M_MgSiO3 * (
                                                                                       M_FeO * (
                                                                                       -M_MgO - M_SiO2) + M_MgO * (
                                                                                       -M_SiO2 + M_m) + M_SiO2 * M_m))) + M_FeO * M_O * (
                                                                               M_MgO * M_SiO2 * M_m + M_MgSiO3 * (
                                                                               M_MgO * (
                                                                               -M_SiO2 + M_m) + M_SiO2 * M_m)) + M_Mg * (
                                                                               M_Fe * (M_MgSiO3 * (
                                                                               M_FeO * (-M_SiO2 + M_m) + M_MgO * (
                                                                               -M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                                                                                       M_MgSiO3 * (
                                                                                       -M_FeO - M_MgO + M_m) + M_SiO2 * (
                                                                                       -M_FeO - M_MgO + M_m)) + M_SiO2 * M_m * (
                                                                                       M_FeO + M_MgO)) + M_FeO * (
                                                                               M_MgO * M_SiO2 * M_m + M_MgSiO3 * (
                                                                               M_MgO * (
                                                                               -M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                                                                               M_MgSiO3 * (-M_MgO + M_m) + M_SiO2 * (
                                                                               -M_MgO + M_m)))) + M_Si * (M_Fe * (
                                                                               M_MgO * (M_FeO * (
                                                                               -M_SiO2 + M_m) + 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                               M_FeO * (
                                                                               -4.0 * M_MgO - M_SiO2 + M_m) + M_MgO * (
                                                                               -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_O * (
                                                                               M_MgO * (
                                                                               -M_FeO - M_SiO2 + M_m) + M_MgSiO3 * (
                                                                               -M_FeO - M_SiO2 + M_m))) + M_FeO * (
                                                                                                          4.0 * M_MgO * M_SiO2 * M_m + M_MgSiO3 * (
                                                                                                          M_MgO * (
                                                                                                          -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_O * (
                                                                                                          M_MgO * (
                                                                                                          -M_SiO2 + M_m) + M_MgSiO3 * (
                                                                                                          -M_SiO2 + M_m))) + M_Mg * (
                                                                                                          M_Fe * (
                                                                                                          M_FeO * (
                                                                                                          -M_SiO2 + M_m) + M_MgO * (
                                                                                                          -M_SiO2 + M_m) + M_MgSiO3 * (
                                                                                                          -4.0 * M_FeO - 4.0 * M_MgO - M_SiO2 + 9.0 * M_m) + M_O * (
                                                                                                          -M_FeO - M_MgO - M_MgSiO3 - M_SiO2 + M_m) + 4.0 * M_SiO2 * M_m) + M_FeO * (
                                                                                                          M_MgO * (
                                                                                                          -M_SiO2 + M_m) + M_MgSiO3 * (
                                                                                                          -4.0 * M_MgO - M_SiO2 + 9.0 * M_m) + M_O * (
                                                                                                          -M_MgO - M_MgSiO3 - M_SiO2 + M_m) + 4.0 * M_SiO2 * M_m)))))) / (
               M_O * (M_Fe * (M_MgO * (4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                              M_FeO * (M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                              M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                              -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m))) + M_Mg * (M_Fe * (M_FeSiO3 * (
               M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
               -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO)) + M_MgO * (
                                                                                           4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (M_MgO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)))) + M_Si * (
               M_Fe * (
               M_MgO * (M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (M_MgO * (
               M_FeO * (-M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
               -9.0 * M_FeO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               -9.0 * M_MgO - M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                                                                                             -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                                                                                             -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m))) + M_Mg * (
               M_Fe * (M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                       M_FeO * (-M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                       -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                       -M_SiO2 + 4.0 * M_m) + M_MgSiO3 * (
                       -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                       M_FeO + M_MgO)) + M_MgO * (
               M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (
               M_FeO * (M_MgO * (-M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-9.0 * M_MgO - 9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (
               -M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (-9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_FeSiO3 * (
               -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m)))) + M_O * (M_MgO * (
               9.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)) + M_MgSiO3 * (M_FeO * (
               M_MgO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (M_FeO * (
               -9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)))) + M_c * (
               M_Fe * (
               M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
               M_MgO * (M_FeSiO3 * (M_FeO - M_m) + M_SiO2 * (M_FeO - M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO + M_SiO2) + M_FeSiO3 * (M_FeO + M_MgO - M_m) + M_MgO * (
               M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Mg * (M_Fe * (
               M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
               M_FeSiO3 * (M_FeO + M_MgO - M_m) + M_MgSiO3 * (M_FeO + M_MgO - M_m) + M_SiO2 * (
               M_FeO + M_MgO - M_m)) + M_SiO2 * M_m * (-M_FeO - M_MgO)) + M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                         M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                         M_FeO * (M_SiO2 - M_m) + M_MgO * (
                                                         M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
                                                         M_FeO * M_SiO2 * (M_MgO - M_m) + M_FeSiO3 * (
                                                         M_FeO * (M_MgO + M_SiO2 - M_m) + M_SiO2 * (
                                                         M_MgO - M_m)) + M_MgSiO3 * (
                                                         M_FeO * (M_MgO - M_m) + M_FeSiO3 * (
                                                         M_FeO + M_MgO - M_m)))) + M_O * (
               M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Si * (M_Fe * (M_MgO * (
               M_FeO * (M_SiO2 - M_m) + M_FeSiO3 * (
               4.0 * M_FeO + M_SiO2 - 9.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               4.0 * M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                     4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_O * (
                                                                                                   M_MgO * (
                                                                                                   M_FeO + M_FeSiO3 + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                   M_FeO + M_FeSiO3 + M_SiO2 - M_m))) + M_Mg * (
                                                                                           M_Fe * (M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                                   M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_O * (
                                                                                                   M_FeO + M_FeSiO3 + M_MgO + M_MgSiO3 + M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeO * (
                                                                                           M_MgO * (
                                                                                           M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_MgO + 4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                           M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_FeSiO3 * (
                                                                                           4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m)) + M_O * (
                                                                                           M_FeO * (
                                                                                           M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                           M_MgO + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                           M_FeO + M_FeSiO3))) + M_MgO * (
                                                                                           -4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (M_MgO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_O * (
                                                                                           M_MgO * (M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                    M_SiO2 - M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                           M_SiO2 - M_m))))))

    def dM_FeO_dTc(self, Moles, dKs, dMi_b):
        dM_Mg_er, dM_Si_er, dM_Fe_er, dM_O_er, dM_c_er, dM_MgO_er, dM_SiO2_er, dM_FeO_er, dM_MgSiO3_er, dM_FeSiO3_er, dM_m_er = self.unwrap_Moles(
            dMi_b)
        M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = Moles
        dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs
        return M_FeO * (M_Fe * (M_O * (M_MgO * (M_FeSiO3 * (
        M_SiO2 * (-4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 8.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_m * (
        -4.0 * dM_FeO_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er)) + M_MgSiO3 * (M_FeSiO3 * (
        M_SiO2 * (-4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_m * (
        -4.0 * dM_FeO_er - 4.0 * dM_MgO_er + 4.0 * dM_SiO2_er)) + M_MgO * (M_SiO2 * (
        4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er) + M_m * (-4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                                      -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er))) + dKFeO_KFeO * (
                                M_Mg * M_O * (
                                M_MgO * (M_FeSiO3 * (4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                M_MgO * (4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_Si * (M_Mg * (
                                M_MgO * (M_FeSiO3 * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_MgSiO3 * (
                                M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
                                M_FeSiO3 * (9.0 * M_MgO + 6.0 * M_SiO2 - 15.0 * M_m) + M_MgO * (
                                M_SiO2 - 4.0 * M_m) + M_MgSiO3 * (
                                9.0 * M_MgO + 4.0 * M_SiO2 - 25.0 * M_m) - 9.0 * M_SiO2 * M_m)) + M_O * (M_MgO * (
                                M_FeSiO3 * (6.0 * M_SiO2 - 15.0 * M_m) - 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeSiO3 * (
                                6.0 * M_SiO2 - 15.0 * M_m) + M_MgO * (
                                                                                                           9.0 * M_SiO2 - 9.0 * M_m) - 9.0 * M_SiO2 * M_m))) + M_c * (
                                M_Mg * (M_MgO * (M_FeSiO3 * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_MgSiO3 * (
                                M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                                        M_FeSiO3 * (-M_MgO + M_m) + M_MgSiO3 * (-M_MgO + M_m) + M_SiO2 * (
                                        -M_MgO + M_m))) + M_O * (M_MgO * M_m * (M_FeSiO3 + M_SiO2) + M_MgSiO3 * (
                                M_MgO * (-M_SiO2 + M_m) + M_m * (M_FeSiO3 + M_SiO2))) + M_Si * (M_Mg * (
                                M_FeSiO3 * (-4.0 * M_MgO - 2.0 * M_SiO2 + 6.0 * M_m) + M_MgO * (
                                -M_SiO2 + M_m) + M_MgSiO3 * (-4.0 * M_MgO - M_SiO2 + 9.0 * M_m) + M_O * (
                                -M_MgO - M_MgSiO3 - M_SiO2 + M_m) + 4.0 * M_SiO2 * M_m) + M_MgO * (M_FeSiO3 * (
                                -2.0 * M_SiO2 + 6.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeSiO3 * (
                                -2.0 * M_SiO2 + 6.0 * M_m) + M_MgO * (
                                                                                               -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_O * (
                                                                                                M_MgO * (
                                                                                                -M_SiO2 + M_m) + M_MgSiO3 * (
                                                                                                -M_SiO2 + M_m)))))) + M_FeSiO3 * dKFeSiO3_KFeSiO3 * (
                        M_O * (M_Fe * (4.0 * M_MgO * M_SiO2 * M_m + M_MgSiO3 * (
                        M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_Mg * (
                               M_MgSiO3 * (M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                               4.0 * M_Fe + 4.0 * M_MgO))) + M_Si * (M_Fe * (
                        M_MgO * M_SiO2 * M_m + M_MgSiO3 * (M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                        M_MgO * (2.0 * M_SiO2 + 10.0 * M_m) + M_MgSiO3 * (
                        -9.0 * M_MgO + 2.0 * M_SiO2 + 10.0 * M_m))) + M_Mg * (M_Fe * (
                        M_O * (2.0 * M_SiO2 + 10.0 * M_m) + M_SiO2 * M_m) + M_MgO * M_SiO2 * M_m + M_MgSiO3 * (M_MgO * (
                        -M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (M_MgO * (-M_SiO2 + 4.0 * M_m) + M_MgSiO3 * (
                        -9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m)) + M_O * (
                                                                     9.0 * M_MgO * M_SiO2 * M_m + M_MgSiO3 * (M_MgO * (
                                                                     -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m))) + M_c * (
                        M_Fe * (
                        -M_MgO * M_SiO2 * M_m + M_MgSiO3 * (M_MgO * (M_O + M_SiO2 - M_m) - M_SiO2 * M_m)) + M_Mg * (
                        M_MgSiO3 * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
                        M_MgSiO3 * (M_MgO - M_m) + M_SiO2 * (M_MgO - M_m)) + M_SiO2 * M_m * (-M_Fe - M_MgO)) + M_O * (
                        -M_MgO * M_SiO2 * M_m + M_MgSiO3 * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_Si * (M_Fe * (
                        M_MgO * (-M_SiO2 - 3.0 * M_m) + M_MgSiO3 * (4.0 * M_MgO - M_SiO2 - 3.0 * M_m) + M_O * (
                        M_MgO + M_MgSiO3)) + M_Mg * (M_Fe * (M_O - M_SiO2 - 3.0 * M_m) + M_MgO * (
                        M_SiO2 - M_m) + M_MgSiO3 * (4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_O * (
                                                     M_MgO + M_MgSiO3 + M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) - 4.0 * M_MgO * M_SiO2 * M_m + M_MgSiO3 * (
                                                                                                              M_MgO * (
                                                                                                              4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_O * (
                                                                                                              M_MgO * (
                                                                                                              M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                              M_SiO2 - M_m))))) + M_Mg * (
                        M_O * (M_Fe * (M_FeSiO3 * (
                        M_SiO2 * (-4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_m * (
                        -4.0 * dM_FeO_er - 4.0 * dM_MgO_er + 4.0 * dM_SiO2_er)) + M_MgSiO3 * (
                                       M_SiO2 * (-4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_m * (
                                       -4.0 * dM_FeO_er - 4.0 * dM_MgO_er + 4.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                       -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er)) + M_MgO * (
                               M_FeSiO3 * (M_SiO2 * (
                               -4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 8.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_m * (
                                           -4.0 * dM_FeO_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                               -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er)) + M_MgSiO3 * (M_FeSiO3 * (
                        M_SiO2 * (-4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_m * (
                        -4.0 * dM_FeO_er - 4.0 * dM_MgO_er + 4.0 * dM_SiO2_er)) + M_MgO * (M_SiO2 * (
                        4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er) + M_m * (
                                                                                           -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                     -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er))) + dKMgO_KMgO * (
                        M_Fe * M_O * (
                        M_MgO * (M_FeSiO3 * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                        M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_Si * (M_Fe * (
                        M_MgO * (M_FeSiO3 * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_MgSiO3 * (
                        M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                        M_MgO * (-9.0 * M_FeSiO3 - M_SiO2 + 4.0 * M_m) + M_MgSiO3 * (
                        -9.0 * M_MgO + 2.0 * M_SiO2 + 10.0 * M_m))) + M_FeSiO3 * M_O * (M_MgO * (
                        -3.0 * M_SiO2 - 6.0 * M_m) + M_MgSiO3 * (6.0 * M_SiO2 - 15.0 * M_m))) + M_c * (M_Fe * (M_MgO * (
                        M_FeSiO3 * (M_SiO2 - M_m) + M_O * (M_FeSiO3 + M_MgSiO3 + M_SiO2) - M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                               M_MgO * (
                                                                                                               M_SiO2 - M_m) - M_SiO2 * M_m)) + M_FeSiO3 * M_O * (
                                                                                                       M_MgO * M_SiO2 + M_MgSiO3 * M_m) + M_Si * (
                                                                                                       M_Fe * (M_MgO * (
                                                                                                       4.0 * M_FeSiO3 + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                               4.0 * M_MgO - M_SiO2 - 3.0 * M_m) + M_O * (
                                                                                                               M_MgO + M_MgSiO3)) + M_FeSiO3 * (
                                                                                                       M_MgO * (
                                                                                                       2.0 * M_SiO2 + 2.0 * M_m) + M_MgSiO3 * (
                                                                                                       -2.0 * M_SiO2 + 6.0 * M_m)))))) + M_MgSiO3 * dKMgSiO3_KMgSiO3 * (
                        M_O * (M_Fe * M_FeSiO3 * M_MgO * (4.0 * M_SiO2 - 4.0 * M_m) + M_Mg * (
                        4.0 * M_Fe * M_SiO2 * M_m + M_FeSiO3 * M_MgO * (4.0 * M_SiO2 - 4.0 * M_m))) + M_Si * (M_Mg * (
                        M_Fe * (M_O * (2.0 * M_SiO2 + 10.0 * M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
                        M_MgO * (M_SiO2 - M_m) + M_O * (9.0 * M_MgO + 6.0 * M_SiO2 - 15.0 * M_m))) + M_MgO * (M_Fe * (
                        M_FeSiO3 * (M_SiO2 - M_m) + M_O * (
                        9.0 * M_FeSiO3 + 3.0 * M_SiO2 + 6.0 * M_m)) + M_FeSiO3 * M_O * (
                                                                                                              9.0 * M_SiO2 - 9.0 * M_m))) + M_c * (
                        M_Mg * (
                        -M_Fe * M_SiO2 * M_m + M_FeSiO3 * (M_MgO * (-M_SiO2 + M_m) + M_O * (-M_MgO + M_m))) + M_MgO * (
                        M_Fe * (M_FeSiO3 * (-M_SiO2 + M_m) + M_O * (-M_FeSiO3 - M_SiO2)) + M_FeSiO3 * M_O * (
                        -M_SiO2 + M_m)) + M_Si * (M_Mg * (M_Fe * (M_O - M_SiO2 - 3.0 * M_m) + M_FeSiO3 * (
                        -4.0 * M_MgO - 2.0 * M_SiO2 + 6.0 * M_m)) + M_MgO * (
                                                  M_Fe * (-4.0 * M_FeSiO3 - 2.0 * M_SiO2 - 2.0 * M_m) + M_FeSiO3 * (
                                                  -4.0 * M_SiO2 + 4.0 * M_m))))) + M_Si * (M_Fe * (M_MgO * (M_FeSiO3 * (
        M_SiO2 * (-1.0 * dM_FeSiO3_er - dM_MgO_er - 2.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
        -1.0 * dM_FeO_er + 1.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                                            -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er)) + M_MgSiO3 * (
                                                                                                   M_FeSiO3 * (
                                                                                                   M_SiO2 * (
                                                                                                   -1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
                                                                                                   -1.0 * dM_FeO_er - 1.0 * dM_MgO_er + 1.0 * dM_SiO2_er)) + M_MgO * (
                                                                                                   M_SiO2 * (
                                                                                                   1.0 * dM_FeO_er + 1.0 * dM_FeSiO3_er) + M_m * (
                                                                                                   -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                                   -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er)) + M_O * (
                                                                                                   M_MgO * (M_FeSiO3 * (
                                                                                                   -6.0 * dM_FeO_er - 15.0 * dM_FeSiO3_er - 15.0 * dM_MgO_er - 24.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                                            -2 * dM_FeO_er - 5.0 * dM_FeSiO3_er - 3.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_m * (
                                                                                                            -4.0 * dM_FeO_er - 10.0 * dM_FeSiO3_er - 6.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er)) + M_MgSiO3 * (
                                                                                                   M_FeSiO3 * (
                                                                                                   -6.0 * dM_FeO_er - 15.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 15.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er) + M_MgO * (
                                                                                                   -9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                   -2 * dM_FeO_er - 5.0 * dM_FeSiO3_er - 3.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_m * (
                                                                                                   -4.0 * dM_FeO_er - 10.0 * dM_FeSiO3_er + 6.0 * dM_MgO_er - 6.0 * dM_SiO2_er)))) + M_Mg * (
                                                                                           M_Fe * (M_FeSiO3 * (
                                                                                           M_SiO2 * (
                                                                                           -1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
                                                                                           -1.0 * dM_FeO_er - 1.0 * dM_MgO_er + 1.0 * dM_SiO2_er)) + M_MgSiO3 * (
                                                                                                   M_SiO2 * (
                                                                                                   -1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
                                                                                                   -1.0 * dM_FeO_er - 1.0 * dM_MgO_er + 1.0 * dM_SiO2_er)) + M_O * (
                                                                                                   M_FeSiO3 * (
                                                                                                   -6.0 * dM_FeO_er - 15.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 15.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er) + M_MgSiO3 * (
                                                                                                   -6.0 * dM_FeO_er - 15.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 15.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                                   -2 * dM_FeO_er - 5.0 * dM_FeSiO3_er - 2 * dM_MgO_er - 5.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_m * (
                                                                                                   -4.0 * dM_FeO_er - 10.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 10.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                                   -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er)) + M_MgO * (
                                                                                           M_FeSiO3 * (M_SiO2 * (
                                                                                           -dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                                                       -1.0 * dM_FeO_er + 1.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                           -dM_FeO_er - dM_FeSiO3_er)) + M_MgSiO3 * (
                                                                                           M_FeSiO3 * (M_SiO2 * (
                                                                                           -1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
                                                                                                       -1.0 * dM_FeO_er - 1.0 * dM_MgO_er + 1.0 * dM_SiO2_er)) + M_MgO * (
                                                                                           M_SiO2 * (
                                                                                           dM_FeO_er + dM_FeSiO3_er) + M_m * (
                                                                                           -dM_FeO_er - dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                           -dM_FeO_er - dM_FeSiO3_er)) + M_O * (
                                                                                           M_FeSiO3 * (M_MgO * (
                                                                                           3.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 15.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                                       -9.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 15.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er) + M_m * (
                                                                                                       -9.0 * dM_FeO_er + 6.0 * dM_MgO_er + 15.0 * dM_MgSiO3_er + 9.0 * dM_SiO2_er)) + M_MgO * (
                                                                                           M_SiO2 * (
                                                                                           dM_FeO_er + dM_FeSiO3_er) + M_m * (
                                                                                           -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er)) + M_MgSiO3 * (
                                                                                           M_FeSiO3 * (
                                                                                           -6.0 * dM_FeO_er - 15.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 15.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er) + M_MgO * (
                                                                                           9.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er) + M_SiO2 * (
                                                                                           4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er) + M_m * (
                                                                                           -25.0 * dM_FeO_er - 25.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                           -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er))) + M_O * (
                                                                                           M_MgO * (M_FeSiO3 * (
                                                                                           M_SiO2 * (
                                                                                           -9.0 * dM_FeSiO3_er - 9.0 * dM_MgO_er - 18.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er) + M_m * (
                                                                                           -9.0 * dM_FeO_er + 9.0 * dM_MgSiO3_er + 9.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                                    -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er)) + M_MgSiO3 * (
                                                                                           M_FeSiO3 * (M_SiO2 * (
                                                                                           -9.0 * dM_FeSiO3_er - 9.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er) + M_m * (
                                                                                                       -9.0 * dM_FeO_er - 9.0 * dM_MgO_er + 9.0 * dM_SiO2_er)) + M_MgO * (
                                                                                           M_SiO2 * (
                                                                                           9.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er) + M_m * (
                                                                                           -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                           -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er))) + dKSiO2_KSiO2 * (
                                                                                           M_O * (M_Fe * (M_MgO * (
                                                                                           M_FeSiO3 * (
                                                                                           -4.0 * M_SiO2 + 10.0 * M_m) + 6.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                          M_FeSiO3 * (
                                                                                                          -4.0 * M_SiO2 + 10.0 * M_m) + M_MgO * (
                                                                                                          -6.0 * M_SiO2 + 6.0 * M_m) + 6.0 * M_SiO2 * M_m)) + M_Mg * (
                                                                                                  M_Fe * (M_FeSiO3 * (
                                                                                                  -4.0 * M_SiO2 + 10.0 * M_m) + M_MgSiO3 * (
                                                                                                          -4.0 * M_SiO2 + 10.0 * M_m) + 6.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                  M_MgO * (
                                                                                                  2.0 * M_SiO2 + 4.0 * M_m) + M_MgSiO3 * (
                                                                                                  -4.0 * M_SiO2 + 10.0 * M_m)))) + M_c * (
                                                                                           M_Fe * (M_MgO * (M_FeSiO3 * (
                                                                                           M_SiO2 - 3.0 * M_m) - 2.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                   M_FeSiO3 * (
                                                                                                   M_SiO2 - 3.0 * M_m) + M_MgO * (
                                                                                                   2.0 * M_SiO2 - 2.0 * M_m) - 2.0 * M_SiO2 * M_m) + M_O * (
                                                                                                   M_MgO * (
                                                                                                   M_FeSiO3 + M_SiO2) + M_MgSiO3 * (
                                                                                                   M_FeSiO3 + M_SiO2))) + M_FeSiO3 * M_O * (
                                                                                           M_MgO * (
                                                                                           M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                           M_SiO2 - M_m)) + M_Mg * (
                                                                                           M_Fe * (M_FeSiO3 * (
                                                                                           M_SiO2 - 3.0 * M_m) + M_MgSiO3 * (
                                                                                                   M_SiO2 - 3.0 * M_m) + M_O * (
                                                                                                   M_FeSiO3 + M_MgSiO3 + M_SiO2) - 2.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_MgO * (
                                                                                           -M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                           M_SiO2 - 3.0 * M_m) + M_O * (
                                                                                           M_MgO + M_MgSiO3 + M_SiO2 - M_m)))))) + M_c * (
                        M_Fe * (M_MgO * (M_FeSiO3 * (
                        M_SiO2 * (dM_FeSiO3_er + dM_MgO_er + 2 * dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                        dM_FeO_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                         dM_FeO_er + 1.0 * dM_FeSiO3_er)) + M_MgSiO3 * (M_FeSiO3 * (
                        M_SiO2 * (dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                        dM_FeO_er + dM_MgO_er - 1.0 * dM_SiO2_er)) + M_MgO * (M_SiO2 * (
                        -dM_FeO_er - 1.0 * dM_FeSiO3_er) + M_m * (dM_FeO_er + 1.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                        dM_FeO_er + 1.0 * dM_FeSiO3_er)) + M_O * (
                                M_MgO * (
                                M_FeSiO3 * (dM_FeSiO3_er + dM_MgO_er + 2 * dM_MgSiO3_er + dM_SiO2_er) + M_SiO2 * (
                                dM_FeSiO3_er + dM_MgO_er + 2 * dM_MgSiO3_er + dM_SiO2_er)) + M_MgSiO3 * (
                                M_FeSiO3 * (dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_MgO * (
                                dM_MgO_er + dM_MgSiO3_er) + M_SiO2 * (
                                dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)))) + M_Mg * (M_Fe * (M_FeSiO3 * (
                        M_SiO2 * (dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                        dM_FeO_er + dM_MgO_er - 1.0 * dM_SiO2_er)) + M_MgSiO3 * (M_SiO2 * (
                        dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                                                                                 dM_FeO_er + dM_MgO_er - 1.0 * dM_SiO2_er)) + M_O * (
                                                                                               M_FeSiO3 * (
                                                                                               dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_MgSiO3 * (
                                                                                               dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_SiO2 * (
                                                                                               dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                               dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er)) + M_MgO * (
                                                                                       M_FeSiO3 * (M_SiO2 * (
                                                                                       1.0 * dM_FeSiO3_er + dM_MgO_er + 2.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er) + M_m * (
                                                                                                   dM_FeO_er - dM_MgSiO3_er - dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                       dM_FeO_er + dM_FeSiO3_er)) + M_MgSiO3 * (
                                                                                       M_FeSiO3 * (M_SiO2 * (
                                                                                       dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                                                                                                   dM_FeO_er + dM_MgO_er - 1.0 * dM_SiO2_er)) + M_MgO * (
                                                                                       M_SiO2 * (
                                                                                       -dM_FeO_er - dM_FeSiO3_er) + M_m * (
                                                                                       dM_FeO_er + dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                       dM_FeO_er + dM_FeSiO3_er)) + M_O * (
                                                                                       M_FeSiO3 * (M_MgO * (
                                                                                       -dM_FeO_er + dM_MgSiO3_er + dM_SiO2_er) + M_SiO2 * (
                                                                                                   dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                                                                                                   dM_FeO_er - dM_MgSiO3_er - dM_SiO2_er)) + M_MgSiO3 * (
                                                                                       M_FeSiO3 * (
                                                                                       dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_MgO * (
                                                                                       -dM_FeO_er - dM_FeSiO3_er) + M_m * (
                                                                                       dM_FeO_er + dM_FeSiO3_er)) + M_SiO2 * (
                                                                                       M_MgO * (
                                                                                       -dM_FeO_er - dM_FeSiO3_er) + M_m * (
                                                                                       dM_FeO_er + dM_FeSiO3_er)))) + M_O * (
                        M_MgO * (M_FeSiO3 * (
                        M_SiO2 * (dM_FeSiO3_er + dM_MgO_er + 2 * dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                        dM_FeO_er - dM_MgSiO3_er - dM_SiO2_er)) + M_SiO2 * M_m * (
                                 dM_FeO_er + dM_FeSiO3_er)) + M_MgSiO3 * (M_FeSiO3 * (
                        M_SiO2 * (dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                        dM_FeO_er + dM_MgO_er - dM_SiO2_er)) + M_MgO * (M_SiO2 * (-dM_FeO_er - dM_FeSiO3_er) + M_m * (
                        dM_FeO_er + dM_FeSiO3_er)) + M_SiO2 * M_m * (dM_FeO_er + dM_FeSiO3_er))) + M_Si * (M_Fe * (
                        M_MgO * (M_FeSiO3 * (
                        2.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er + 6.0 * dM_MgO_er + 10.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_SiO2 * (
                                 dM_FeO_er + 3.0 * dM_FeSiO3_er + 2 * dM_MgO_er + 4.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_m * (
                                 dM_FeO_er + 3.0 * dM_FeSiO3_er + 2.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er)) + M_MgSiO3 * (
                        M_FeSiO3 * (
                        2.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_MgO * (
                        4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er) + M_SiO2 * (
                        dM_FeO_er + 3.0 * dM_FeSiO3_er + 2.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_m * (
                        dM_FeO_er + 3.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er + 2.0 * dM_SiO2_er)) + M_O * (
                        M_MgO * (dM_MgO_er + dM_MgSiO3_er) + M_MgSiO3 * (dM_MgO_er + dM_MgSiO3_er))) + M_Mg * (M_Fe * (
                        M_FeSiO3 * (
                        2.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_MgSiO3 * (
                        2.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_SiO2 * (
                        dM_FeO_er + 3.0 * dM_FeSiO3_er + dM_MgO_er + 3.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_m * (
                        dM_FeO_er + 3.0 * dM_FeSiO3_er + dM_MgO_er + 3.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                                                               M_MgO * (
                                                                                                               -2.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                                               4.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
                                                                                                               4.0 * dM_FeO_er - 2.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)) + M_MgO * (
                                                                                                               M_SiO2 * (
                                                                                                               -dM_FeO_er - dM_FeSiO3_er) + M_m * (
                                                                                                               dM_FeO_er + dM_FeSiO3_er)) + M_MgSiO3 * (
                                                                                                               M_FeSiO3 * (
                                                                                                               2.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_MgO * (
                                                                                                               -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er) + M_SiO2 * (
                                                                                                               -dM_FeO_er - dM_FeSiO3_er) + M_m * (
                                                                                                               9.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er)) + M_O * (
                                                                                                               M_MgO * (
                                                                                                               -dM_FeO_er - dM_FeSiO3_er) + M_MgSiO3 * (
                                                                                                               -dM_FeO_er - dM_FeSiO3_er) + M_SiO2 * (
                                                                                                               -dM_FeO_er - dM_FeSiO3_er) + M_m * (
                                                                                                               dM_FeO_er + dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                                               4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er)) + M_MgO * (
                                                                                                           M_FeSiO3 * (
                                                                                                           M_SiO2 * (
                                                                                                           4.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 8.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
                                                                                                           4.0 * dM_FeO_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                                           4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er)) + M_MgSiO3 * (
                                                                                                           M_FeSiO3 * (
                                                                                                           M_SiO2 * (
                                                                                                           4.0 * dM_FeSiO3_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
                                                                                                           4.0 * dM_FeO_er + 4.0 * dM_MgO_er - 4.0 * dM_SiO2_er)) + M_MgO * (
                                                                                                           M_SiO2 * (
                                                                                                           -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er) + M_m * (
                                                                                                           4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                                           4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er)) + M_O * (
                                                                                                           M_MgO * (
                                                                                                           M_SiO2 * (
                                                                                                           -dM_FeO_er - dM_FeSiO3_er) + M_m * (
                                                                                                           dM_FeO_er + dM_FeSiO3_er)) + M_MgSiO3 * (
                                                                                                           M_SiO2 * (
                                                                                                           -dM_FeO_er - dM_FeSiO3_er) + M_m * (
                                                                                                           dM_FeO_er + dM_FeSiO3_er)))))) / (
               M_O * (M_Fe * (M_MgO * (4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                              M_FeO * (M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                              M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                              -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m))) + M_Mg * (M_Fe * (M_FeSiO3 * (
               M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
               -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO)) + M_MgO * (
                                                                                           4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (M_MgO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)))) + M_Si * (
               M_Fe * (
               M_MgO * (M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (M_MgO * (
               M_FeO * (-M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
               -9.0 * M_FeO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               -9.0 * M_MgO - M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                                                                                             -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                                                                                             -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m))) + M_Mg * (
               M_Fe * (M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                       M_FeO * (-M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                       -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                       -M_SiO2 + 4.0 * M_m) + M_MgSiO3 * (
                       -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                       M_FeO + M_MgO)) + M_MgO * (
               M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (
               M_FeO * (M_MgO * (-M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-9.0 * M_MgO - 9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (
               -M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (-9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_FeSiO3 * (
               -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m)))) + M_O * (M_MgO * (
               9.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)) + M_MgSiO3 * (M_FeO * (
               M_MgO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (M_FeO * (
               -9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)))) + M_c * (
               M_Fe * (
               M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
               M_MgO * (M_FeSiO3 * (M_FeO - M_m) + M_SiO2 * (M_FeO - M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO + M_SiO2) + M_FeSiO3 * (M_FeO + M_MgO - M_m) + M_MgO * (
               M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Mg * (M_Fe * (
               M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
               M_FeSiO3 * (M_FeO + M_MgO - M_m) + M_MgSiO3 * (M_FeO + M_MgO - M_m) + M_SiO2 * (
               M_FeO + M_MgO - M_m)) + M_SiO2 * M_m * (-M_FeO - M_MgO)) + M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                         M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                         M_FeO * (M_SiO2 - M_m) + M_MgO * (
                                                         M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
                                                         M_FeO * M_SiO2 * (M_MgO - M_m) + M_FeSiO3 * (
                                                         M_FeO * (M_MgO + M_SiO2 - M_m) + M_SiO2 * (
                                                         M_MgO - M_m)) + M_MgSiO3 * (
                                                         M_FeO * (M_MgO - M_m) + M_FeSiO3 * (
                                                         M_FeO + M_MgO - M_m)))) + M_O * (
               M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Si * (M_Fe * (M_MgO * (
               M_FeO * (M_SiO2 - M_m) + M_FeSiO3 * (
               4.0 * M_FeO + M_SiO2 - 9.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               4.0 * M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                     4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_O * (
                                                                                                   M_MgO * (
                                                                                                   M_FeO + M_FeSiO3 + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                   M_FeO + M_FeSiO3 + M_SiO2 - M_m))) + M_Mg * (
                                                                                           M_Fe * (M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                                   M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_O * (
                                                                                                   M_FeO + M_FeSiO3 + M_MgO + M_MgSiO3 + M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeO * (
                                                                                           M_MgO * (
                                                                                           M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_MgO + 4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                           M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_FeSiO3 * (
                                                                                           4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m)) + M_O * (
                                                                                           M_FeO * (
                                                                                           M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                           M_MgO + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                           M_FeO + M_FeSiO3))) + M_MgO * (
                                                                                           -4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (M_MgO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_O * (
                                                                                           M_MgO * (M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                    M_SiO2 - M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                           M_SiO2 - M_m))))))

    def dM_Si_dTc(self, Moles, dKs, dMi_b):
        dM_Mg_er, dM_Si_er, dM_Fe_er, dM_O_er, dM_c_er, dM_MgO_er, dM_SiO2_er, dM_FeO_er, dM_MgSiO3_er, dM_FeSiO3_er, dM_m_er = self.unwrap_Moles(
            dMi_b)
        M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = Moles
        dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs
        return M_Si * (M_Fe * (M_O * (M_MgO * (M_FeO * (
        M_SiO2 * (-2 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 2 * dM_MgO_er - 4.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
        -4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)) + M_FeSiO3 * (M_FeO * (
        -6.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                    4.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 8.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
                                                                                    10.0 * dM_FeO_er - 10.0 * dM_MgSiO3_er - 10.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                               6.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er)) + M_MgSiO3 * (M_FeO * (
        M_MgO * (-6.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er) + M_SiO2 * (
        -2 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 2.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
        -4.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er - 4.0 * dM_SiO2_er)) + M_FeSiO3 * (M_FeO * (
        -6.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                 -6.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                 4.0 * dM_FeSiO3_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
                                                                                 10.0 * dM_FeO_er + 10.0 * dM_MgO_er - 10.0 * dM_SiO2_er)) + M_MgO * (
                                                                                                    M_SiO2 * (
                                                                                                    -6.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er) + M_m * (
                                                                                                    6.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                                    6.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er))) + dKFeO_KFeO * (
                               M_O * (M_Mg * (
                               M_FeO * (M_MgSiO3 * (4.0 * M_SiO2 - 10.0 * M_m) - 6.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                               M_FeO * (6.0 * M_SiO2 - 6.0 * M_m) + M_MgO * (
                               2.0 * M_SiO2 + 4.0 * M_m) - 6.0 * M_SiO2 * M_m)) + M_MgO * (
                                      -6.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                      M_FeO * (6.0 * M_SiO2 - 6.0 * M_m) - 6.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                      M_FeO * (M_MgO * (6.0 * M_SiO2 - 6.0 * M_m) - 6.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                      M_FeO * (6.0 * M_SiO2 - 6.0 * M_m) + M_MgO * (
                                      6.0 * M_SiO2 - 6.0 * M_m) - 6.0 * M_SiO2 * M_m))) + M_c * (M_Mg * (
                               M_FeO * (M_MgSiO3 * (-M_SiO2 + 3.0 * M_m) + 2.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                               M_FeO * (-2.0 * M_SiO2 + 2.0 * M_m) + M_MgO * (
                               -M_SiO2 - M_m) + 2.0 * M_SiO2 * M_m) + M_O * (
                               M_FeO * (-M_MgSiO3 - M_SiO2) + M_FeSiO3 * (M_MgO - M_m))) + M_MgO * (
                                                                                                 2.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                 M_FeO * (
                                                                                                 -2.0 * M_SiO2 + 2.0 * M_m) + 2.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                 M_FeO * (M_MgO * (
                                                                                                 -2.0 * M_SiO2 + 2.0 * M_m) + 2.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                 M_FeO * (
                                                                                                 -2.0 * M_SiO2 + 2.0 * M_m) + M_MgO * (
                                                                                                 -2.0 * M_SiO2 + 2.0 * M_m) + 2.0 * M_SiO2 * M_m)) + M_O * (
                                                                                                 M_MgO * (
                                                                                                 -M_FeO * M_SiO2 - M_FeSiO3 * M_m) + M_MgSiO3 * (
                                                                                                 -M_FeO * M_SiO2 - M_FeSiO3 * M_m))))) + M_FeSiO3 * dKFeSiO3_KFeSiO3 * (
                       M_O * (M_Fe * (M_MgO * (M_FeO * (2.0 * M_SiO2 + 4.0 * M_m) - 6.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                       M_FeO * (2.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                       6.0 * M_SiO2 - 6.0 * M_m) - 6.0 * M_SiO2 * M_m)) + M_Mg * (M_Fe * (
                       M_FeO * (2.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                       2.0 * M_SiO2 + 4.0 * M_m) - 6.0 * M_SiO2 * M_m) + M_FeO * (M_MgO * (
                       2.0 * M_SiO2 + 4.0 * M_m) + M_MgSiO3 * (-4.0 * M_SiO2 + 10.0 * M_m)))) + M_c * (M_Fe * (
                       M_MgO * (M_FeO * (-M_SiO2 - M_m) + 2.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                       M_FeO * (-M_SiO2 - M_m) + M_MgO * (-2.0 * M_SiO2 + 2.0 * M_m) + 2.0 * M_SiO2 * M_m) + M_O * (
                       M_MgO * (M_FeO - M_m) + M_MgSiO3 * (M_FeO - M_m))) + M_FeO * M_O * (M_MgO * (
                       M_SiO2 - M_m) + M_MgSiO3 * (M_SiO2 - M_m)) + M_Mg * (M_Fe * (
                       M_FeO * (-M_SiO2 - M_m) + M_MgO * (-M_SiO2 - M_m) + M_O * (
                       M_FeO + M_MgO - M_m) + 2.0 * M_SiO2 * M_m) + M_FeO * (M_MgO * (-M_SiO2 - M_m) + M_MgSiO3 * (
                       M_SiO2 - 3.0 * M_m) + M_O * (M_MgO + M_MgSiO3 + M_SiO2 - M_m))))) + M_Mg * (M_O * (M_Fe * (
        M_FeO * (
        M_SiO2 * (-2 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 2 * dM_MgO_er - 4.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
        -4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)) + M_FeSiO3 * (
        M_FeO * (-6.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er) + M_MgO * (
        -6.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er) + M_SiO2 * (
        4.0 * dM_FeSiO3_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
        10.0 * dM_FeO_er + 10.0 * dM_MgO_er - 10.0 * dM_SiO2_er)) + M_MgO * (
        M_SiO2 * (-2 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 2 * dM_MgO_er - 4.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
        -4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)) + M_MgSiO3 * (
        M_FeO * (-6.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er) + M_MgO * (
        -6.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er) + M_SiO2 * (
        4.0 * dM_FeSiO3_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
        10.0 * dM_FeO_er + 10.0 * dM_MgO_er - 10.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
        6.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er + 6.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er)) + M_FeO * (M_MgO * (
        M_SiO2 * (-2 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 2 * dM_MgO_er - 4.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
        -4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                                 6.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er)) + M_FeSiO3 * (
                                                                                                          M_FeO * (
                                                                                                          M_MgO * (
                                                                                                          -6.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                          -6.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er) + M_m * (
                                                                                                          6.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er)) + M_MgO * (
                                                                                                          M_SiO2 * (
                                                                                                          -2.0 * dM_FeSiO3_er - 2 * dM_MgO_er - 4.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                                                                                                          4.0 * dM_FeO_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                                          6.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er)) + M_MgSiO3 * (
                                                                                                          M_FeO * (
                                                                                                          M_MgO * (
                                                                                                          -6.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                          4.0 * dM_FeO_er + 8.0 * dM_FeSiO3_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
                                                                                                          -10.0 * dM_FeSiO3_er + 10.0 * dM_MgO_er - 10.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                                                          M_FeO * (
                                                                                                          -6.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                                          -6.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                          4.0 * dM_FeSiO3_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
                                                                                                          10.0 * dM_FeO_er + 10.0 * dM_MgO_er - 10.0 * dM_SiO2_er)))) + dKMgO_KMgO * (
                                                                                                   M_O * (M_Fe * (
                                                                                                   M_MgO * (M_FeSiO3 * (
                                                                                                   4.0 * M_SiO2 - 10.0 * M_m) - 6.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                   M_FeO * (
                                                                                                   2.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                                                                                                   6.0 * M_SiO2 - 6.0 * M_m) - 6.0 * M_SiO2 * M_m)) + M_MgO * (
                                                                                                          -6.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                          M_FeO * (
                                                                                                          6.0 * M_SiO2 - 6.0 * M_m) - 6.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                          M_FeO * (
                                                                                                          M_MgO * (
                                                                                                          6.0 * M_SiO2 - 6.0 * M_m) - 6.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                          M_FeO * (
                                                                                                          6.0 * M_SiO2 - 6.0 * M_m) + M_MgO * (
                                                                                                          6.0 * M_SiO2 - 6.0 * M_m) - 6.0 * M_SiO2 * M_m))) + M_c * (
                                                                                                   M_Fe * (M_MgO * (
                                                                                                   M_FeSiO3 * (
                                                                                                   -M_SiO2 + 3.0 * M_m) + 2.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           -M_SiO2 - M_m) + M_MgO * (
                                                                                                           -2.0 * M_SiO2 + 2.0 * M_m) + 2.0 * M_SiO2 * M_m) + M_O * (
                                                                                                           M_MgO * (
                                                                                                           -M_FeSiO3 - M_SiO2) + M_MgSiO3 * (
                                                                                                           M_FeO - M_m))) + M_MgO * (
                                                                                                   2.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                   M_FeO * (
                                                                                                   -2.0 * M_SiO2 + 2.0 * M_m) + 2.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                   M_FeO * (M_MgO * (
                                                                                                   -2.0 * M_SiO2 + 2.0 * M_m) + 2.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                   M_FeO * (
                                                                                                   -2.0 * M_SiO2 + 2.0 * M_m) + M_MgO * (
                                                                                                   -2.0 * M_SiO2 + 2.0 * M_m) + 2.0 * M_SiO2 * M_m)) + M_O * (
                                                                                                   M_MgO * M_SiO2 * (
                                                                                                   -M_FeO - M_FeSiO3) + M_MgSiO3 * M_m * (
                                                                                                   -M_FeO - M_FeSiO3))))) + M_MgSiO3 * dKMgSiO3_KMgSiO3 * (
                       M_O * (M_Fe * M_MgO * (
                       M_FeO * (2.0 * M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (-4.0 * M_SiO2 + 10.0 * M_m)) + M_Mg * (M_Fe * (
                       M_FeO * (2.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                       2.0 * M_SiO2 + 4.0 * M_m) - 6.0 * M_SiO2 * M_m) + M_FeO * (M_MgO * (
                       2.0 * M_SiO2 + 4.0 * M_m) - 6.0 * M_SiO2 * M_m) + M_FeSiO3 * (M_FeO * (
                       6.0 * M_SiO2 - 6.0 * M_m) + M_MgO * (2.0 * M_SiO2 + 4.0 * M_m) - 6.0 * M_SiO2 * M_m))) + M_c * (
                       M_Mg * (M_Fe * (M_FeO * (-M_SiO2 - M_m) + M_MgO * (-M_SiO2 - M_m) + M_O * (
                       M_FeO + M_MgO - M_m) + 2.0 * M_SiO2 * M_m) + M_FeO * (
                               M_MgO * (-M_SiO2 - M_m) + 2.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                               M_FeO * (-2.0 * M_SiO2 + 2.0 * M_m) + M_MgO * (
                               -M_SiO2 - M_m) + 2.0 * M_SiO2 * M_m) + M_O * (
                               M_FeO * (M_MgO - M_m) + M_FeSiO3 * (M_MgO - M_m))) + M_MgO * (M_Fe * (
                       M_FeO * (-M_SiO2 - M_m) + M_FeSiO3 * (M_SiO2 - 3.0 * M_m) + M_O * (
                       M_FeO + M_FeSiO3 + M_SiO2 - M_m)) + M_O * (M_FeO * (M_SiO2 - M_m) + M_FeSiO3 * (
                       M_SiO2 - M_m))))) + M_c * (M_Fe * (M_MgO * (M_FeO * (
        M_SiO2 * (dM_FeO_er + 2.0 * dM_FeSiO3_er + dM_MgO_er + 2.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er) + M_m * (
        dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_FeSiO3 * (M_FeO * (
        2.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                 -dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                 -3.0 * dM_FeO_er + 3.0 * dM_MgSiO3_er + 3.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                   -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er)) + M_MgSiO3 * (
                                                          M_FeO * (M_MgO * (
                                                          2.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                   dM_FeO_er + 2.0 * dM_FeSiO3_er + 1.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er) + M_m * (
                                                                   dM_FeSiO3_er - dM_MgO_er + dM_SiO2_er)) + M_FeSiO3 * (
                                                          M_FeO * (
                                                          2.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er) + M_MgO * (
                                                          2.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                          -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                          -3.0 * dM_FeO_er - 3.0 * dM_MgO_er + 3.0 * dM_SiO2_er)) + M_MgO * (
                                                          M_SiO2 * (2.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er) + M_m * (
                                                          -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                          -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er)) + M_O * (M_MgO * (
        M_FeO * (-dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_FeSiO3 * (
        -dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_SiO2 * (
        -dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
        dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_MgSiO3 * (M_FeO * (
        -dM_FeSiO3_er + dM_MgO_er - dM_SiO2_er) + M_FeSiO3 * (-dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_SiO2 * (
                                                                 -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                 dM_FeSiO3_er - dM_MgO_er + dM_SiO2_er)))) + M_Mg * (
                                                  M_Fe * (M_FeO * (M_SiO2 * (
                                                  dM_FeO_er + 2.0 * dM_FeSiO3_er + dM_MgO_er + 2.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er) + M_m * (
                                                                   dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_FeSiO3 * (
                                                          M_FeO * (
                                                          2.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er) + M_MgO * (
                                                          2.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                          -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                          -3.0 * dM_FeO_er - 3.0 * dM_MgO_er + 3.0 * dM_SiO2_er)) + M_MgO * (
                                                          M_SiO2 * (
                                                          dM_FeO_er + 2.0 * dM_FeSiO3_er + dM_MgO_er + 2.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er) + M_m * (
                                                          dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_MgSiO3 * (
                                                          M_FeO * (
                                                          2.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er) + M_MgO * (
                                                          2.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                          -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                          -3.0 * dM_FeO_er - 3.0 * dM_MgO_er + 3.0 * dM_SiO2_er)) + M_O * (
                                                          M_FeO * (
                                                          -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_FeSiO3 * (
                                                          -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_MgO * (
                                                          -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_MgSiO3 * (
                                                          -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_SiO2 * (
                                                          -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                          dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_SiO2 * M_m * (
                                                          -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er)) + M_FeO * (
                                                  M_MgO * (M_SiO2 * (
                                                  dM_FeO_er + 2.0 * dM_FeSiO3_er + dM_MgO_er + 2.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er) + M_m * (
                                                           dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_SiO2 * M_m * (
                                                  -2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er)) + M_FeSiO3 * (M_FeO * (
                                                  M_MgO * (
                                                  2.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                  2.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er) + M_m * (
                                                  -2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er)) + M_MgO * (M_SiO2 * (
                                                  1.0 * dM_FeSiO3_er + dM_MgO_er + 2.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er) + M_m * (
                                                                                                     -dM_FeO_er + dM_MgSiO3_er + dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                                        -2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er)) + M_MgSiO3 * (
                                                  M_FeO * (M_MgO * (
                                                  2.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                           -dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                           3.0 * dM_FeSiO3_er - 3.0 * dM_MgO_er + 3.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                  M_FeO * (
                                                  2.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er) + M_MgO * (
                                                  2.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                  -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                  -3.0 * dM_FeO_er - 3.0 * dM_MgO_er + 3.0 * dM_SiO2_er))) + M_O * (
                                                  M_FeO * (
                                                  M_MgO * (-dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_SiO2 * (
                                                  -dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                  dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_FeSiO3 * (
                                                  M_MgO * (dM_FeO_er - dM_MgSiO3_er - dM_SiO2_er) + M_SiO2 * (
                                                  -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                  -dM_FeO_er + dM_MgSiO3_er + dM_SiO2_er)) + M_MgSiO3 * (M_FeO * (
                                                  -dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_FeSiO3 * (
                                                                                                         -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er)))) + M_O * (
                                                  M_MgO * (M_FeO * (M_SiO2 * (
                                                  -dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                    dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_FeSiO3 * (
                                                           M_SiO2 * (
                                                           -dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                           -dM_FeO_er + dM_MgSiO3_er + dM_SiO2_er))) + M_MgSiO3 * (
                                                  M_FeO * (M_SiO2 * (
                                                  -dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                           dM_FeSiO3_er - dM_MgO_er + dM_SiO2_er)) + M_FeSiO3 * (
                                                  M_SiO2 * (-dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                  -dM_FeO_er - dM_MgO_er + dM_SiO2_er))))) + dKSiO2_KSiO2 * (M_O * (
        M_Fe * (M_MgO * (4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
        M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                M_FeO * (M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m))) + M_Mg * (M_Fe * (M_FeSiO3 * (
        M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                     M_FeO * (
                                                                                     -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                                                                                     -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                                                                                     4.0 * M_FeO + 4.0 * M_MgO)) + M_MgO * (
                                                                             4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                             M_FeO * (
                                                                             -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                             M_FeO * (M_MgO * (
                                                                             -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                             M_FeO * (
                                                                             -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                                                                             -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)))) + M_c * (
                                                                                                             M_Fe * (
                                                                                                             M_MgO * (
                                                                                                             -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                             M_FeO * (
                                                                                                             M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                             M_FeO * (
                                                                                                             M_MgO * (
                                                                                                             M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                             M_FeO * (
                                                                                                             M_SiO2 - M_m) + M_MgO * (
                                                                                                             M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
                                                                                                             M_MgO * (
                                                                                                             M_FeSiO3 * (
                                                                                                             M_FeO - M_m) + M_SiO2 * (
                                                                                                             M_FeO - M_m)) + M_MgSiO3 * (
                                                                                                             M_FeO * (
                                                                                                             M_MgO + M_SiO2) + M_FeSiO3 * (
                                                                                                             M_FeO + M_MgO - M_m) + M_MgO * (
                                                                                                             M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Mg * (
                                                                                                             M_Fe * (
                                                                                                             M_FeSiO3 * (
                                                                                                             M_FeO * (
                                                                                                             M_SiO2 - M_m) + M_MgO * (
                                                                                                             M_SiO2 - M_m) - M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                             M_FeO * (
                                                                                                             M_SiO2 - M_m) + M_MgO * (
                                                                                                             M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
                                                                                                             M_FeSiO3 * (
                                                                                                             M_FeO + M_MgO - M_m) + M_MgSiO3 * (
                                                                                                             M_FeO + M_MgO - M_m) + M_SiO2 * (
                                                                                                             M_FeO + M_MgO - M_m)) + M_SiO2 * M_m * (
                                                                                                             -M_FeO - M_MgO)) + M_MgO * (
                                                                                                             -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                             M_FeO * (
                                                                                                             M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                             M_FeO * (
                                                                                                             M_MgO * (
                                                                                                             M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                             M_FeO * (
                                                                                                             M_SiO2 - M_m) + M_MgO * (
                                                                                                             M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
                                                                                                             M_FeO * M_SiO2 * (
                                                                                                             M_MgO - M_m) + M_FeSiO3 * (
                                                                                                             M_FeO * (
                                                                                                             M_MgO + M_SiO2 - M_m) + M_SiO2 * (
                                                                                                             M_MgO - M_m)) + M_MgSiO3 * (
                                                                                                             M_FeO * (
                                                                                                             M_MgO - M_m) + M_FeSiO3 * (
                                                                                                             M_FeO + M_MgO - M_m)))) + M_O * (
                                                                                                             M_MgO * (
                                                                                                             -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                             M_FeO * (
                                                                                                             M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                             M_FeO * (
                                                                                                             M_MgO * (
                                                                                                             M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                             M_FeO * (
                                                                                                             M_SiO2 - M_m) + M_MgO * (
                                                                                                             M_SiO2 - M_m) - M_SiO2 * M_m)))))) / (
               M_O * (M_Fe * (M_MgO * (4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                              M_FeO * (M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                              M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                              -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m))) + M_Mg * (M_Fe * (M_FeSiO3 * (
               M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
               -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO)) + M_MgO * (
                                                                                           4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (M_MgO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)))) + M_Si * (
               M_Fe * (
               M_MgO * (M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (M_MgO * (
               M_FeO * (-M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
               -9.0 * M_FeO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               -9.0 * M_MgO - M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                                                                                             -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                                                                                             -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m))) + M_Mg * (
               M_Fe * (M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                       M_FeO * (-M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                       -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                       -M_SiO2 + 4.0 * M_m) + M_MgSiO3 * (
                       -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                       M_FeO + M_MgO)) + M_MgO * (
               M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (
               M_FeO * (M_MgO * (-M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-9.0 * M_MgO - 9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (
               -M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (-9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_FeSiO3 * (
               -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m)))) + M_O * (M_MgO * (
               9.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)) + M_MgSiO3 * (M_FeO * (
               M_MgO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (M_FeO * (
               -9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)))) + M_c * (
               M_Fe * (
               M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
               M_MgO * (M_FeSiO3 * (M_FeO - M_m) + M_SiO2 * (M_FeO - M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO + M_SiO2) + M_FeSiO3 * (M_FeO + M_MgO - M_m) + M_MgO * (
               M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Mg * (M_Fe * (
               M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
               M_FeSiO3 * (M_FeO + M_MgO - M_m) + M_MgSiO3 * (M_FeO + M_MgO - M_m) + M_SiO2 * (
               M_FeO + M_MgO - M_m)) + M_SiO2 * M_m * (-M_FeO - M_MgO)) + M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                         M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                         M_FeO * (M_SiO2 - M_m) + M_MgO * (
                                                         M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
                                                         M_FeO * M_SiO2 * (M_MgO - M_m) + M_FeSiO3 * (
                                                         M_FeO * (M_MgO + M_SiO2 - M_m) + M_SiO2 * (
                                                         M_MgO - M_m)) + M_MgSiO3 * (
                                                         M_FeO * (M_MgO - M_m) + M_FeSiO3 * (
                                                         M_FeO + M_MgO - M_m)))) + M_O * (
               M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Si * (M_Fe * (M_MgO * (
               M_FeO * (M_SiO2 - M_m) + M_FeSiO3 * (
               4.0 * M_FeO + M_SiO2 - 9.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               4.0 * M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                     4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_O * (
                                                                                                   M_MgO * (
                                                                                                   M_FeO + M_FeSiO3 + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                   M_FeO + M_FeSiO3 + M_SiO2 - M_m))) + M_Mg * (
                                                                                           M_Fe * (M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                                   M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_O * (
                                                                                                   M_FeO + M_FeSiO3 + M_MgO + M_MgSiO3 + M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeO * (
                                                                                           M_MgO * (
                                                                                           M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_MgO + 4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                           M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_FeSiO3 * (
                                                                                           4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m)) + M_O * (
                                                                                           M_FeO * (
                                                                                           M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                           M_MgO + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                           M_FeO + M_FeSiO3))) + M_MgO * (
                                                                                           -4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (M_MgO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_O * (
                                                                                           M_MgO * (M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                    M_SiO2 - M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                           M_SiO2 - M_m))))))

    def dM_Fe_dTc(self, Moles, dKs, dMi_b):
        dM_Mg_er, dM_Si_er, dM_Fe_er, dM_O_er, dM_c_er, dM_MgO_er, dM_SiO2_er, dM_FeO_er, dM_MgSiO3_er, dM_FeSiO3_er, dM_m_er = self.unwrap_Moles(
            dMi_b)
        M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = Moles
        dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs
        return M_Fe * (M_FeSiO3 * dKFeSiO3_KFeSiO3 * (M_Mg * M_O * (4.0 * M_MgO * M_SiO2 * M_m + M_MgSiO3 * (
        M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_Si * (
                                                      M_Mg * (M_MgO * M_SiO2 * M_m + M_MgSiO3 * (
                                                      M_FeO * (-M_SiO2 + M_m) + M_MgO * (
                                                      -M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                                                              M_FeO * (-3.0 * M_SiO2 - 6.0 * M_m) + M_MgO * (
                                                              -M_SiO2 + 4.0 * M_m) + M_MgSiO3 * (
                                                              -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m)) + M_O * (
                                                      M_MgO * (M_FeO * (
                                                      -3.0 * M_SiO2 - 6.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                      M_FeO * (-3.0 * M_SiO2 - 6.0 * M_m) + M_MgO * (
                                                      -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m))) + M_c * (
                                                      M_Mg * (-M_MgO * M_SiO2 * M_m + M_MgSiO3 * (
                                                      M_FeO * (M_SiO2 - M_m) + M_MgO * (
                                                      M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
                                                              M_MgSiO3 * (M_FeO + M_MgO - M_m) + M_SiO2 * (
                                                              M_FeO + M_MgO - M_m))) + M_O * (
                                                      M_MgO * M_SiO2 * (M_FeO - M_m) + M_MgSiO3 * (
                                                      M_MgO * (M_SiO2 - M_m) + M_SiO2 * (M_FeO - M_m))) + M_Si * (
                                                      M_Mg * (M_FeO * (2.0 * M_SiO2 + 2.0 * M_m) + M_MgO * (
                                                      M_SiO2 - M_m) + M_MgSiO3 * (
                                                              4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_O * (
                                                              M_MgO + M_MgSiO3 + M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_MgO * (
                                                      M_FeO * (
                                                      2.0 * M_SiO2 + 2.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                      M_FeO * (2.0 * M_SiO2 + 2.0 * M_m) + M_MgO * (
                                                      4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_O * (
                                                      M_MgO * (M_SiO2 - M_m) + M_MgSiO3 * (M_SiO2 - M_m))))) + M_Mg * (
                       M_O * (M_FeSiO3 * (M_FeO * (M_SiO2 * (-4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_m * (
                       4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er)) + M_MgO * (M_SiO2 * (
                       -4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 8.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_m * (
                                                                         -4.0 * dM_FeO_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                          4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er)) + M_MgSiO3 * (M_FeO * (
                       M_SiO2 * (4.0 * dM_FeO_er + 8.0 * dM_FeSiO3_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
                       -4.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er - 4.0 * dM_SiO2_er)) + M_MgO * (M_SiO2 * (
                       4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er) + M_m * (
                                                                                             -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                               -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                              M_FeO * (4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er) + M_MgO * (
                              -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er))) + dKMgO_KMgO * (M_O * (M_MgO * (
                       -4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                       M_FeO * (4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (M_FeO * (
                       M_MgO * (4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (M_FeO * (
                       4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m))) + M_Si * (
                                                                                       M_MgO * (
                                                                                       -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                       M_FeO * (
                                                                                       M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                       M_FeO * (M_MgO * (
                                                                                       M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                       M_FeO * (
                                                                                       M_SiO2 - M_m) + M_MgO * (
                                                                                       M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
                                                                                       M_MgO * (M_FeO * (
                                                                                       M_SiO2 - 4.0 * M_m) + M_FeSiO3 * (
                                                                                                9.0 * M_FeO - 2.0 * M_SiO2 - 10.0 * M_m)) + M_MgSiO3 * (
                                                                                       M_FeO * (
                                                                                       9.0 * M_MgO - 2.0 * M_SiO2 - 10.0 * M_m) + M_FeSiO3 * (
                                                                                       9.0 * M_FeO + 9.0 * M_MgO + 4.0 * M_SiO2 - 25.0 * M_m)))) + M_c * (
                                                                                       M_MgO * (
                                                                                       M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                       M_FeO * (
                                                                                       -M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                       M_FeO * (M_MgO * (
                                                                                       -M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                       M_FeO * (
                                                                                       -M_SiO2 + M_m) + M_MgO * (
                                                                                       -M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (
                                                                                       M_FeO * M_MgO * (
                                                                                       -M_FeSiO3 - M_SiO2) + M_MgSiO3 * (
                                                                                       -M_FeO * M_MgO + M_FeSiO3 * (
                                                                                       -M_FeO - M_MgO + M_m))) + M_Si * (
                                                                                       M_MgO * (M_FeO * (
                                                                                       -M_SiO2 + M_m) + M_FeSiO3 * (
                                                                                                -4.0 * M_FeO + M_SiO2 + 3.0 * M_m)) + M_MgSiO3 * (
                                                                                       M_FeO * (
                                                                                       -4.0 * M_MgO + M_SiO2 + 3.0 * M_m) + M_FeSiO3 * (
                                                                                       -4.0 * M_FeO - 4.0 * M_MgO - M_SiO2 + 9.0 * M_m)) + M_O * (
                                                                                       M_MgO * (
                                                                                       -M_FeO - M_FeSiO3) + M_MgSiO3 * (
                                                                                       -M_FeO - M_FeSiO3)))))) + M_MgSiO3 * dKMgSiO3_KMgSiO3 * (
                       M_Mg * M_O * (-4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                       M_FeO * (4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                       4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_Si * (M_Mg * (
                       -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                       M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
                       M_FeO * (-2.0 * M_SiO2 - 10.0 * M_m) + M_FeSiO3 * (
                       9.0 * M_FeO + 9.0 * M_MgO + 4.0 * M_SiO2 - 25.0 * M_m))) + M_MgO * M_O * (M_FeO * (
                       -3.0 * M_SiO2 - 6.0 * M_m) + M_FeSiO3 * (6.0 * M_SiO2 - 15.0 * M_m))) + M_c * (M_Mg * (
                       M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_O * (
                       -M_FeO - M_MgO + M_m) + M_SiO2 * M_m)) + M_MgO * M_O * (
                                                                                                      M_FeO * M_SiO2 + M_FeSiO3 * M_m) + M_Si * (
                                                                                                      M_Mg * (M_FeO * (
                                                                                                      M_SiO2 + 3.0 * M_m) + M_FeSiO3 * (
                                                                                                              -4.0 * M_FeO - 4.0 * M_MgO - M_SiO2 + 9.0 * M_m) + M_O * (
                                                                                                              -M_FeO - M_FeSiO3)) + M_MgO * (
                                                                                                      M_FeO * (
                                                                                                      2.0 * M_SiO2 + 2.0 * M_m) + M_FeSiO3 * (
                                                                                                      -2.0 * M_SiO2 + 6.0 * M_m))))) + M_Si * (
                       M_Mg * (M_FeSiO3 * (M_FeO * (M_SiO2 * (-1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_m * (
                       1.0 * dM_MgO_er + 1.0 * dM_MgSiO3_er)) + M_MgO * (M_SiO2 * (
                       -dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                         -1.0 * dM_FeO_er + 1.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                           1.0 * dM_MgO_er + 1.0 * dM_MgSiO3_er)) + M_MgSiO3 * (M_FeO * (
                       M_SiO2 * (dM_FeO_er + 2.0 * dM_FeSiO3_er + 1.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er) + M_m * (
                       -1.0 * dM_FeSiO3_er + 1.0 * dM_MgO_er - 1.0 * dM_SiO2_er)) + M_MgO * (M_SiO2 * (
                       dM_FeO_er + dM_FeSiO3_er) + M_m * (-dM_FeO_er - dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                                -dM_FeO_er - dM_FeSiO3_er)) + M_O * (
                               M_FeO * (M_SiO2 * (
                               3.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er + 2 * dM_MgO_er + 5.0 * dM_MgSiO3_er + 3.0 * dM_SiO2_er) + M_m * (
                                        6.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 10.0 * dM_MgSiO3_er + 6.0 * dM_SiO2_er)) + M_FeSiO3 * (
                               M_FeO * (9.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er) + M_MgO * (
                               3.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 15.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er) + M_SiO2 * (
                               -6.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 10.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er) + M_m * (
                               -15.0 * dM_FeO_er + 10.0 * dM_MgO_er + 25.0 * dM_MgSiO3_er + 15.0 * dM_SiO2_er)) + M_MgO * (
                               M_SiO2 * (dM_FeO_er + dM_FeSiO3_er) + M_m * (
                               -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er)) + M_MgSiO3 * (M_FeO * (
                               15.0 * dM_FeO_er + 24.0 * dM_FeSiO3_er + 6.0 * dM_MgO_er + 15.0 * dM_MgSiO3_er + 9.0 * dM_SiO2_er) + M_MgO * (
                                                                                     9.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er) + M_SiO2 * (
                                                                                     4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er) + M_m * (
                                                                                     -25.0 * dM_FeO_er - 25.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                               -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                               M_FeO * (1.0 * dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_MgO * (
                               -dM_FeO_er - dM_FeSiO3_er))) + M_O * (M_MgO * (M_FeO * (M_SiO2 * (
                       3.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er + 3.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 3.0 * dM_SiO2_er) + M_m * (
                                                                                       6.0 * dM_FeSiO3_er + 6.0 * dM_MgSiO3_er + 6.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                              M_FeO * (
                                                                              9.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er + 9.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                              -6.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 12.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er) + M_m * (
                                                                              -15.0 * dM_FeO_er + 15.0 * dM_MgSiO3_er + 15.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                              -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er)) + M_MgSiO3 * (
                                                                     M_FeO * (M_MgO * (
                                                                     9.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er + 9.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                              3.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er + 3.0 * dM_MgSiO3_er + 3.0 * dM_SiO2_er) + M_m * (
                                                                              6.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er + 6.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                     M_FeO * (
                                                                     9.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er + 9.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er) + M_MgO * (
                                                                     9.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er + 9.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                     -6.0 * dM_FeSiO3_er - 6.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er) + M_m * (
                                                                     -15.0 * dM_FeO_er - 15.0 * dM_MgO_er + 15.0 * dM_SiO2_er)) + M_MgO * (
                                                                     M_SiO2 * (
                                                                     9.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er) + M_m * (
                                                                     -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                     -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er))) + dKSiO2_KSiO2 * (
                       M_O * (M_Mg * (
                       M_FeO * (M_MgSiO3 * (4.0 * M_SiO2 - 10.0 * M_m) - 6.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                       M_FeO * (6.0 * M_SiO2 - 6.0 * M_m) + M_MgO * (
                       2.0 * M_SiO2 + 4.0 * M_m) - 6.0 * M_SiO2 * M_m)) + M_MgO * (
                              -6.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                              M_FeO * (6.0 * M_SiO2 - 6.0 * M_m) - 6.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                              M_FeO * (M_MgO * (6.0 * M_SiO2 - 6.0 * M_m) - 6.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                              M_FeO * (6.0 * M_SiO2 - 6.0 * M_m) + M_MgO * (
                              6.0 * M_SiO2 - 6.0 * M_m) - 6.0 * M_SiO2 * M_m))) + M_c * (M_Mg * (
                       M_FeO * (M_MgSiO3 * (-M_SiO2 + 3.0 * M_m) + 2.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                       M_FeO * (-2.0 * M_SiO2 + 2.0 * M_m) + M_MgO * (-M_SiO2 - M_m) + 2.0 * M_SiO2 * M_m) + M_O * (
                       M_FeO * (-M_MgSiO3 - M_SiO2) + M_FeSiO3 * (M_MgO - M_m))) + M_MgO * (
                                                                                         2.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                         M_FeO * (
                                                                                         -2.0 * M_SiO2 + 2.0 * M_m) + 2.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                         M_FeO * (M_MgO * (
                                                                                         -2.0 * M_SiO2 + 2.0 * M_m) + 2.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                         M_FeO * (
                                                                                         -2.0 * M_SiO2 + 2.0 * M_m) + M_MgO * (
                                                                                         -2.0 * M_SiO2 + 2.0 * M_m) + 2.0 * M_SiO2 * M_m)) + M_O * (
                                                                                         M_MgO * (
                                                                                         -M_FeO * M_SiO2 - M_FeSiO3 * M_m) + M_MgSiO3 * (
                                                                                         -M_FeO * M_SiO2 - M_FeSiO3 * M_m))))) + M_c * (
                       M_Mg * (M_FeSiO3 * (M_FeO * (
                       M_SiO2 * (dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_m * (-dM_MgO_er - 1.0 * dM_MgSiO3_er)) + M_MgO * (
                                           M_SiO2 * (
                                           1.0 * dM_FeSiO3_er + dM_MgO_er + 2.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er) + M_m * (
                                           dM_FeO_er - dM_MgSiO3_er - dM_SiO2_er)) + M_SiO2 * M_m * (
                                           -dM_MgO_er - 1.0 * dM_MgSiO3_er)) + M_MgSiO3 * (M_FeO * (
                       M_SiO2 * (-dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                       1.0 * dM_FeSiO3_er - dM_MgO_er + 1.0 * dM_SiO2_er)) + M_MgO * (M_SiO2 * (
                       -dM_FeO_er - dM_FeSiO3_er) + M_m * (dM_FeO_er + dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                           dM_FeO_er + dM_FeSiO3_er)) + M_O * (
                               M_FeSiO3 * (M_FeO * (-dM_FeO_er - dM_FeSiO3_er) + M_MgO * (
                               -dM_FeO_er + dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                                           dM_FeO_er - dM_MgSiO3_er - dM_SiO2_er)) + M_MgSiO3 * (
                               M_FeO * (-dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_MgO * (
                               -dM_FeO_er - dM_FeSiO3_er) + M_m * (dM_FeO_er + dM_FeSiO3_er)) + M_SiO2 * (
                               M_FeO * (-dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_MgO * (
                               -dM_FeO_er - dM_FeSiO3_er) + M_m * (dM_FeO_er + dM_FeSiO3_er))) + M_SiO2 * M_m * (
                               M_FeO * (-dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_MgO * (
                               dM_FeO_er + dM_FeSiO3_er))) + M_O * (M_MgO * (M_FeSiO3 * (
                       M_FeO * (-dM_FeO_er - dM_FeSiO3_er - dM_MgO_er - dM_MgSiO3_er) + M_m * (
                       dM_FeO_er - dM_MgSiO3_er - dM_SiO2_er)) + M_SiO2 * (M_FeO * (
                       -dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                           dM_FeO_er + dM_FeSiO3_er))) + M_MgSiO3 * (
                                                                    M_FeO * (M_MgO * (
                                                                    -dM_FeO_er - dM_FeSiO3_er - dM_MgO_er - dM_MgSiO3_er) + M_SiO2 * (
                                                                             -dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er)) + M_FeSiO3 * (
                                                                    M_FeO * (
                                                                    -dM_FeO_er - dM_FeSiO3_er - dM_MgO_er - dM_MgSiO3_er) + M_MgO * (
                                                                    -dM_FeO_er - dM_FeSiO3_er - dM_MgO_er - dM_MgSiO3_er) + M_m * (
                                                                    dM_FeO_er + dM_MgO_er - dM_SiO2_er)) + M_MgO * (
                                                                    M_SiO2 * (-dM_FeO_er - dM_FeSiO3_er) + M_m * (
                                                                    dM_FeO_er + dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                    dM_FeO_er + dM_FeSiO3_er))) + M_Si * (M_Mg * (
                       M_FeO * (M_SiO2 * (
                       -2 * dM_FeO_er - 4.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                                -2.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er)) + M_FeSiO3 * (
                       M_FeO * (-4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er) + M_MgO * (
                       -2.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_SiO2 * (
                       2.0 * dM_FeSiO3_er + dM_MgO_er + 3.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_m * (
                       6.0 * dM_FeO_er - 3.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er)) + M_MgO * (
                       M_SiO2 * (-dM_FeO_er - dM_FeSiO3_er) + M_m * (dM_FeO_er + dM_FeSiO3_er)) + M_MgSiO3 * (M_FeO * (
                       -6.0 * dM_FeO_er - 10.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_MgO * (
                                                                                                              -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er) + M_SiO2 * (
                                                                                                              -dM_FeO_er - dM_FeSiO3_er) + M_m * (
                                                                                                              9.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er)) + M_O * (
                       M_FeO * (-dM_FeO_er - dM_FeSiO3_er) + M_FeSiO3 * (-dM_FeO_er - dM_FeSiO3_er) + M_MgO * (
                       -dM_FeO_er - dM_FeSiO3_er) + M_MgSiO3 * (-dM_FeO_er - dM_FeSiO3_er) + M_SiO2 * (
                       -dM_FeO_er - dM_FeSiO3_er) + M_m * (dM_FeO_er + dM_FeSiO3_er)) + M_SiO2 * M_m * (
                       4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er)) + M_MgO * (M_FeO * (M_SiO2 * (
                       -2 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 2 * dM_MgO_er - 4.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                                                                                  -2.0 * dM_FeSiO3_er - 2.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                         M_FeO * (
                                                                         -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                         2.0 * dM_FeSiO3_er + 2 * dM_MgO_er + 4.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_m * (
                                                                         6.0 * dM_FeO_er - 6.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                         4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er)) + M_MgSiO3 * (
                                                                                                          M_FeO * (
                                                                                                          M_MgO * (
                                                                                                          -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                          -2 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 2.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                                                                                                          -2.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er - 2.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                                                          M_FeO * (
                                                                                                          -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                                          -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                          2.0 * dM_FeSiO3_er + 2.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_m * (
                                                                                                          6.0 * dM_FeO_er + 6.0 * dM_MgO_er - 6.0 * dM_SiO2_er)) + M_MgO * (
                                                                                                          M_SiO2 * (
                                                                                                          -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er) + M_m * (
                                                                                                          4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                                          4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er)) + M_O * (
                                                                                                          M_MgO * (
                                                                                                          M_FeO * (
                                                                                                          -dM_FeO_er - dM_FeSiO3_er - dM_MgO_er - dM_MgSiO3_er) + M_FeSiO3 * (
                                                                                                          -dM_FeO_er - dM_FeSiO3_er - dM_MgO_er - dM_MgSiO3_er) + M_SiO2 * (
                                                                                                          -dM_FeO_er - dM_FeSiO3_er) + M_m * (
                                                                                                          dM_FeO_er + dM_FeSiO3_er)) + M_MgSiO3 * (
                                                                                                          M_FeO * (
                                                                                                          -dM_FeO_er - dM_FeSiO3_er - dM_MgO_er - dM_MgSiO3_er) + M_FeSiO3 * (
                                                                                                          -dM_FeO_er - dM_FeSiO3_er - dM_MgO_er - dM_MgSiO3_er) + M_SiO2 * (
                                                                                                          -dM_FeO_er - dM_FeSiO3_er) + M_m * (
                                                                                                          dM_FeO_er + dM_FeSiO3_er))))) + dKFeO_KFeO * (
                       M_Mg * M_O * (M_MgO * (4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                       M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                     M_FeO * (M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                     M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                                     -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m))) + M_Si * (M_Mg * (M_MgO * (
                       M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                          M_FeO * (
                                                                                                          M_MgO * (
                                                                                                          -M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                          M_FeO * (
                                                                                                          -M_SiO2 + M_m) + M_MgO * (
                                                                                                          -M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (
                                                                                                          M_FeO * (
                                                                                                          M_MgO * (
                                                                                                          -M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                          M_FeO * (
                                                                                                          -9.0 * M_MgO - 9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (
                                                                                                          -M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                          M_FeO * (
                                                                                                          -9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_FeSiO3 * (
                                                                                                          -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m)))) + M_O * (
                                                                                                  M_MgO * (
                                                                                                  9.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                  M_FeO * (
                                                                                                  -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                  M_FeO * (M_MgO * (
                                                                                                  -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                  M_FeO * (
                                                                                                  -9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (
                                                                                                  -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)))) + M_c * (
                       M_Mg * (M_MgO * (
                       -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                               M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
                               M_FeO * M_SiO2 * (M_MgO - M_m) + M_FeSiO3 * (
                               M_FeO * (M_MgO + M_SiO2 - M_m) + M_SiO2 * (M_MgO - M_m)) + M_MgSiO3 * (
                               M_FeO * (M_MgO - M_m) + M_FeSiO3 * (M_FeO + M_MgO - M_m)))) + M_O * (M_MgO * (
                       -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                    M_FeO * (M_MgO * (
                                                                                                    M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                    M_FeO * (
                                                                                                    M_SiO2 - M_m) + M_MgO * (
                                                                                                    M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Si * (
                       M_Mg * (M_FeO * (M_MgO * (M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                       M_FeO * (4.0 * M_MgO + 4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                       M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                               M_FeO * (4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_FeSiO3 * (
                               4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m)) + M_O * (
                               M_FeO * (M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (M_MgO + M_SiO2 - M_m) + M_MgSiO3 * (
                               M_FeO + M_FeSiO3))) + M_MgO * (-4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                       M_FeO * (4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                       M_FeO * (M_MgO * (4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                       M_FeO * (4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                       4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_O * (
                       M_MgO * (M_FeO * (M_SiO2 - M_m) + M_FeSiO3 * (M_SiO2 - M_m)) + M_MgSiO3 * (
                       M_FeO * (M_SiO2 - M_m) + M_FeSiO3 * (M_SiO2 - M_m))))))) / (M_O * (M_Fe * (M_MgO * (
        4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
        M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (M_FeO * (
        M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (M_FeO * (
        -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m))) + M_Mg * (M_Fe * (
        M_FeSiO3 * (
        M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
        M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
        -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (4.0 * M_FeO + 4.0 * M_MgO)) + M_MgO * (
                                                                                                           4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           M_MgO * (
                                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)))) + M_Si * (
                                                                                   M_Fe * (M_MgO * (
                                                                                   M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                   M_FeO * (
                                                                                   -M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (M_MgO * (
                                                                                           -M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           -M_SiO2 + M_m) + M_MgO * (
                                                                                           -M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (
                                                                                           M_MgO * (M_FeO * (
                                                                                           -M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                                                                                                    -9.0 * M_FeO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                           M_FeO * (
                                                                                           -9.0 * M_MgO - M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                                                                                           -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                                                                                           -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m))) + M_Mg * (
                                                                                   M_Fe * (M_FeSiO3 * (
                                                                                   M_FeO * (-M_SiO2 + M_m) + M_MgO * (
                                                                                   -M_SiO2 + M_m) + M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                           M_FeO * (
                                                                                           -M_SiO2 + M_m) + M_MgO * (
                                                                                           -M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                                                                                           M_FeO * (
                                                                                           -M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                                                                                           -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                                                                                           -M_SiO2 + 4.0 * M_m) + M_MgSiO3 * (
                                                                                           -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                                                                                           M_FeO + M_MgO)) + M_MgO * (
                                                                                   M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                   M_FeO * (
                                                                                   -M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                   M_FeO * (M_MgO * (
                                                                                   -M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                   M_FeO * (-M_SiO2 + M_m) + M_MgO * (
                                                                                   -M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (
                                                                                   M_FeO * (M_MgO * (
                                                                                   -M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                   M_FeO * (
                                                                                   -9.0 * M_MgO - 9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (
                                                                                   -M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                   M_FeO * (
                                                                                   -9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_FeSiO3 * (
                                                                                   -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m)))) + M_O * (
                                                                                   M_MgO * (
                                                                                   9.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                   M_FeO * (
                                                                                   -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                   M_FeO * (M_MgO * (
                                                                                   -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                   M_FeO * (
                                                                                   -9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (
                                                                                   -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)))) + M_c * (
                                                                                   M_Fe * (M_MgO * (
                                                                                   -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                   M_FeO * (
                                                                                   M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (M_MgO * (
                                                                                           M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           M_SiO2 - M_m) + M_MgO * (
                                                                                           M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
                                                                                           M_MgO * (M_FeSiO3 * (
                                                                                           M_FeO - M_m) + M_SiO2 * (
                                                                                                    M_FeO - M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (
                                                                                           M_MgO + M_SiO2) + M_FeSiO3 * (
                                                                                           M_FeO + M_MgO - M_m) + M_MgO * (
                                                                                           M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Mg * (
                                                                                   M_Fe * (M_FeSiO3 * (
                                                                                   M_FeO * (M_SiO2 - M_m) + M_MgO * (
                                                                                   M_SiO2 - M_m) - M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                           M_FeO * (
                                                                                           M_SiO2 - M_m) + M_MgO * (
                                                                                           M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
                                                                                           M_FeSiO3 * (
                                                                                           M_FeO + M_MgO - M_m) + M_MgSiO3 * (
                                                                                           M_FeO + M_MgO - M_m) + M_SiO2 * (
                                                                                           M_FeO + M_MgO - M_m)) + M_SiO2 * M_m * (
                                                                                           -M_FeO - M_MgO)) + M_MgO * (
                                                                                   -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                   M_FeO * (
                                                                                   M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                   M_FeO * (M_MgO * (
                                                                                   M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                   M_FeO * (M_SiO2 - M_m) + M_MgO * (
                                                                                   M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
                                                                                   M_FeO * M_SiO2 * (
                                                                                   M_MgO - M_m) + M_FeSiO3 * (M_FeO * (
                                                                                   M_MgO + M_SiO2 - M_m) + M_SiO2 * (
                                                                                                              M_MgO - M_m)) + M_MgSiO3 * (
                                                                                   M_FeO * (M_MgO - M_m) + M_FeSiO3 * (
                                                                                   M_FeO + M_MgO - M_m)))) + M_O * (
                                                                                   M_MgO * (
                                                                                   -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                   M_FeO * (
                                                                                   M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                   M_FeO * (M_MgO * (
                                                                                   M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                   M_FeO * (M_SiO2 - M_m) + M_MgO * (
                                                                                   M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Si * (
                                                                                   M_Fe * (M_MgO * (
                                                                                   M_FeO * (M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                   4.0 * M_FeO + M_SiO2 - 9.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                           4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_O * (
                                                                                           M_MgO * (
                                                                                           M_FeO + M_FeSiO3 + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                           M_FeO + M_FeSiO3 + M_SiO2 - M_m))) + M_Mg * (
                                                                                   M_Fe * (
                                                                                   M_FeO * (M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                   4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                   M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                   4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_O * (
                                                                                   M_FeO + M_FeSiO3 + M_MgO + M_MgSiO3 + M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeO * (
                                                                                   M_MgO * (
                                                                                   M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                   M_FeO * (
                                                                                   4.0 * M_MgO + 4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                   M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                   M_FeO * (
                                                                                   4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_FeSiO3 * (
                                                                                   4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m)) + M_O * (
                                                                                   M_FeO * (
                                                                                   M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                   M_MgO + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                   M_FeO + M_FeSiO3))) + M_MgO * (
                                                                                   -4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                   M_FeO * (
                                                                                   4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                   M_FeO * (M_MgO * (
                                                                                   4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                   M_FeO * (
                                                                                   4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                   4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_O * (
                                                                                   M_MgO * (
                                                                                   M_FeO * (M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                   M_SiO2 - M_m)) + M_MgSiO3 * (
                                                                                   M_FeO * (M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                   M_SiO2 - M_m))))))

    def dM_Mg_dTc(self, Moles, dKs, dMi_b):
        dM_Mg_er, dM_Si_er, dM_Fe_er, dM_O_er, dM_c_er, dM_MgO_er, dM_SiO2_er, dM_FeO_er, dM_MgSiO3_er, dM_FeSiO3_er, dM_m_er = self.unwrap_Moles(
            dMi_b)
        M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = Moles
        dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs
        return M_Mg * (M_Fe * (M_O * (M_FeSiO3 * (M_FeO * (
        M_SiO2 * (4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er) + M_m * (-4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er)) + M_MgO * (
                                                  M_SiO2 * (
                                                  4.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 8.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_m * (
                                                  4.0 * dM_FeO_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                  -4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er)) + M_MgSiO3 * (M_FeO * (
        M_SiO2 * (-4.0 * dM_FeO_er - 8.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_m * (
        4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er + 4.0 * dM_SiO2_er)) + M_MgO * (M_SiO2 * (
        -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er) + M_m * (4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                                        4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                      M_FeO * (-4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_MgO * (
                                      4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er))) + dKFeO_KFeO * (M_O * (M_MgO * (
        -4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
        M_FeO * (4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (M_FeO * (
        M_MgO * (4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (M_FeO * (
        4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m))) + M_Si * (M_MgO * (
        -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (M_FeO * (
        M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) + M_MgO * (
        M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (M_MgO * (
        M_FeO * (M_SiO2 - 4.0 * M_m) + M_FeSiO3 * (9.0 * M_FeO - 2.0 * M_SiO2 - 10.0 * M_m)) + M_MgSiO3 * (
                                                M_FeO * (9.0 * M_MgO - 2.0 * M_SiO2 - 10.0 * M_m) + M_FeSiO3 * (
                                                9.0 * M_FeO + 9.0 * M_MgO + 4.0 * M_SiO2 - 25.0 * M_m)))) + M_c * (
                                                                                              M_MgO * (
                                                                                              M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                              M_FeO * (
                                                                                              -M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                              M_FeO * (M_MgO * (
                                                                                              -M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                              M_FeO * (
                                                                                              -M_SiO2 + M_m) + M_MgO * (
                                                                                              -M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (
                                                                                              M_FeO * M_MgO * (
                                                                                              -M_FeSiO3 - M_SiO2) + M_MgSiO3 * (
                                                                                              -M_FeO * M_MgO + M_FeSiO3 * (
                                                                                              -M_FeO - M_MgO + M_m))) + M_Si * (
                                                                                              M_MgO * (M_FeO * (
                                                                                              -M_SiO2 + M_m) + M_FeSiO3 * (
                                                                                                       -4.0 * M_FeO + M_SiO2 + 3.0 * M_m)) + M_MgSiO3 * (
                                                                                              M_FeO * (
                                                                                              -4.0 * M_MgO + M_SiO2 + 3.0 * M_m) + M_FeSiO3 * (
                                                                                              -4.0 * M_FeO - 4.0 * M_MgO - M_SiO2 + 9.0 * M_m)) + M_O * (
                                                                                              M_MgO * (
                                                                                              -M_FeO - M_FeSiO3) + M_MgSiO3 * (
                                                                                              -M_FeO - M_FeSiO3)))))) + M_FeSiO3 * dKFeSiO3_KFeSiO3 * (
                       M_Fe * M_O * (-4.0 * M_MgO * M_SiO2 * M_m + M_MgSiO3 * (
                       M_FeO * (4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                       4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_Si * (M_Fe * (
                       -M_MgO * M_SiO2 * M_m + M_MgSiO3 * (
                       M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
                       M_MgO * (-2.0 * M_SiO2 - 10.0 * M_m) + M_MgSiO3 * (
                       9.0 * M_FeO + 9.0 * M_MgO + 4.0 * M_SiO2 - 25.0 * M_m))) + M_FeO * M_O * (M_MgO * (
                       -3.0 * M_SiO2 - 6.0 * M_m) + M_MgSiO3 * (6.0 * M_SiO2 - 15.0 * M_m))) + M_c * (M_Fe * (
                       M_MgO * M_SiO2 * M_m + M_MgSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_O * (
                       -M_FeO - M_MgO + M_m) + M_SiO2 * M_m)) + M_FeO * M_O * (
                                                                                                      M_MgO * M_SiO2 + M_MgSiO3 * M_m) + M_Si * (
                                                                                                      M_Fe * (M_MgO * (
                                                                                                      M_SiO2 + 3.0 * M_m) + M_MgSiO3 * (
                                                                                                              -4.0 * M_FeO - 4.0 * M_MgO - M_SiO2 + 9.0 * M_m) + M_O * (
                                                                                                              -M_MgO - M_MgSiO3)) + M_FeO * (
                                                                                                      M_MgO * (
                                                                                                      2.0 * M_SiO2 + 2.0 * M_m) + M_MgSiO3 * (
                                                                                                      -2.0 * M_SiO2 + 6.0 * M_m))))) + M_MgSiO3 * dKMgSiO3_KMgSiO3 * (
                       M_Fe * M_O * (4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                       M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                       -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_Si * (M_Fe * (
                       M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                       M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                       M_FeO * (-M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                       -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                       -3.0 * M_SiO2 - 6.0 * M_m) + 9.0 * M_SiO2 * M_m)) + M_O * (M_FeO * (
                       M_MgO * (-3.0 * M_SiO2 - 6.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (M_FeO * (
                       -9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (
                                                                                               -3.0 * M_SiO2 - 6.0 * M_m) + 9.0 * M_SiO2 * M_m))) + M_c * (
                       M_Fe * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                       M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
                               M_FeSiO3 * (M_FeO + M_MgO - M_m) + M_SiO2 * (M_FeO + M_MgO - M_m))) + M_O * (
                       M_FeO * M_SiO2 * (M_MgO - M_m) + M_FeSiO3 * (
                       M_FeO * (M_SiO2 - M_m) + M_SiO2 * (M_MgO - M_m))) + M_Si * (M_Fe * (
                       M_FeO * (M_SiO2 - M_m) + M_FeSiO3 * (4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                       2.0 * M_SiO2 + 2.0 * M_m) + M_O * (
                       M_FeO + M_FeSiO3 + M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeO * (M_MgO * (
                       2.0 * M_SiO2 + 2.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (M_FeO * (
                       4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (2.0 * M_SiO2 + 2.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_O * (
                                                                                   M_FeO * (M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                   M_SiO2 - M_m))))) + M_Si * (M_Fe * (
        M_FeSiO3 * (M_FeO * (M_SiO2 * (dM_MgO_er + dM_MgSiO3_er) + M_m * (-dM_MgO_er - dM_MgSiO3_er)) + M_MgO * (
        M_SiO2 * (1.0 * dM_FeSiO3_er + dM_MgO_er + 2.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er) + M_m * (
        1.0 * dM_FeO_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                    -dM_MgO_er - dM_MgSiO3_er)) + M_MgSiO3 * (M_FeO * (
        M_SiO2 * (-dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
        1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er + 1.0 * dM_SiO2_er)) + M_MgO * (
                                                              M_SiO2 * (-1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er) + M_m * (
                                                              1.0 * dM_FeO_er + 1.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                              1.0 * dM_FeO_er + 1.0 * dM_FeSiO3_er)) + M_O * (
        M_FeO * (M_SiO2 * (dM_MgO_er + dM_MgSiO3_er) + M_m * (-4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er)) + M_FeSiO3 * (
        M_FeO * (9.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er) + M_MgO * (
        6.0 * dM_FeO_er + 15.0 * dM_FeSiO3_er + 15.0 * dM_MgO_er + 24.0 * dM_MgSiO3_er + 9.0 * dM_SiO2_er) + M_SiO2 * (
        4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er) + M_m * (-25.0 * dM_MgO_er - 25.0 * dM_MgSiO3_er)) + M_MgO * (M_SiO2 * (
        2 * dM_FeO_er + 5.0 * dM_FeSiO3_er + 3.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 3.0 * dM_SiO2_er) + M_m * (
                                                                                                            4.0 * dM_FeO_er + 10.0 * dM_FeSiO3_er + 6.0 * dM_MgSiO3_er + 6.0 * dM_SiO2_er)) + M_MgSiO3 * (
        M_FeO * (
        -6.0 * dM_FeO_er - 15.0 * dM_FeSiO3_er + 3.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er) + M_MgO * (
        9.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er) + M_SiO2 * (
        -4.0 * dM_FeO_er - 10.0 * dM_FeSiO3_er - 6.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er) + M_m * (
        10.0 * dM_FeO_er + 25.0 * dM_FeSiO3_er - 15.0 * dM_MgO_er + 15.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
        -9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er)) + M_SiO2 * M_m * (
        M_FeO * (-dM_MgO_er - dM_MgSiO3_er) + M_MgO * (1.0 * dM_FeO_er + 1.0 * dM_FeSiO3_er))) + M_O * (M_FeO * (
        M_MgO * (M_SiO2 * (
        3.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er + 3.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 3.0 * dM_SiO2_er) + M_m * (
                 6.0 * dM_FeSiO3_er + 6.0 * dM_MgSiO3_er + 6.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
        -9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er)) + M_FeSiO3 * (M_FeO * (
        M_MgO * (9.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er + 9.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er) + M_SiO2 * (
        9.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er) + M_m * (-9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er)) + M_MgO * (M_SiO2 * (
        3.0 * dM_FeSiO3_er + 3.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 3.0 * dM_SiO2_er) + M_m * (
                                                                                                          -6.0 * dM_FeO_er + 6.0 * dM_MgSiO3_er + 6.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                              -9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er)) + M_MgSiO3 * (
                                                                                                        M_FeO * (
                                                                                                        M_MgO * (
                                                                                                        9.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er + 9.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                        -6.0 * dM_FeO_er - 12.0 * dM_FeSiO3_er - 6.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er) + M_m * (
                                                                                                        15.0 * dM_FeSiO3_er - 15.0 * dM_MgO_er + 15.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                                                        M_FeO * (
                                                                                                        9.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er + 9.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                                        9.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er + 9.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                        -6.0 * dM_FeSiO3_er - 6.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er) + M_m * (
                                                                                                        -15.0 * dM_FeO_er - 15.0 * dM_MgO_er + 15.0 * dM_SiO2_er)))) + dKSiO2_KSiO2 * (
                                                                                                               M_O * (
                                                                                                               M_Fe * (
                                                                                                               M_MgO * (
                                                                                                               M_FeSiO3 * (
                                                                                                               4.0 * M_SiO2 - 10.0 * M_m) - 6.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                               M_FeO * (
                                                                                                               2.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                                                                                                               6.0 * M_SiO2 - 6.0 * M_m) - 6.0 * M_SiO2 * M_m)) + M_MgO * (
                                                                                                               -6.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                               M_FeO * (
                                                                                                               6.0 * M_SiO2 - 6.0 * M_m) - 6.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                               M_FeO * (
                                                                                                               M_MgO * (
                                                                                                               6.0 * M_SiO2 - 6.0 * M_m) - 6.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                               M_FeO * (
                                                                                                               6.0 * M_SiO2 - 6.0 * M_m) + M_MgO * (
                                                                                                               6.0 * M_SiO2 - 6.0 * M_m) - 6.0 * M_SiO2 * M_m))) + M_c * (
                                                                                                               M_Fe * (
                                                                                                               M_MgO * (
                                                                                                               M_FeSiO3 * (
                                                                                                               -M_SiO2 + 3.0 * M_m) + 2.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                               M_FeO * (
                                                                                                               -M_SiO2 - M_m) + M_MgO * (
                                                                                                               -2.0 * M_SiO2 + 2.0 * M_m) + 2.0 * M_SiO2 * M_m) + M_O * (
                                                                                                               M_MgO * (
                                                                                                               -M_FeSiO3 - M_SiO2) + M_MgSiO3 * (
                                                                                                               M_FeO - M_m))) + M_MgO * (
                                                                                                               2.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                               M_FeO * (
                                                                                                               -2.0 * M_SiO2 + 2.0 * M_m) + 2.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                               M_FeO * (
                                                                                                               M_MgO * (
                                                                                                               -2.0 * M_SiO2 + 2.0 * M_m) + 2.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                               M_FeO * (
                                                                                                               -2.0 * M_SiO2 + 2.0 * M_m) + M_MgO * (
                                                                                                               -2.0 * M_SiO2 + 2.0 * M_m) + 2.0 * M_SiO2 * M_m)) + M_O * (
                                                                                                               M_MgO * M_SiO2 * (
                                                                                                               -M_FeO - M_FeSiO3) + M_MgSiO3 * M_m * (
                                                                                                               -M_FeO - M_FeSiO3))))) + M_c * (
                       M_Fe * (M_FeSiO3 * (
                       M_FeO * (M_SiO2 * (-dM_MgO_er - dM_MgSiO3_er) + M_m * (dM_MgO_er + dM_MgSiO3_er)) + M_MgO * (
                       M_SiO2 * (-dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                       -dM_FeO_er + 1.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                       dM_MgO_er + dM_MgSiO3_er)) + M_MgSiO3 * (M_FeO * (
                       M_SiO2 * (dM_FeO_er + 2.0 * dM_FeSiO3_er + 1.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er) + M_m * (
                       -dM_FeSiO3_er + dM_MgO_er - dM_SiO2_er)) + M_MgO * (
                                                                M_SiO2 * (dM_FeO_er + 1.0 * dM_FeSiO3_er) + M_m * (
                                                                -dM_FeO_er - 1.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                -dM_FeO_er - 1.0 * dM_FeSiO3_er)) + M_O * (M_FeSiO3 * (
                       M_FeO * (-dM_MgO_er - dM_MgSiO3_er) + M_MgO * (
                       -dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                       dM_MgO_er + dM_MgSiO3_er)) + M_MgSiO3 * (M_FeO * (
                       dM_FeSiO3_er - dM_MgO_er + dM_SiO2_er) + M_MgO * (-dM_MgO_er - dM_MgSiO3_er) + M_m * (
                                                                -dM_FeSiO3_er + dM_MgO_er - dM_SiO2_er)) + M_SiO2 * (
                                                                                                           M_FeO * (
                                                                                                           -dM_MgO_er - dM_MgSiO3_er) + M_MgO * (
                                                                                                           -dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                                                           dM_MgO_er + dM_MgSiO3_er))) + M_SiO2 * M_m * (
                               M_FeO * (dM_MgO_er + dM_MgSiO3_er) + M_MgO * (
                               -dM_FeO_er - 1.0 * dM_FeSiO3_er))) + M_O * (M_FeO * M_SiO2 * (
                       M_MgO * (-dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                       dM_MgO_er + dM_MgSiO3_er)) + M_FeSiO3 * (M_FeO * (
                       M_MgO * (-dM_FeO_er - dM_FeSiO3_er - dM_MgO_er - dM_MgSiO3_er) + M_SiO2 * (
                       -dM_MgO_er - dM_MgSiO3_er) + M_m * (dM_MgO_er + dM_MgSiO3_er)) + M_SiO2 * (M_MgO * (
                       -dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                                                  dM_MgO_er + dM_MgSiO3_er))) + M_MgSiO3 * (
                                                                           M_FeO * (M_MgO * (
                                                                           -dM_FeO_er - dM_FeSiO3_er - dM_MgO_er - dM_MgSiO3_er) + M_m * (
                                                                                    -dM_FeSiO3_er + dM_MgO_er - dM_SiO2_er)) + M_FeSiO3 * (
                                                                           M_FeO * (
                                                                           -dM_FeO_er - dM_FeSiO3_er - dM_MgO_er - dM_MgSiO3_er) + M_MgO * (
                                                                           -dM_FeO_er - dM_FeSiO3_er - dM_MgO_er - dM_MgSiO3_er) + M_m * (
                                                                           dM_FeO_er + dM_MgO_er - dM_SiO2_er)))) + M_Si * (
                       M_Fe * (
                       M_FeO * (M_SiO2 * (-dM_MgO_er - dM_MgSiO3_er) + M_m * (dM_MgO_er + dM_MgSiO3_er)) + M_FeSiO3 * (
                       M_FeO * (-4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_MgO * (
                       -2.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 10.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_SiO2 * (
                       -dM_MgO_er - dM_MgSiO3_er) + M_m * (9.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er)) + M_MgO * (M_SiO2 * (
                       -dM_FeO_er - 3.0 * dM_FeSiO3_er - 2 * dM_MgO_er - 4.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                                                                                                             -dM_FeO_er - 3.0 * dM_FeSiO3_er - 2.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er)) + M_MgSiO3 * (
                       M_FeO * (
                       2.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_MgO * (
                       -4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_SiO2 * (
                       dM_FeO_er + 3.0 * dM_FeSiO3_er + 2.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_m * (
                       -3.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er + 6.0 * dM_MgO_er - 6.0 * dM_SiO2_er)) + M_O * (
                       M_FeO * (-dM_MgO_er - dM_MgSiO3_er) + M_FeSiO3 * (-dM_MgO_er - dM_MgSiO3_er) + M_MgO * (
                       -dM_MgO_er - dM_MgSiO3_er) + M_MgSiO3 * (-dM_MgO_er - dM_MgSiO3_er) + M_SiO2 * (
                       -dM_MgO_er - dM_MgSiO3_er) + M_m * (dM_MgO_er + dM_MgSiO3_er)) + M_SiO2 * M_m * (
                       4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er)) + M_FeO * (M_MgO * (M_SiO2 * (
                       -2 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 2 * dM_MgO_er - 4.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                                                                                  -2.0 * dM_FeSiO3_er - 2.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                         4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er)) + M_FeSiO3 * (
                       M_FeO * (M_MgO * (
                       -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_SiO2 * (
                                -4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_m * (
                                4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er)) + M_MgO * (
                       M_SiO2 * (-2.0 * dM_FeSiO3_er - 2 * dM_MgO_er - 4.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                       2.0 * dM_FeO_er - 2.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                       4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er)) + M_MgSiO3 * (M_FeO * (M_MgO * (
                       -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                     2 * dM_FeO_er + 4.0 * dM_FeSiO3_er + 2.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_m * (
                                                                                     -6.0 * dM_FeSiO3_er + 6.0 * dM_MgO_er - 6.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                            M_FeO * (
                                                                            -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_MgO * (
                                                                            -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                            2.0 * dM_FeSiO3_er + 2.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_m * (
                                                                            6.0 * dM_FeO_er + 6.0 * dM_MgO_er - 6.0 * dM_SiO2_er))) + M_O * (
                       M_FeO * (M_MgO * (-dM_FeO_er - dM_FeSiO3_er - dM_MgO_er - dM_MgSiO3_er) + M_SiO2 * (
                       -dM_MgO_er - dM_MgSiO3_er) + M_m * (dM_MgO_er + dM_MgSiO3_er)) + M_FeSiO3 * (
                       M_MgO * (-dM_FeO_er - dM_FeSiO3_er - dM_MgO_er - dM_MgSiO3_er) + M_SiO2 * (
                       -dM_MgO_er - dM_MgSiO3_er) + M_m * (dM_MgO_er + dM_MgSiO3_er)) + M_MgSiO3 * (
                       M_FeO * (-dM_FeO_er - dM_FeSiO3_er - dM_MgO_er - dM_MgSiO3_er) + M_FeSiO3 * (
                       -dM_FeO_er - dM_FeSiO3_er - dM_MgO_er - dM_MgSiO3_er))))) + dKMgO_KMgO * (M_Fe * M_O * (M_MgO * (
        4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
        M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (M_FeO * (
        M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (M_FeO * (
        -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m))) + M_Si * (M_Fe * (
        M_MgO * (M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
        M_FeO * (M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
        M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (M_MgO * (
        M_FeO * (-M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
        -9.0 * M_FeO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
        -9.0 * M_MgO - M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                                                                                      -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                                                                                      -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m))) + M_O * (
                                                                                                           M_MgO * (
                                                                                                           9.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           M_MgO * (
                                                                                                           -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                           M_FeO * (
                                                                                                           -9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (
                                                                                                           -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)))) + M_c * (
                                                                                                 M_Fe * (M_MgO * (
                                                                                                 -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                 M_FeO * (
                                                                                                 M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                         M_FeO * (
                                                                                                         M_MgO * (
                                                                                                         M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                         M_FeO * (
                                                                                                         M_SiO2 - M_m) + M_MgO * (
                                                                                                         M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
                                                                                                         M_MgO * (
                                                                                                         M_FeSiO3 * (
                                                                                                         M_FeO - M_m) + M_SiO2 * (
                                                                                                         M_FeO - M_m)) + M_MgSiO3 * (
                                                                                                         M_FeO * (
                                                                                                         M_MgO + M_SiO2) + M_FeSiO3 * (
                                                                                                         M_FeO + M_MgO - M_m) + M_MgO * (
                                                                                                         M_SiO2 - M_m) - M_SiO2 * M_m))) + M_O * (
                                                                                                 M_MgO * (
                                                                                                 -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                 M_FeO * (
                                                                                                 M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                 M_FeO * (M_MgO * (
                                                                                                 M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                 M_FeO * (
                                                                                                 M_SiO2 - M_m) + M_MgO * (
                                                                                                 M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Si * (
                                                                                                 M_Fe * (M_MgO * (
                                                                                                 M_FeO * (
                                                                                                 M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                 4.0 * M_FeO + M_SiO2 - 9.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                         M_FeO * (
                                                                                                         4.0 * M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                         4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                                         4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_O * (
                                                                                                         M_MgO * (
                                                                                                         M_FeO + M_FeSiO3 + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                         M_FeO + M_FeSiO3 + M_SiO2 - M_m))) + M_MgO * (
                                                                                                 -4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                 M_FeO * (
                                                                                                 4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                 M_FeO * (M_MgO * (
                                                                                                 4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                 M_FeO * (
                                                                                                 4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                                 4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_O * (
                                                                                                 M_MgO * (M_FeO * (
                                                                                                 M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                          M_SiO2 - M_m)) + M_MgSiO3 * (
                                                                                                 M_FeO * (
                                                                                                 M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                 M_SiO2 - M_m))))))) / (
               M_O * (M_Fe * (M_MgO * (4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                              M_FeO * (M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                              M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                              -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m))) + M_Mg * (M_Fe * (M_FeSiO3 * (
               M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
               -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO)) + M_MgO * (
                                                                                           4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (M_MgO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)))) + M_Si * (
               M_Fe * (
               M_MgO * (M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (M_MgO * (
               M_FeO * (-M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
               -9.0 * M_FeO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               -9.0 * M_MgO - M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                                                                                             -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                                                                                             -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m))) + M_Mg * (
               M_Fe * (M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                       M_FeO * (-M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                       -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                       -M_SiO2 + 4.0 * M_m) + M_MgSiO3 * (
                       -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                       M_FeO + M_MgO)) + M_MgO * (
               M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (
               M_FeO * (M_MgO * (-M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-9.0 * M_MgO - 9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (
               -M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (-9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_FeSiO3 * (
               -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m)))) + M_O * (M_MgO * (
               9.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)) + M_MgSiO3 * (M_FeO * (
               M_MgO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (M_FeO * (
               -9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)))) + M_c * (
               M_Fe * (
               M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
               M_MgO * (M_FeSiO3 * (M_FeO - M_m) + M_SiO2 * (M_FeO - M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO + M_SiO2) + M_FeSiO3 * (M_FeO + M_MgO - M_m) + M_MgO * (
               M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Mg * (M_Fe * (
               M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
               M_FeSiO3 * (M_FeO + M_MgO - M_m) + M_MgSiO3 * (M_FeO + M_MgO - M_m) + M_SiO2 * (
               M_FeO + M_MgO - M_m)) + M_SiO2 * M_m * (-M_FeO - M_MgO)) + M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                         M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                         M_FeO * (M_SiO2 - M_m) + M_MgO * (
                                                         M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
                                                         M_FeO * M_SiO2 * (M_MgO - M_m) + M_FeSiO3 * (
                                                         M_FeO * (M_MgO + M_SiO2 - M_m) + M_SiO2 * (
                                                         M_MgO - M_m)) + M_MgSiO3 * (
                                                         M_FeO * (M_MgO - M_m) + M_FeSiO3 * (
                                                         M_FeO + M_MgO - M_m)))) + M_O * (
               M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Si * (M_Fe * (M_MgO * (
               M_FeO * (M_SiO2 - M_m) + M_FeSiO3 * (
               4.0 * M_FeO + M_SiO2 - 9.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               4.0 * M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                     4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_O * (
                                                                                                   M_MgO * (
                                                                                                   M_FeO + M_FeSiO3 + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                   M_FeO + M_FeSiO3 + M_SiO2 - M_m))) + M_Mg * (
                                                                                           M_Fe * (M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                                   M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_O * (
                                                                                                   M_FeO + M_FeSiO3 + M_MgO + M_MgSiO3 + M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeO * (
                                                                                           M_MgO * (
                                                                                           M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_MgO + 4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                           M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_FeSiO3 * (
                                                                                           4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m)) + M_O * (
                                                                                           M_FeO * (
                                                                                           M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                           M_MgO + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                           M_FeO + M_FeSiO3))) + M_MgO * (
                                                                                           -4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (M_MgO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_O * (
                                                                                           M_MgO * (M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                    M_SiO2 - M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                           M_SiO2 - M_m))))))

    def dM_m_dTc(self, Moles, dKs, dMi_b):
        dM_Mg_er, dM_Si_er, dM_Fe_er, dM_O_er, dM_c_er, dM_MgO_er, dM_SiO2_er, dM_FeO_er, dM_MgSiO3_er, dM_FeSiO3_er, dM_m_er = self.unwrap_Moles(
            dMi_b)
        M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = Moles
        dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs
        return M_m * (M_Fe * (M_O * (M_MgO * (M_FeO * M_SiO2 * (
        -4.0 * dM_FeO_er - 8.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 8.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_FeSiO3 * (
                                              M_FeO * (
                                              -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_SiO2 * (
                                              -4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 8.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er))) + M_MgSiO3 * (
                                     M_FeO * (M_MgO * (
                                     -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_SiO2 * (
                                              -4.0 * dM_FeO_er - 8.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                     M_FeO * (
                                     -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_MgO * (
                                     -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_SiO2 * (
                                     -4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)))) + dKFeO_KFeO * (
                              M_Mg * M_O * M_SiO2 * (-4.0 * M_FeO * M_MgSiO3 + 4.0 * M_FeSiO3 * M_MgO) + M_Si * (
                              M_Mg * (M_O * (M_FeO * (-15.0 * M_MgSiO3 - 3.0 * M_SiO2) + M_FeSiO3 * (
                              -9.0 * M_FeO + 6.0 * M_MgO + 6.0 * M_SiO2)) + M_SiO2 * (
                                      -M_FeO * M_MgSiO3 + M_FeSiO3 * M_MgO)) + M_O * (
                              M_MgO * (-3.0 * M_FeO * M_SiO2 + M_FeSiO3 * (-9.0 * M_FeO + 6.0 * M_SiO2)) + M_MgSiO3 * (
                              M_FeO * (-9.0 * M_MgO - 3.0 * M_SiO2) + M_FeSiO3 * (
                              -9.0 * M_FeO - 9.0 * M_MgO + 6.0 * M_SiO2)))) + M_c * (M_Mg * (M_FeO * (
                              M_MgSiO3 * M_SiO2 + M_O * (
                              M_FeSiO3 + M_MgSiO3 + M_SiO2)) - M_FeSiO3 * M_MgO * M_SiO2) + M_O * (M_FeO * M_MgO * (
                              M_FeSiO3 + M_SiO2) + M_MgSiO3 * (M_FeO * (M_MgO + M_SiO2) + M_FeSiO3 * (
                              M_FeO + M_MgO))) + M_Si * (M_Mg * (M_FeO * (6.0 * M_MgSiO3 + 2.0 * M_SiO2) + M_FeSiO3 * (
                              4.0 * M_FeO - 2.0 * M_MgO - 2.0 * M_SiO2) + M_O * (M_FeO + M_FeSiO3)) + M_MgO * (
                                                         2.0 * M_FeO * M_SiO2 + M_FeSiO3 * (
                                                         4.0 * M_FeO - 2.0 * M_SiO2)) + M_MgSiO3 * (
                                                         M_FeO * (4.0 * M_MgO + 2.0 * M_SiO2) + M_FeSiO3 * (
                                                         4.0 * M_FeO + 4.0 * M_MgO - 2.0 * M_SiO2)) + M_O * (
                                                         M_MgO * (M_FeO + M_FeSiO3) + M_MgSiO3 * (
                                                         M_FeO + M_FeSiO3)))))) + M_FeSiO3 * dKFeSiO3_KFeSiO3 * (
                      M_O * M_SiO2 * (M_Fe * M_FeO * (4.0 * M_MgO + 4.0 * M_MgSiO3) + M_Mg * (
                      M_Fe * (4.0 * M_FeO + 4.0 * M_MgO) + M_FeO * (4.0 * M_MgO + 4.0 * M_MgSiO3))) + M_Si * (M_Fe * (
                      M_FeO * M_SiO2 * (M_MgO + M_MgSiO3) + M_O * (M_MgO * (6.0 * M_FeO + 6.0 * M_SiO2) + M_MgSiO3 * (
                      6.0 * M_FeO - 9.0 * M_MgO + 6.0 * M_SiO2))) + M_FeO * M_O * M_SiO2 * (
                                                                                                              9.0 * M_MgO + 9.0 * M_MgSiO3) + M_Mg * (
                                                                                                              M_Fe * (
                                                                                                              M_O * (
                                                                                                              6.0 * M_FeO + 6.0 * M_MgO + 6.0 * M_SiO2) + M_SiO2 * (
                                                                                                              M_FeO + M_MgO)) + M_FeO * (
                                                                                                              M_O * (
                                                                                                              6.0 * M_MgO + 15.0 * M_MgSiO3 + 9.0 * M_SiO2) + M_SiO2 * (
                                                                                                              M_MgO + M_MgSiO3)))) + M_c * (
                      M_Fe * (
                      -M_FeO * M_MgO * M_SiO2 + M_MgSiO3 * (-M_FeO * M_SiO2 + M_MgO * M_O)) + M_FeO * M_O * M_SiO2 * (
                      -M_MgO - M_MgSiO3) + M_Mg * (M_Fe * M_SiO2 * (-M_FeO - M_MgO) + M_FeO * (
                      M_O * (-M_MgSiO3 - M_SiO2) + M_SiO2 * (-M_MgO - M_MgSiO3))) + M_Si * (M_Fe * (
                      M_MgO * (-2.0 * M_FeO - 2.0 * M_SiO2) + M_MgSiO3 * (
                      -2.0 * M_FeO + 4.0 * M_MgO - 2.0 * M_SiO2) + M_O * (M_MgO + M_MgSiO3)) + M_FeO * M_SiO2 * (
                                                                                            -4.0 * M_MgO - 4.0 * M_MgSiO3) + M_Mg * (
                                                                                            M_Fe * (
                                                                                            -2.0 * M_FeO - 2.0 * M_MgO + M_O - 2.0 * M_SiO2) + M_FeO * (
                                                                                            -2.0 * M_MgO - 6.0 * M_MgSiO3 - 4.0 * M_SiO2))))) + M_Mg * (
                      M_O * (M_Fe * (M_FeSiO3 * (
                      M_FeO * (-4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_MgO * (
                      -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_SiO2 * (
                      -4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)) + M_MgSiO3 * (M_FeO * (
                      -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                                  -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                  -4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)) + M_SiO2 * (
                                     M_FeO * (
                                     -4.0 * dM_FeO_er - 8.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 8.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_MgO * (
                                     -4.0 * dM_FeO_er - 8.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 8.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er))) + M_MgO * (
                             M_FeO * M_SiO2 * (
                             -4.0 * dM_FeO_er - 8.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 8.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er) + M_FeSiO3 * (
                             M_FeO * (
                             -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_SiO2 * (
                             -4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 8.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er))) + M_MgSiO3 * (
                             M_FeO * (M_MgO * (
                             -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_SiO2 * (
                                      -4.0 * dM_FeO_er - 8.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)) + M_FeSiO3 * (
                             M_FeO * (
                             -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_MgO * (
                             -4.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er) + M_SiO2 * (
                             -4.0 * dM_FeSiO3_er - 4.0 * dM_MgSiO3_er - 4.0 * dM_SiO2_er)))) + dKMgO_KMgO * (
                      M_Fe * M_O * M_SiO2 * (4.0 * M_FeO * M_MgSiO3 - 4.0 * M_FeSiO3 * M_MgO) + M_Si * (M_Fe * (M_O * (
                      M_MgO * (-15.0 * M_FeSiO3 - 3.0 * M_SiO2) + M_MgSiO3 * (
                      6.0 * M_FeO - 9.0 * M_MgO + 6.0 * M_SiO2)) + M_SiO2 * (
                                                                                                                M_FeO * M_MgSiO3 - M_FeSiO3 * M_MgO)) + M_O * (
                                                                                                        M_MgO * (
                                                                                                        -3.0 * M_FeO * M_SiO2 + M_FeSiO3 * (
                                                                                                        -9.0 * M_FeO - 3.0 * M_SiO2)) + M_MgSiO3 * (
                                                                                                        M_FeO * (
                                                                                                        -9.0 * M_MgO + 6.0 * M_SiO2) + M_FeSiO3 * (
                                                                                                        -9.0 * M_FeO - 9.0 * M_MgO + 6.0 * M_SiO2)))) + M_c * (
                      M_Fe * (-M_FeO * M_MgSiO3 * M_SiO2 + M_MgO * (
                      M_FeSiO3 * M_SiO2 + M_O * (M_FeSiO3 + M_MgSiO3 + M_SiO2))) + M_O * (
                      M_MgO * (M_FeO * M_SiO2 + M_FeSiO3 * (M_FeO + M_SiO2)) + M_MgSiO3 * (
                      M_FeO * M_MgO + M_FeSiO3 * (M_FeO + M_MgO))) + M_Si * (M_Fe * (
                      M_MgO * (6.0 * M_FeSiO3 + 2.0 * M_SiO2) + M_MgSiO3 * (
                      -2.0 * M_FeO + 4.0 * M_MgO - 2.0 * M_SiO2) + M_O * (M_MgO + M_MgSiO3)) + M_MgO * (
                                                                             2.0 * M_FeO * M_SiO2 + M_FeSiO3 * (
                                                                             4.0 * M_FeO + 2.0 * M_SiO2)) + M_MgSiO3 * (
                                                                             M_FeO * (
                                                                             4.0 * M_MgO - 2.0 * M_SiO2) + M_FeSiO3 * (
                                                                             4.0 * M_FeO + 4.0 * M_MgO - 2.0 * M_SiO2)) + M_O * (
                                                                             M_MgO * (M_FeO + M_FeSiO3) + M_MgSiO3 * (
                                                                             M_FeO + M_FeSiO3)))))) + M_MgSiO3 * dKMgSiO3_KMgSiO3 * (
                      M_O * M_SiO2 * (M_Fe * M_MgO * (4.0 * M_FeO + 4.0 * M_FeSiO3) + M_Mg * (
                      M_Fe * (4.0 * M_FeO + 4.0 * M_MgO) + M_MgO * (4.0 * M_FeO + 4.0 * M_FeSiO3))) + M_Si * (M_Mg * (
                      M_Fe * (
                      M_O * (6.0 * M_FeO + 6.0 * M_MgO + 6.0 * M_SiO2) + M_SiO2 * (M_FeO + M_MgO)) + M_MgO * M_SiO2 * (
                      M_FeO + M_FeSiO3) + M_O * (M_FeO * (6.0 * M_MgO + 6.0 * M_SiO2) + M_FeSiO3 * (
                      -9.0 * M_FeO + 6.0 * M_MgO + 6.0 * M_SiO2))) + M_MgO * (M_Fe * (
                      M_O * (6.0 * M_FeO + 15.0 * M_FeSiO3 + 9.0 * M_SiO2) + M_SiO2 * (
                      M_FeO + M_FeSiO3)) + M_O * M_SiO2 * (9.0 * M_FeO + 9.0 * M_FeSiO3))) + M_c * (M_Mg * (
                      M_FeSiO3 * (M_FeO * M_O - M_MgO * M_SiO2) + M_SiO2 * (
                      M_Fe * (-M_FeO - M_MgO) - M_FeO * M_MgO)) + M_MgO * (M_Fe * (
                      M_O * (-M_FeSiO3 - M_SiO2) + M_SiO2 * (-M_FeO - M_FeSiO3)) + M_O * M_SiO2 * (
                                                                           -M_FeO - M_FeSiO3)) + M_Si * (M_Mg * (
                      M_Fe * (-2.0 * M_FeO - 2.0 * M_MgO + M_O - 2.0 * M_SiO2) + M_FeO * (
                      -2.0 * M_MgO - 2.0 * M_SiO2) + M_FeSiO3 * (4.0 * M_FeO - 2.0 * M_MgO - 2.0 * M_SiO2) + M_O * (
                      M_FeO + M_FeSiO3)) + M_MgO * (M_Fe * (-2.0 * M_FeO - 6.0 * M_FeSiO3 - 4.0 * M_SiO2) + M_SiO2 * (
                      -4.0 * M_FeO - 4.0 * M_FeSiO3))))) + M_Si * (M_Fe * (M_MgO * (
        M_FeO * M_SiO2 * (-dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_FeSiO3 * (
        M_FeO * (-1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_SiO2 * (
        -1.0 * dM_FeSiO3_er - dM_MgO_er - 2.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er))) + M_MgSiO3 * (M_FeO * (
        M_MgO * (-1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_SiO2 * (
        -dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er)) + M_FeSiO3 * (M_FeO * (
        -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                  -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                  -1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er))) + M_O * (
                                                                           M_MgO * (M_FeO * (
                                                                           -4.0 * dM_FeO_er - 10.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 10.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                                    -10.0 * dM_FeO_er - 25.0 * dM_FeSiO3_er - 25.0 * dM_MgO_er - 40.0 * dM_MgSiO3_er - 15.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                    -6.0 * dM_FeO_er - 15.0 * dM_FeSiO3_er - 9.0 * dM_MgO_er - 18.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er)) + M_MgSiO3 * (
                                                                           M_FeO * (
                                                                           -4.0 * dM_FeO_er - 10.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                           -10.0 * dM_FeO_er - 25.0 * dM_FeSiO3_er - 10.0 * dM_MgO_er - 25.0 * dM_MgSiO3_er - 15.0 * dM_SiO2_er) + M_MgO * (
                                                                           -9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                           -6.0 * dM_FeO_er - 15.0 * dM_FeSiO3_er - 9.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er)))) + M_Mg * (
                                                                   M_Fe * (M_FeSiO3 * (M_FeO * (
                                                                   -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                       -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                       -1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er)) + M_MgSiO3 * (
                                                                           M_FeO * (
                                                                           -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_MgO * (
                                                                           -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                           -1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er)) + M_O * (
                                                                           M_FeO * (
                                                                           -4.0 * dM_FeO_er - 10.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 10.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                           -10.0 * dM_FeO_er - 25.0 * dM_FeSiO3_er - 10.0 * dM_MgO_er - 25.0 * dM_MgSiO3_er - 15.0 * dM_SiO2_er) + M_MgO * (
                                                                           -4.0 * dM_FeO_er - 10.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 10.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er) + M_MgSiO3 * (
                                                                           -10.0 * dM_FeO_er - 25.0 * dM_FeSiO3_er - 10.0 * dM_MgO_er - 25.0 * dM_MgSiO3_er - 15.0 * dM_SiO2_er) + M_SiO2 * (
                                                                           -6.0 * dM_FeO_er - 15.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 15.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er)) + M_SiO2 * (
                                                                           M_FeO * (
                                                                           -dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_MgO * (
                                                                           -dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er))) + M_MgO * (
                                                                   M_FeO * M_SiO2 * (
                                                                   -dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_FeSiO3 * (
                                                                   M_FeO * (
                                                                   -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                   -dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er))) + M_MgSiO3 * (
                                                                   M_FeO * (M_MgO * (
                                                                   -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                            -dM_FeO_er - 2.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                   M_FeO * (
                                                                   -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_MgO * (
                                                                   -1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                   -1.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er))) + M_O * (
                                                                   M_FeO * (M_MgO * (
                                                                   -4.0 * dM_FeO_er - 10.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 10.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er) + M_SiO2 * (
                                                                            -9.0 * dM_FeO_er - 18.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 15.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                   M_FeO * (
                                                                   -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er) + M_MgO * (
                                                                   2.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 4.0 * dM_MgO_er - 10.0 * dM_MgSiO3_er - 6.0 * dM_SiO2_er) + M_SiO2 * (
                                                                   -9.0 * dM_FeSiO3_er - 6.0 * dM_MgO_er - 15.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er)) + M_MgSiO3 * (
                                                                   M_FeO * (
                                                                   -25.0 * dM_FeO_er - 40.0 * dM_FeSiO3_er - 10.0 * dM_MgO_er - 25.0 * dM_MgSiO3_er - 15.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                   -10.0 * dM_FeO_er - 25.0 * dM_FeSiO3_er - 10.0 * dM_MgO_er - 25.0 * dM_MgSiO3_er - 15.0 * dM_SiO2_er)))) + M_O * (
                                                                   M_MgO * (M_FeO * M_SiO2 * (
                                                                   -9.0 * dM_FeO_er - 18.0 * dM_FeSiO3_er - 9.0 * dM_MgO_er - 18.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                            M_FeO * (
                                                                            -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er - 9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                            -9.0 * dM_FeSiO3_er - 9.0 * dM_MgO_er - 18.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er))) + M_MgSiO3 * (
                                                                   M_FeO * (M_MgO * (
                                                                   -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er - 9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                            -9.0 * dM_FeO_er - 18.0 * dM_FeSiO3_er - 9.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                   M_FeO * (
                                                                   -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er - 9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er) + M_MgO * (
                                                                   -9.0 * dM_FeO_er - 9.0 * dM_FeSiO3_er - 9.0 * dM_MgO_er - 9.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                   -9.0 * dM_FeSiO3_er - 9.0 * dM_MgSiO3_er - 9.0 * dM_SiO2_er)))) + dKSiO2_KSiO2 * (
                                                                   M_O * (M_Fe * (M_MgO * (
                                                                   2.0 * M_FeO * M_SiO2 + M_FeSiO3 * (
                                                                   6.0 * M_FeO - 4.0 * M_SiO2)) + M_MgSiO3 * (M_FeO * (
                                                                   6.0 * M_MgO + 2.0 * M_SiO2) + M_FeSiO3 * (
                                                                                                              6.0 * M_FeO + 6.0 * M_MgO - 4.0 * M_SiO2))) + M_Mg * (
                                                                          M_Fe * (M_FeSiO3 * (
                                                                          6.0 * M_FeO + 6.0 * M_MgO - 4.0 * M_SiO2) + M_MgSiO3 * (
                                                                                  6.0 * M_FeO + 6.0 * M_MgO - 4.0 * M_SiO2) + M_SiO2 * (
                                                                                  2.0 * M_FeO + 2.0 * M_MgO)) + M_MgO * (
                                                                          2.0 * M_FeO * M_SiO2 + M_FeSiO3 * (
                                                                          6.0 * M_FeO + 2.0 * M_SiO2)) + M_MgSiO3 * (
                                                                          M_FeO * (
                                                                          6.0 * M_MgO - 4.0 * M_SiO2) + M_FeSiO3 * (
                                                                          6.0 * M_FeO + 6.0 * M_MgO - 4.0 * M_SiO2)))) + M_c * (
                                                                   M_Fe * (M_MgO * (-M_FeO * M_SiO2 + M_FeSiO3 * (
                                                                   -2.0 * M_FeO + M_SiO2)) + M_MgSiO3 * (M_FeO * (
                                                                   -2.0 * M_MgO - M_SiO2) + M_FeSiO3 * (
                                                                                                         -2.0 * M_FeO - 2.0 * M_MgO + M_SiO2)) + M_O * (
                                                                           M_MgO * (M_FeSiO3 + M_SiO2) + M_MgSiO3 * (
                                                                           M_FeSiO3 + M_SiO2))) + M_Mg * (M_Fe * (
                                                                   M_FeSiO3 * (
                                                                   -2.0 * M_FeO - 2.0 * M_MgO + M_SiO2) + M_MgSiO3 * (
                                                                   -2.0 * M_FeO - 2.0 * M_MgO + M_SiO2) + M_O * (
                                                                   M_FeSiO3 + M_MgSiO3 + M_SiO2) + M_SiO2 * (
                                                                   -M_FeO - M_MgO)) + M_MgO * (
                                                                                                          -M_FeO * M_SiO2 + M_FeSiO3 * (
                                                                                                          -2.0 * M_FeO - M_SiO2)) + M_MgSiO3 * (
                                                                                                          M_FeO * (
                                                                                                          -2.0 * M_MgO + M_SiO2) + M_FeSiO3 * (
                                                                                                          -2.0 * M_FeO - 2.0 * M_MgO + M_SiO2)) + M_O * (
                                                                                                          M_MgSiO3 * (
                                                                                                          M_FeO + M_FeSiO3) + M_SiO2 * (
                                                                                                          M_FeO + M_FeSiO3))) + M_O * M_SiO2 * (
                                                                   M_MgO * (M_FeO + M_FeSiO3) + M_MgSiO3 * (
                                                                   M_FeO + M_FeSiO3))))) + M_c * (M_Fe * (M_MgO * (
        M_FeO * M_SiO2 * (
        dM_FeO_er + 2.0 * dM_FeSiO3_er + dM_MgO_er + 2.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er) + M_FeSiO3 * (
        M_FeO * (dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_SiO2 * (
        dM_FeSiO3_er + dM_MgO_er + 2 * dM_MgSiO3_er + dM_SiO2_er))) + M_MgSiO3 * (M_FeO * (
        M_MgO * (dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_SiO2 * (
        dM_FeO_er + 2.0 * dM_FeSiO3_er + 1.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er)) + M_FeSiO3 * (M_FeO * (
        dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                               dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                               dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er))) + M_O * (
                                                                                                          M_MgO * (
                                                                                                          M_FeSiO3 * (
                                                                                                          dM_FeSiO3_er + dM_MgO_er + 2 * dM_MgSiO3_er + dM_SiO2_er) + M_SiO2 * (
                                                                                                          dM_FeSiO3_er + dM_MgO_er + 2 * dM_MgSiO3_er + dM_SiO2_er)) + M_MgSiO3 * (
                                                                                                          M_FeSiO3 * (
                                                                                                          dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_MgO * (
                                                                                                          dM_MgO_er + dM_MgSiO3_er) + M_SiO2 * (
                                                                                                          dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)))) + M_Mg * (
                                                                                                  M_Fe * (M_FeSiO3 * (
                                                                                                  M_FeO * (
                                                                                                  dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                                  dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                  dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_MgSiO3 * (
                                                                                                          M_FeO * (
                                                                                                          dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                                          dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                          dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_O * (
                                                                                                          M_FeSiO3 * (
                                                                                                          dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_MgSiO3 * (
                                                                                                          dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_SiO2 * (
                                                                                                          dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_SiO2 * (
                                                                                                          M_FeO * (
                                                                                                          dM_FeO_er + 2.0 * dM_FeSiO3_er + dM_MgO_er + 2.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er) + M_MgO * (
                                                                                                          dM_FeO_er + 2.0 * dM_FeSiO3_er + dM_MgO_er + 2.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er))) + M_MgO * (
                                                                                                  M_FeO * M_SiO2 * (
                                                                                                  dM_FeO_er + 2.0 * dM_FeSiO3_er + dM_MgO_er + 2.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                                                  M_FeO * (
                                                                                                  dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                  1.0 * dM_FeSiO3_er + dM_MgO_er + 2.0 * dM_MgSiO3_er + 1.0 * dM_SiO2_er))) + M_MgSiO3 * (
                                                                                                  M_FeO * (M_MgO * (
                                                                                                  dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                           dM_FeO_er + 2 * dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_FeSiO3 * (
                                                                                                  M_FeO * (
                                                                                                  dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                                  dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                  dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er))) + M_O * (
                                                                                                  M_FeO * M_SiO2 * (
                                                                                                  dM_FeO_er + 2 * dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_FeSiO3 * (
                                                                                                  M_FeO * (
                                                                                                  dM_FeO_er + dM_FeSiO3_er) + M_SiO2 * (
                                                                                                  dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_MgSiO3 * (
                                                                                                  M_FeO * (
                                                                                                  dM_FeO_er + 2 * dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_FeSiO3 * (
                                                                                                  dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)))) + M_O * (
                                                                                                  M_MgO * (
                                                                                                  M_FeO * M_SiO2 * (
                                                                                                  dM_FeO_er + 2 * dM_FeSiO3_er + dM_MgO_er + 2 * dM_MgSiO3_er + dM_SiO2_er) + M_FeSiO3 * (
                                                                                                  M_FeO * (
                                                                                                  dM_FeO_er + dM_FeSiO3_er + dM_MgO_er + dM_MgSiO3_er) + M_SiO2 * (
                                                                                                  dM_FeSiO3_er + dM_MgO_er + 2 * dM_MgSiO3_er + dM_SiO2_er))) + M_MgSiO3 * (
                                                                                                  M_FeO * (M_MgO * (
                                                                                                  dM_FeO_er + dM_FeSiO3_er + dM_MgO_er + dM_MgSiO3_er) + M_SiO2 * (
                                                                                                           dM_FeO_er + 2 * dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)) + M_FeSiO3 * (
                                                                                                  M_FeO * (
                                                                                                  dM_FeO_er + dM_FeSiO3_er + dM_MgO_er + dM_MgSiO3_er) + M_MgO * (
                                                                                                  dM_FeO_er + dM_FeSiO3_er + dM_MgO_er + dM_MgSiO3_er) + M_SiO2 * (
                                                                                                  dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er)))) + M_Si * (
                                                                                                  M_Fe * (M_MgO * (
                                                                                                  M_FeO * (
                                                                                                  dM_FeO_er + 3.0 * dM_FeSiO3_er + dM_MgO_er + 3.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                                                  3.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er + 9.0 * dM_MgO_er + 15.0 * dM_MgSiO3_er + 6.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                                  2.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 8.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er)) + M_MgSiO3 * (
                                                                                                          M_FeO * (
                                                                                                          dM_FeO_er + 3.0 * dM_FeSiO3_er - 1.0 * dM_MgO_er + 1.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                                                          3.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er + 3.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er + 6.0 * dM_SiO2_er) + M_MgO * (
                                                                                                          4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                          2.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er)) + M_O * (
                                                                                                          M_MgO * (
                                                                                                          dM_MgO_er + dM_MgSiO3_er) + M_MgSiO3 * (
                                                                                                          dM_MgO_er + dM_MgSiO3_er))) + M_Mg * (
                                                                                                  M_Fe * (M_FeO * (
                                                                                                  dM_FeO_er + 3.0 * dM_FeSiO3_er + dM_MgO_er + 3.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                                                          3.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er + 3.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er + 6.0 * dM_SiO2_er) + M_MgO * (
                                                                                                          dM_FeO_er + 3.0 * dM_FeSiO3_er + dM_MgO_er + 3.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_MgSiO3 * (
                                                                                                          3.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er + 3.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er + 6.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                                          2.0 * dM_FeO_er + 6.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er)) + M_FeO * (
                                                                                                  M_MgO * (
                                                                                                  dM_FeO_er + 3.0 * dM_FeSiO3_er + dM_MgO_er + 3.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                                  4.0 * dM_FeO_er + 8.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                                                  M_FeO * (
                                                                                                  4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er) + M_MgO * (
                                                                                                  -1.0 * dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 3.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                                  4.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 6.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er)) + M_MgSiO3 * (
                                                                                                  M_FeO * (
                                                                                                  9.0 * dM_FeO_er + 15.0 * dM_FeSiO3_er + 3.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er + 6.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                                                  3.0 * dM_FeO_er + 9.0 * dM_FeSiO3_er + 3.0 * dM_MgO_er + 9.0 * dM_MgSiO3_er + 6.0 * dM_SiO2_er)) + M_O * (
                                                                                                  M_FeO * (
                                                                                                  dM_FeO_er + dM_FeSiO3_er) + M_FeSiO3 * (
                                                                                                  dM_FeO_er + dM_FeSiO3_er))) + M_MgO * (
                                                                                                  M_FeO * M_SiO2 * (
                                                                                                  4.0 * dM_FeO_er + 8.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 8.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                                                  M_FeO * (
                                                                                                  4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                  4.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 8.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er))) + M_MgSiO3 * (
                                                                                                  M_FeO * (M_MgO * (
                                                                                                  4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                           4.0 * dM_FeO_er + 8.0 * dM_FeSiO3_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                                                  M_FeO * (
                                                                                                  4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                                  4.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er + 4.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                  4.0 * dM_FeSiO3_er + 4.0 * dM_MgSiO3_er + 4.0 * dM_SiO2_er))) + M_O * (
                                                                                                  M_MgO * (M_FeO * (
                                                                                                  dM_FeO_er + dM_FeSiO3_er + dM_MgO_er + dM_MgSiO3_er) + M_FeSiO3 * (
                                                                                                           dM_FeO_er + dM_FeSiO3_er + dM_MgO_er + dM_MgSiO3_er)) + M_MgSiO3 * (
                                                                                                  M_FeO * (
                                                                                                  dM_FeO_er + dM_FeSiO3_er + dM_MgO_er + dM_MgSiO3_er) + M_FeSiO3 * (
                                                                                                  dM_FeO_er + dM_FeSiO3_er + dM_MgO_er + dM_MgSiO3_er)))))) / (
               M_O * (M_Fe * (M_MgO * (4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                              M_FeO * (M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                              M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                              -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m))) + M_Mg * (M_Fe * (M_FeSiO3 * (
               M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
               -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO)) + M_MgO * (
                                                                                           4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (M_MgO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)))) + M_Si * (
               M_Fe * (
               M_MgO * (M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (M_MgO * (
               M_FeO * (-M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
               -9.0 * M_FeO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               -9.0 * M_MgO - M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                                                                                             -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                                                                                             -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m))) + M_Mg * (
               M_Fe * (M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                       M_FeO * (-M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                       -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                       -M_SiO2 + 4.0 * M_m) + M_MgSiO3 * (
                       -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                       M_FeO + M_MgO)) + M_MgO * (
               M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (
               M_FeO * (M_MgO * (-M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-9.0 * M_MgO - 9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (
               -M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (-9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_FeSiO3 * (
               -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m)))) + M_O * (M_MgO * (
               9.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)) + M_MgSiO3 * (M_FeO * (
               M_MgO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (M_FeO * (
               -9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)))) + M_c * (
               M_Fe * (
               M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
               M_MgO * (M_FeSiO3 * (M_FeO - M_m) + M_SiO2 * (M_FeO - M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO + M_SiO2) + M_FeSiO3 * (M_FeO + M_MgO - M_m) + M_MgO * (
               M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Mg * (M_Fe * (
               M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
               M_FeSiO3 * (M_FeO + M_MgO - M_m) + M_MgSiO3 * (M_FeO + M_MgO - M_m) + M_SiO2 * (
               M_FeO + M_MgO - M_m)) + M_SiO2 * M_m * (-M_FeO - M_MgO)) + M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                         M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                         M_FeO * (M_SiO2 - M_m) + M_MgO * (
                                                         M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
                                                         M_FeO * M_SiO2 * (M_MgO - M_m) + M_FeSiO3 * (
                                                         M_FeO * (M_MgO + M_SiO2 - M_m) + M_SiO2 * (
                                                         M_MgO - M_m)) + M_MgSiO3 * (
                                                         M_FeO * (M_MgO - M_m) + M_FeSiO3 * (
                                                         M_FeO + M_MgO - M_m)))) + M_O * (
               M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Si * (M_Fe * (M_MgO * (
               M_FeO * (M_SiO2 - M_m) + M_FeSiO3 * (
               4.0 * M_FeO + M_SiO2 - 9.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               4.0 * M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                     4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_O * (
                                                                                                   M_MgO * (
                                                                                                   M_FeO + M_FeSiO3 + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                   M_FeO + M_FeSiO3 + M_SiO2 - M_m))) + M_Mg * (
                                                                                           M_Fe * (M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                                   M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_O * (
                                                                                                   M_FeO + M_FeSiO3 + M_MgO + M_MgSiO3 + M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeO * (
                                                                                           M_MgO * (
                                                                                           M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_MgO + 4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                           M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_FeSiO3 * (
                                                                                           4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m)) + M_O * (
                                                                                           M_FeO * (
                                                                                           M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                           M_MgO + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                           M_FeO + M_FeSiO3))) + M_MgO * (
                                                                                           -4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (M_MgO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_O * (
                                                                                           M_MgO * (M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                    M_SiO2 - M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                           M_SiO2 - M_m))))))

    def dM_c_dTc(self, Moles, dKs, dMi_b):
        dM_Mg_er, dM_Si_er, dM_Fe_er, dM_O_er, dM_c_er, dM_MgO_er, dM_SiO2_er, dM_FeO_er, dM_MgSiO3_er, dM_FeSiO3_er, dM_m_er = self.unwrap_Moles(
            dMi_b)
        M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = Moles
        dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs
        return M_c * (M_Fe * (M_O * (M_MgO * (M_FeSiO3 * (
        M_FeO * (-2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_m * (
        2.0 * dM_FeO_er - 2.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er)) + M_SiO2 * (M_FeO * (
        -2.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                                                                              2.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er))) + M_MgSiO3 * (
                                     M_FeO * (M_MgO * (
                                     -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_SiO2 * (
                                              -2.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 2.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                     M_FeO * (
                                     -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_MgO * (
                                     -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_m * (
                                     2.0 * dM_FeO_er + 2.0 * dM_MgO_er - 2.0 * dM_SiO2_er)) + M_MgO * (
                                     M_SiO2 * (-2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er) + M_m * (
                                     2.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                     2.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er))) + dKFeO_KFeO * (M_O * (M_Mg * (
        M_FeO * M_m * (-2.0 * M_MgSiO3 - 2.0 * M_SiO2) + M_FeSiO3 * (
        M_FeO * (2.0 * M_SiO2 - 2.0 * M_m) + M_SiO2 * (2.0 * M_MgO - 2.0 * M_m))) + M_MgO * (
                                                                                                    -2.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                    M_FeO * (
                                                                                                    2.0 * M_SiO2 - 2.0 * M_m) - 2.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                    M_FeO * (M_MgO * (
                                                                                                    2.0 * M_SiO2 - 2.0 * M_m) - 2.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                    M_FeO * (
                                                                                                    2.0 * M_SiO2 - 2.0 * M_m) + M_MgO * (
                                                                                                    2.0 * M_SiO2 - 2.0 * M_m) - 2.0 * M_SiO2 * M_m))) + M_Si * (
                                                                                             M_Mg * (M_FeO * (
                                                                                             M_MgSiO3 * (
                                                                                             M_SiO2 - 3.0 * M_m) - 2.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                     M_FeO * (
                                                                                                     2.0 * M_SiO2 - 2.0 * M_m) + M_MgO * (
                                                                                                     M_SiO2 + M_m) - 2.0 * M_SiO2 * M_m) + M_O * (
                                                                                                     M_FeO * (
                                                                                                     -3.0 * M_MgSiO3 - M_SiO2 - 2.0 * M_m) + M_FeSiO3 * (
                                                                                                     3.0 * M_MgO + 2.0 * M_SiO2 - 5.0 * M_m))) + M_MgO * (
                                                                                             -2.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                             M_FeO * (
                                                                                             2.0 * M_SiO2 - 2.0 * M_m) - 2.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                             M_FeO * (M_MgO * (
                                                                                             2.0 * M_SiO2 - 2.0 * M_m) - 2.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                             M_FeO * (
                                                                                             2.0 * M_SiO2 - 2.0 * M_m) + M_MgO * (
                                                                                             2.0 * M_SiO2 - 2.0 * M_m) - 2.0 * M_SiO2 * M_m)) + M_O * (
                                                                                             M_MgO * (M_FeO * (
                                                                                             -M_SiO2 - 2.0 * M_m) + M_FeSiO3 * (
                                                                                                      2.0 * M_SiO2 - 5.0 * M_m)) + M_MgSiO3 * (
                                                                                             M_FeO * (
                                                                                             -M_SiO2 - 2.0 * M_m) + M_FeSiO3 * (
                                                                                             2.0 * M_SiO2 - 5.0 * M_m)))))) + M_FeSiO3 * dKFeSiO3_KFeSiO3 * (
                      M_O * (M_Fe * (M_MgO * M_SiO2 * (2.0 * M_FeO - 2.0 * M_m) + M_MgSiO3 * (
                      M_MgO * (2.0 * M_SiO2 - 2.0 * M_m) + M_SiO2 * (2.0 * M_FeO - 2.0 * M_m))) + M_Mg * (
                             M_Fe * M_SiO2 * (2.0 * M_FeO + 2.0 * M_MgO - 2.0 * M_m) + M_FeO * (
                             2.0 * M_MgO * M_SiO2 + 2.0 * M_MgSiO3 * M_m))) + M_Si * (M_Fe * (
                      M_MgO * (M_FeO * (M_SiO2 + M_m) - 2.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                      M_FeO * (M_SiO2 + M_m) + M_MgO * (2.0 * M_SiO2 - 2.0 * M_m) - 2.0 * M_SiO2 * M_m) + M_O * (
                      M_MgO * (3.0 * M_FeO + 2.0 * M_SiO2 - 5.0 * M_m) + M_MgSiO3 * (
                      3.0 * M_FeO + 2.0 * M_SiO2 - 5.0 * M_m))) + M_FeO * M_O * (M_MgO * (
                      3.0 * M_SiO2 - 3.0 * M_m) + M_MgSiO3 * (3.0 * M_SiO2 - 3.0 * M_m)) + M_Mg * (M_Fe * (
                      M_FeO * (M_SiO2 + M_m) + M_MgO * (M_SiO2 + M_m) + M_O * (
                      3.0 * M_FeO + 3.0 * M_MgO + 2.0 * M_SiO2 - 5.0 * M_m) - 2.0 * M_SiO2 * M_m) + M_FeO * (M_MgO * (
                      M_SiO2 + M_m) + M_MgSiO3 * (-M_SiO2 + 3.0 * M_m) + M_O * (
                                                                                                             3.0 * M_MgO + 3.0 * M_MgSiO3 + 3.0 * M_SiO2 - 3.0 * M_m))))) + M_Mg * (
                      M_O * (M_Fe * (M_FeSiO3 * (
                      M_FeO * (-2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_MgO * (
                      -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_m * (
                      2.0 * dM_FeO_er + 2.0 * dM_MgO_er - 2.0 * dM_SiO2_er)) + M_MgSiO3 * (M_FeO * (
                      -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                           -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_m * (
                                                                                           2.0 * dM_FeO_er + 2.0 * dM_MgO_er - 2.0 * dM_SiO2_er)) + M_SiO2 * (
                                     M_FeO * (
                                     -2.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_MgO * (
                                     -2.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                                     2.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er))) + M_FeO * M_SiO2 * (
                             M_MgO * (
                             -2.0 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                             2.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er)) + M_FeSiO3 * (M_FeO * (M_MgO * (
                      -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                           -2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_m * (
                                                                                           2.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er)) + M_SiO2 * (
                                                                                  M_MgO * (
                                                                                  -2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 4.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                                                                                  2.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er))) + M_MgSiO3 * (
                             M_FeO * (M_MgO * (
                             -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_m * (
                                      -2.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er - 2.0 * dM_SiO2_er)) + M_FeSiO3 * (M_FeO * (
                             -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                                               -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_m * (
                                                                                                               2.0 * dM_FeO_er + 2.0 * dM_MgO_er - 2.0 * dM_SiO2_er)))) + dKMgO_KMgO * (
                      M_O * (M_Fe * (M_MgO * M_m * (-2.0 * M_FeSiO3 - 2.0 * M_SiO2) + M_MgSiO3 * (
                      M_MgO * (2.0 * M_SiO2 - 2.0 * M_m) + M_SiO2 * (2.0 * M_FeO - 2.0 * M_m))) + M_MgO * (
                             -2.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                             M_FeO * (2.0 * M_SiO2 - 2.0 * M_m) - 2.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                             M_FeO * (M_MgO * (2.0 * M_SiO2 - 2.0 * M_m) - 2.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                             M_FeO * (2.0 * M_SiO2 - 2.0 * M_m) + M_MgO * (
                             2.0 * M_SiO2 - 2.0 * M_m) - 2.0 * M_SiO2 * M_m))) + M_Si * (M_Fe * (
                      M_MgO * (M_FeSiO3 * (M_SiO2 - 3.0 * M_m) - 2.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                      M_FeO * (M_SiO2 + M_m) + M_MgO * (2.0 * M_SiO2 - 2.0 * M_m) - 2.0 * M_SiO2 * M_m) + M_O * (
                      M_MgO * (-3.0 * M_FeSiO3 - M_SiO2 - 2.0 * M_m) + M_MgSiO3 * (
                      3.0 * M_FeO + 2.0 * M_SiO2 - 5.0 * M_m))) + M_MgO * (-2.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                      M_FeO * (2.0 * M_SiO2 - 2.0 * M_m) - 2.0 * M_SiO2 * M_m)) + M_MgSiO3 * (M_FeO * (
                      M_MgO * (2.0 * M_SiO2 - 2.0 * M_m) - 2.0 * M_SiO2 * M_m) + M_FeSiO3 * (M_FeO * (
                      2.0 * M_SiO2 - 2.0 * M_m) + M_MgO * (2.0 * M_SiO2 - 2.0 * M_m) - 2.0 * M_SiO2 * M_m)) + M_O * (
                                                                                         M_MgO * (M_FeO * (
                                                                                         -M_SiO2 - 2.0 * M_m) + M_FeSiO3 * (
                                                                                                  -M_SiO2 - 2.0 * M_m)) + M_MgSiO3 * (
                                                                                         M_FeO * (
                                                                                         2.0 * M_SiO2 - 5.0 * M_m) + M_FeSiO3 * (
                                                                                         2.0 * M_SiO2 - 5.0 * M_m)))))) + M_MgSiO3 * dKMgSiO3_KMgSiO3 * (
                      M_O * (M_Fe * M_MgO * (2.0 * M_FeO * M_SiO2 + 2.0 * M_FeSiO3 * M_m) + M_Mg * (
                      M_FeSiO3 * (M_FeO * (2.0 * M_SiO2 - 2.0 * M_m) + M_SiO2 * (2.0 * M_MgO - 2.0 * M_m)) + M_SiO2 * (
                      M_Fe * (2.0 * M_FeO + 2.0 * M_MgO - 2.0 * M_m) + M_FeO * (2.0 * M_MgO - 2.0 * M_m)))) + M_Si * (
                      M_Mg * (M_Fe * (M_FeO * (M_SiO2 + M_m) + M_MgO * (M_SiO2 + M_m) + M_O * (
                      3.0 * M_FeO + 3.0 * M_MgO + 2.0 * M_SiO2 - 5.0 * M_m) - 2.0 * M_SiO2 * M_m) + M_FeO * (
                              M_MgO * (M_SiO2 + M_m) - 2.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                              M_FeO * (2.0 * M_SiO2 - 2.0 * M_m) + M_MgO * (
                              M_SiO2 + M_m) - 2.0 * M_SiO2 * M_m) + M_O * (
                              M_FeO * (3.0 * M_MgO + 2.0 * M_SiO2 - 5.0 * M_m) + M_FeSiO3 * (
                              3.0 * M_MgO + 2.0 * M_SiO2 - 5.0 * M_m))) + M_MgO * (M_Fe * (
                      M_FeO * (M_SiO2 + M_m) + M_FeSiO3 * (-M_SiO2 + 3.0 * M_m) + M_O * (
                      3.0 * M_FeO + 3.0 * M_FeSiO3 + 3.0 * M_SiO2 - 3.0 * M_m)) + M_O * (M_FeO * (
                      3.0 * M_SiO2 - 3.0 * M_m) + M_FeSiO3 * (3.0 * M_SiO2 - 3.0 * M_m))))) + M_Si * (M_Fe * (M_MgO * (
        M_FeO * (M_SiO2 * (-dM_FeO_er - 2.0 * dM_FeSiO3_er - dM_MgO_er - 2.0 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
        -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er)) + M_FeSiO3 * (
        M_FeO * (-2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_SiO2 * (
        dM_FeSiO3_er + dM_MgO_er + 2.0 * dM_MgSiO3_er + dM_SiO2_er) + M_m * (
        3.0 * dM_FeO_er - 3.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
        2.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er)) + M_MgSiO3 * (M_FeO * (
        M_MgO * (-2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_SiO2 * (
        -dM_FeO_er - 2.0 * dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
        -dM_FeSiO3_er + dM_MgO_er - dM_SiO2_er)) + M_FeSiO3 * (M_FeO * (
        -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_MgO * (
                                                               -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                               dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                                                               3.0 * dM_FeO_er + 3.0 * dM_MgO_er - 3.0 * dM_SiO2_er)) + M_MgO * (
                                                             M_SiO2 * (-2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er) + M_m * (
                                                             2.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                             2.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er)) + M_O * (M_MgO * (
        M_FeO * (
        -2.0 * dM_FeO_er - 5.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 5.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_FeSiO3 * (
        -2.0 * dM_FeO_er - 5.0 * dM_FeSiO3_er - 5.0 * dM_MgO_er - 8.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_SiO2 * (
        -2.0 * dM_FeO_er - 5.0 * dM_FeSiO3_er - 3.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_m * (
        2.0 * dM_FeO_er + 5.0 * dM_FeSiO3_er + 3.0 * dM_MgSiO3_er + 3.0 * dM_SiO2_er)) + M_MgSiO3 * (M_FeO * (
        -2.0 * dM_FeO_er - 5.0 * dM_FeSiO3_er + dM_MgO_er - 2.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                                                     -2.0 * dM_FeO_er - 5.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 5.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                                     -2.0 * dM_FeO_er - 5.0 * dM_FeSiO3_er - 3.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_m * (
                                                                                                     2.0 * dM_FeO_er + 5.0 * dM_FeSiO3_er - 3.0 * dM_MgO_er + 3.0 * dM_SiO2_er)))) + M_Mg * (
                                                                                                      M_Fe * (M_FeO * (
                                                                                                      M_SiO2 * (
                                                                                                      -dM_FeO_er - 2.0 * dM_FeSiO3_er - dM_MgO_er - 2.0 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                                                      -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er)) + M_FeSiO3 * (
                                                                                                              M_FeO * (
                                                                                                              -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                                              -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                              dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                                                                                                              3.0 * dM_FeO_er + 3.0 * dM_MgO_er - 3.0 * dM_SiO2_er)) + M_MgO * (
                                                                                                              M_SiO2 * (
                                                                                                              -dM_FeO_er - 2.0 * dM_FeSiO3_er - dM_MgO_er - 2.0 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                                                              -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er)) + M_MgSiO3 * (
                                                                                                              M_FeO * (
                                                                                                              -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                                              -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                              dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                                                                                                              3.0 * dM_FeO_er + 3.0 * dM_MgO_er - 3.0 * dM_SiO2_er)) + M_O * (
                                                                                                              M_FeO * (
                                                                                                              -2.0 * dM_FeO_er - 5.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 5.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                                                              -2.0 * dM_FeO_er - 5.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 5.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_MgO * (
                                                                                                              -2.0 * dM_FeO_er - 5.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 5.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_MgSiO3 * (
                                                                                                              -2.0 * dM_FeO_er - 5.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 5.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                                              -2.0 * dM_FeO_er - 5.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 5.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_m * (
                                                                                                              2.0 * dM_FeO_er + 5.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 5.0 * dM_MgSiO3_er + 3.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                                              2.0 * dM_FeO_er + 2.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er)) + M_FeO * (
                                                                                                      M_MgO * (
                                                                                                      M_SiO2 * (
                                                                                                      -dM_FeO_er - 2.0 * dM_FeSiO3_er - dM_MgO_er - 2.0 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                                                      -dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                                      2.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er)) + M_FeSiO3 * (
                                                                                                      M_FeO * (M_MgO * (
                                                                                                      -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                               -2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_m * (
                                                                                                               2.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er)) + M_MgO * (
                                                                                                      M_SiO2 * (
                                                                                                      -dM_FeSiO3_er - dM_MgO_er - 2.0 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                                                      dM_FeO_er - dM_MgSiO3_er - dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                                                      2.0 * dM_MgO_er + 2.0 * dM_MgSiO3_er)) + M_MgSiO3 * (
                                                                                                      M_FeO * (M_MgO * (
                                                                                                      -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                               dM_FeO_er + 2.0 * dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                                                                                                               -3.0 * dM_FeSiO3_er + 3.0 * dM_MgO_er - 3.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                                                      M_FeO * (
                                                                                                      -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                                      -2.0 * dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 2.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                      dM_FeSiO3_er + dM_MgSiO3_er + dM_SiO2_er) + M_m * (
                                                                                                      3.0 * dM_FeO_er + 3.0 * dM_MgO_er - 3.0 * dM_SiO2_er))) + M_O * (
                                                                                                      M_FeO * (M_MgO * (
                                                                                                      -2.0 * dM_FeO_er - 5.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 5.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                                               -3.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 5.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_m * (
                                                                                                               3.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 5.0 * dM_MgSiO3_er + 3.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                                                      M_MgO * (
                                                                                                      dM_FeO_er - 2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 5.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                                      -3.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 5.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_m * (
                                                                                                      -3.0 * dM_FeO_er + 2.0 * dM_MgO_er + 5.0 * dM_MgSiO3_er + 3.0 * dM_SiO2_er)) + M_MgSiO3 * (
                                                                                                      M_FeO * (
                                                                                                      -5.0 * dM_FeO_er - 8.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 5.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                                                      -2.0 * dM_FeO_er - 5.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er - 5.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er)))) + M_O * (
                                                                                                      M_MgO * (M_FeO * (
                                                                                                      M_SiO2 * (
                                                                                                      -3.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 3.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_m * (
                                                                                                      3.0 * dM_FeSiO3_er + 3.0 * dM_MgSiO3_er + 3.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                                                               M_SiO2 * (
                                                                                                               -3.0 * dM_FeSiO3_er - 3.0 * dM_MgO_er - 6.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_m * (
                                                                                                               -3.0 * dM_FeO_er + 3.0 * dM_MgSiO3_er + 3.0 * dM_SiO2_er))) + M_MgSiO3 * (
                                                                                                      M_FeO * (
                                                                                                      M_SiO2 * (
                                                                                                      -3.0 * dM_FeO_er - 6.0 * dM_FeSiO3_er - 3.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_m * (
                                                                                                      3.0 * dM_FeSiO3_er - 3.0 * dM_MgO_er + 3.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                                                      M_SiO2 * (
                                                                                                      -3.0 * dM_FeSiO3_er - 3.0 * dM_MgSiO3_er - 3.0 * dM_SiO2_er) + M_m * (
                                                                                                      -3.0 * dM_FeO_er - 3.0 * dM_MgO_er + 3.0 * dM_SiO2_er)))) + dKSiO2_KSiO2 * (
                                                                                                      M_Fe * (M_MgO * (
                                                                                                      M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                      M_FeO * (
                                                                                                      -M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                              M_FeO * (
                                                                                                              M_MgO * (
                                                                                                              -M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                              M_FeO * (
                                                                                                              -M_SiO2 + M_m) + M_MgO * (
                                                                                                              -M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (
                                                                                                              M_MgO * (
                                                                                                              M_FeSiO3 * (
                                                                                                              3.0 * M_FeO - 5.0 * M_m) + M_SiO2 * (
                                                                                                              M_FeO - 3.0 * M_m)) + M_MgSiO3 * (
                                                                                                              M_FeO * (
                                                                                                              3.0 * M_MgO + M_SiO2) + M_FeSiO3 * (
                                                                                                              3.0 * M_FeO + 3.0 * M_MgO - 5.0 * M_m) + M_MgO * (
                                                                                                              3.0 * M_SiO2 - 3.0 * M_m) - 3.0 * M_SiO2 * M_m))) + M_Mg * (
                                                                                                      M_Fe * (
                                                                                                      M_FeSiO3 * (
                                                                                                      M_FeO * (
                                                                                                      -M_SiO2 + M_m) + M_MgO * (
                                                                                                      -M_SiO2 + M_m) + M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                      M_FeO * (
                                                                                                      -M_SiO2 + M_m) + M_MgO * (
                                                                                                      -M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                                                                                                      M_FeSiO3 * (
                                                                                                      3.0 * M_FeO + 3.0 * M_MgO - 5.0 * M_m) + M_MgSiO3 * (
                                                                                                      3.0 * M_FeO + 3.0 * M_MgO - 5.0 * M_m) + M_SiO2 * (
                                                                                                      M_FeO + M_MgO - 3.0 * M_m)) + M_SiO2 * M_m * (
                                                                                                      M_FeO + M_MgO)) + M_MgO * (
                                                                                                      M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                      M_FeO * (
                                                                                                      -M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                      M_FeO * (M_MgO * (
                                                                                                      -M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                      M_FeO * (
                                                                                                      -M_SiO2 + M_m) + M_MgO * (
                                                                                                      -M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (
                                                                                                      M_FeO * M_SiO2 * (
                                                                                                      M_MgO - 3.0 * M_m) + M_FeSiO3 * (
                                                                                                      M_FeO * (
                                                                                                      3.0 * M_MgO + 3.0 * M_SiO2 - 3.0 * M_m) + M_SiO2 * (
                                                                                                      M_MgO - 3.0 * M_m)) + M_MgSiO3 * (
                                                                                                      M_FeO * (
                                                                                                      3.0 * M_MgO - 5.0 * M_m) + M_FeSiO3 * (
                                                                                                      3.0 * M_FeO + 3.0 * M_MgO - 5.0 * M_m)))) + M_O * (
                                                                                                      M_MgO * (
                                                                                                      -3.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                      M_FeO * (
                                                                                                      3.0 * M_SiO2 - 3.0 * M_m) - 3.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                      M_FeO * (M_MgO * (
                                                                                                      3.0 * M_SiO2 - 3.0 * M_m) - 3.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                      M_FeO * (
                                                                                                      3.0 * M_SiO2 - 3.0 * M_m) + M_MgO * (
                                                                                                      3.0 * M_SiO2 - 3.0 * M_m) - 3.0 * M_SiO2 * M_m)))))) / (
               M_O * (M_Fe * (M_MgO * (4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                              M_FeO * (M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                              M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                              -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m))) + M_Mg * (M_Fe * (M_FeSiO3 * (
               M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
               -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO)) + M_MgO * (
                                                                                           4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (M_MgO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)))) + M_Si * (
               M_Fe * (
               M_MgO * (M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (M_MgO * (
               M_FeO * (-M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
               -9.0 * M_FeO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               -9.0 * M_MgO - M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                                                                                             -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                                                                                             -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m))) + M_Mg * (
               M_Fe * (M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                       M_FeO * (-M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                       -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                       -M_SiO2 + 4.0 * M_m) + M_MgSiO3 * (
                       -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                       M_FeO + M_MgO)) + M_MgO * (
               M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (
               M_FeO * (M_MgO * (-M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-9.0 * M_MgO - 9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (
               -M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (-9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_FeSiO3 * (
               -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m)))) + M_O * (M_MgO * (
               9.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)) + M_MgSiO3 * (M_FeO * (
               M_MgO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (M_FeO * (
               -9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)))) + M_c * (
               M_Fe * (
               M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
               M_MgO * (M_FeSiO3 * (M_FeO - M_m) + M_SiO2 * (M_FeO - M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO + M_SiO2) + M_FeSiO3 * (M_FeO + M_MgO - M_m) + M_MgO * (
               M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Mg * (M_Fe * (
               M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
               M_FeSiO3 * (M_FeO + M_MgO - M_m) + M_MgSiO3 * (M_FeO + M_MgO - M_m) + M_SiO2 * (
               M_FeO + M_MgO - M_m)) + M_SiO2 * M_m * (-M_FeO - M_MgO)) + M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                         M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                         M_FeO * (M_SiO2 - M_m) + M_MgO * (
                                                         M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
                                                         M_FeO * M_SiO2 * (M_MgO - M_m) + M_FeSiO3 * (
                                                         M_FeO * (M_MgO + M_SiO2 - M_m) + M_SiO2 * (
                                                         M_MgO - M_m)) + M_MgSiO3 * (
                                                         M_FeO * (M_MgO - M_m) + M_FeSiO3 * (
                                                         M_FeO + M_MgO - M_m)))) + M_O * (
               M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Si * (M_Fe * (M_MgO * (
               M_FeO * (M_SiO2 - M_m) + M_FeSiO3 * (
               4.0 * M_FeO + M_SiO2 - 9.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               4.0 * M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                     4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_O * (
                                                                                                   M_MgO * (
                                                                                                   M_FeO + M_FeSiO3 + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                   M_FeO + M_FeSiO3 + M_SiO2 - M_m))) + M_Mg * (
                                                                                           M_Fe * (M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                                   M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_O * (
                                                                                                   M_FeO + M_FeSiO3 + M_MgO + M_MgSiO3 + M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeO * (
                                                                                           M_MgO * (
                                                                                           M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_MgO + 4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                           M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_FeSiO3 * (
                                                                                           4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m)) + M_O * (
                                                                                           M_FeO * (
                                                                                           M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                           M_MgO + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                           M_FeO + M_FeSiO3))) + M_MgO * (
                                                                                           -4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (M_MgO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_O * (
                                                                                           M_MgO * (M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                    M_SiO2 - M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                           M_SiO2 - M_m))))))

    def dM_O_dTc(self, Moles, dKs, dMi_b):
        dM_Mg_er, dM_Si_er, dM_Fe_er, dM_O_er, dM_c_er, dM_MgO_er, dM_SiO2_er, dM_FeO_er, dM_MgSiO3_er, dM_FeSiO3_er, dM_m_er = self.unwrap_Moles(
            dMi_b)
        M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = Moles
        dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs
        return M_O * (M_Fe * dKFeO_KFeO * (M_Si * (M_Mg * (
        M_FeO * (M_MgSiO3 * (2.0 * M_SiO2 - 5.0 * M_m) - 3.0 * M_SiO2 * M_m) + M_FeSiO3 * (
        M_FeO * (3.0 * M_SiO2 - 3.0 * M_m) + M_MgO * (M_SiO2 + 2.0 * M_m) - 3.0 * M_SiO2 * M_m)) + M_MgO * (
                                                   -3.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (
                                                   3.0 * M_SiO2 - 3.0 * M_m) - 3.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                   M_FeO * (M_MgO * (
                                                   3.0 * M_SiO2 - 3.0 * M_m) - 3.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                   M_FeO * (3.0 * M_SiO2 - 3.0 * M_m) + M_MgO * (
                                                   3.0 * M_SiO2 - 3.0 * M_m) - 3.0 * M_SiO2 * M_m))) + M_c * (M_Mg * (
        M_FeO * M_m * (-M_MgSiO3 - M_SiO2) + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) + M_SiO2 * (M_MgO - M_m))) + M_MgO * (
                                                                                                              -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                              M_FeO * (
                                                                                                              M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                              M_FeO * (
                                                                                                              M_MgO * (
                                                                                                              M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                              M_FeO * (
                                                                                                              M_SiO2 - M_m) + M_MgO * (
                                                                                                              M_SiO2 - M_m) - M_SiO2 * M_m)) + M_Si * (
                                                                                                              M_Mg * (
                                                                                                              M_FeO * (
                                                                                                              -2.0 * M_MgSiO3 - M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                              2.0 * M_MgO + M_SiO2 - 3.0 * M_m)) + M_MgO * (
                                                                                                              M_FeO * (
                                                                                                              -M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                              M_SiO2 - 3.0 * M_m)) + M_MgSiO3 * (
                                                                                                              M_FeO * (
                                                                                                              -M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                              M_SiO2 - 3.0 * M_m))))) + M_FeSiO3 * dKFeSiO3_KFeSiO3 * (
                      M_Si * (M_Fe * (M_MgO * (M_FeO * (M_SiO2 + 2.0 * M_m) - 3.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                      M_FeO * (M_SiO2 + 2.0 * M_m) + M_MgO * (
                      3.0 * M_SiO2 - 3.0 * M_m) - 3.0 * M_SiO2 * M_m)) + M_Mg * (M_Fe * (
                      M_FeO * (M_SiO2 + 2.0 * M_m) + M_MgO * (M_SiO2 + 2.0 * M_m) - 3.0 * M_SiO2 * M_m) + M_FeO * (
                                                                                 M_MgO * (
                                                                                 M_SiO2 + 2.0 * M_m) + M_MgSiO3 * (
                                                                                 -2.0 * M_SiO2 + 5.0 * M_m)))) + M_c * (
                      M_Fe * (M_MgO * M_SiO2 * (M_FeO - M_m) + M_MgSiO3 * (
                      M_MgO * (M_SiO2 - M_m) + M_SiO2 * (M_FeO - M_m))) + M_Mg * (
                      M_Fe * M_SiO2 * (M_FeO + M_MgO - M_m) + M_FeO * (M_MgO * M_SiO2 + M_MgSiO3 * M_m)) + M_Si * (
                      M_Fe * (M_MgO * (2.0 * M_FeO + M_SiO2 - 3.0 * M_m) + M_MgSiO3 * (
                      2.0 * M_FeO + M_SiO2 - 3.0 * M_m)) + M_FeO * (
                      M_MgO * (2.0 * M_SiO2 - 2.0 * M_m) + M_MgSiO3 * (2.0 * M_SiO2 - 2.0 * M_m)) + M_Mg * (
                      M_Fe * (2.0 * M_FeO + 2.0 * M_MgO + M_SiO2 - 3.0 * M_m) + M_FeO * (
                      2.0 * M_MgO + 2.0 * M_MgSiO3 + 2.0 * M_SiO2 - 2.0 * M_m))))) + M_Mg * dKMgO_KMgO * (M_Si * (
        M_Fe * (M_MgO * (M_FeSiO3 * (2.0 * M_SiO2 - 5.0 * M_m) - 3.0 * M_SiO2 * M_m) + M_MgSiO3 * (
        M_FeO * (M_SiO2 + 2.0 * M_m) + M_MgO * (3.0 * M_SiO2 - 3.0 * M_m) - 3.0 * M_SiO2 * M_m)) + M_MgO * (
        -3.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
        M_FeO * (3.0 * M_SiO2 - 3.0 * M_m) - 3.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
        M_FeO * (M_MgO * (3.0 * M_SiO2 - 3.0 * M_m) - 3.0 * M_SiO2 * M_m) + M_FeSiO3 * (
        M_FeO * (3.0 * M_SiO2 - 3.0 * M_m) + M_MgO * (3.0 * M_SiO2 - 3.0 * M_m) - 3.0 * M_SiO2 * M_m))) + M_c * (
                                                                                                          M_Fe * (
                                                                                                          M_MgO * M_m * (
                                                                                                          -M_FeSiO3 - M_SiO2) + M_MgSiO3 * (
                                                                                                          M_MgO * (
                                                                                                          M_SiO2 - M_m) + M_SiO2 * (
                                                                                                          M_FeO - M_m))) + M_MgO * (
                                                                                                          -M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                          M_FeO * (
                                                                                                          M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                          M_FeO * (
                                                                                                          M_MgO * (
                                                                                                          M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                          M_FeO * (
                                                                                                          M_SiO2 - M_m) + M_MgO * (
                                                                                                          M_SiO2 - M_m) - M_SiO2 * M_m)) + M_Si * (
                                                                                                          M_Fe * (
                                                                                                          M_MgO * (
                                                                                                          -2.0 * M_FeSiO3 - M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                          2.0 * M_FeO + M_SiO2 - 3.0 * M_m)) + M_MgO * (
                                                                                                          M_FeO * (
                                                                                                          -M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                          -M_SiO2 - M_m)) + M_MgSiO3 * (
                                                                                                          M_FeO * (
                                                                                                          M_SiO2 - 3.0 * M_m) + M_FeSiO3 * (
                                                                                                          M_SiO2 - 3.0 * M_m))))) + M_MgSiO3 * dKMgSiO3_KMgSiO3 * (
                      M_Si * (
                      M_Fe * M_MgO * (M_FeO * (M_SiO2 + 2.0 * M_m) + M_FeSiO3 * (-2.0 * M_SiO2 + 5.0 * M_m)) + M_Mg * (
                      M_Fe * (
                      M_FeO * (M_SiO2 + 2.0 * M_m) + M_MgO * (M_SiO2 + 2.0 * M_m) - 3.0 * M_SiO2 * M_m) + M_FeO * (
                      M_MgO * (M_SiO2 + 2.0 * M_m) - 3.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                      M_FeO * (3.0 * M_SiO2 - 3.0 * M_m) + M_MgO * (
                      M_SiO2 + 2.0 * M_m) - 3.0 * M_SiO2 * M_m))) + M_c * (
                      M_Fe * M_MgO * (M_FeO * M_SiO2 + M_FeSiO3 * M_m) + M_Mg * (
                      M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) + M_SiO2 * (M_MgO - M_m)) + M_SiO2 * (
                      M_Fe * (M_FeO + M_MgO - M_m) + M_FeO * (M_MgO - M_m))) + M_Si * (M_Mg * (
                      M_Fe * (2.0 * M_FeO + 2.0 * M_MgO + M_SiO2 - 3.0 * M_m) + M_FeO * (
                      2.0 * M_MgO + M_SiO2 - 3.0 * M_m) + M_FeSiO3 * (2.0 * M_MgO + M_SiO2 - 3.0 * M_m)) + M_MgO * (
                                                                                       M_Fe * (
                                                                                       2.0 * M_FeO + 2.0 * M_FeSiO3 + 2.0 * M_SiO2 - 2.0 * M_m) + M_FeO * (
                                                                                       2.0 * M_SiO2 - 2.0 * M_m) + M_FeSiO3 * (
                                                                                       2.0 * M_SiO2 - 2.0 * M_m))))) + M_Si * (
                      M_Fe * (M_MgO * (M_FeO * (
                      M_SiO2 * (-dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                      -2.0 * dM_FeSiO3_er - 2.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er)) + M_FeSiO3 * (M_FeO * (
                      -3.0 * dM_FeO_er - 3.0 * dM_FeSiO3_er - 3.0 * dM_MgO_er - 3.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                  2.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er + 4.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_m * (
                                                                                                  5.0 * dM_FeO_er - 5.0 * dM_MgSiO3_er - 5.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                       3.0 * dM_FeO_er + 3.0 * dM_FeSiO3_er)) + M_MgSiO3 * (M_FeO * (M_MgO * (
                      -3.0 * dM_FeO_er - 3.0 * dM_FeSiO3_er - 3.0 * dM_MgO_er - 3.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                     -dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                                                     -2.0 * dM_FeSiO3_er + 2.0 * dM_MgO_er - 2.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                                            M_FeO * (
                                                                                            -3.0 * dM_FeO_er - 3.0 * dM_FeSiO3_er - 3.0 * dM_MgO_er - 3.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                            -3.0 * dM_FeO_er - 3.0 * dM_FeSiO3_er - 3.0 * dM_MgO_er - 3.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                            2.0 * dM_FeSiO3_er + 2.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_m * (
                                                                                            5.0 * dM_FeO_er + 5.0 * dM_MgO_er - 5.0 * dM_SiO2_er)) + M_MgO * (
                                                                                            M_SiO2 * (
                                                                                            -3.0 * dM_FeO_er - 3.0 * dM_FeSiO3_er) + M_m * (
                                                                                            3.0 * dM_FeO_er + 3.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                                                                                            3.0 * dM_FeO_er + 3.0 * dM_FeSiO3_er))) + M_Mg * (
                      M_Fe * (M_FeO * (
                      M_SiO2 * (-dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                      -2.0 * dM_FeSiO3_er - 2.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er)) + M_FeSiO3 * (M_FeO * (
                      -3.0 * dM_FeO_er - 3.0 * dM_FeSiO3_er - 3.0 * dM_MgO_er - 3.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                                  -3.0 * dM_FeO_er - 3.0 * dM_FeSiO3_er - 3.0 * dM_MgO_er - 3.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                  2.0 * dM_FeSiO3_er + 2.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_m * (
                                                                                                  5.0 * dM_FeO_er + 5.0 * dM_MgO_er - 5.0 * dM_SiO2_er)) + M_MgO * (
                              M_SiO2 * (
                              -dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                              -2.0 * dM_FeSiO3_er - 2.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er)) + M_MgSiO3 * (M_FeO * (
                      -3.0 * dM_FeO_er - 3.0 * dM_FeSiO3_er - 3.0 * dM_MgO_er - 3.0 * dM_MgSiO3_er) + M_MgO * (
                                                                                                          -3.0 * dM_FeO_er - 3.0 * dM_FeSiO3_er - 3.0 * dM_MgO_er - 3.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                                          2.0 * dM_FeSiO3_er + 2.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_m * (
                                                                                                          5.0 * dM_FeO_er + 5.0 * dM_MgO_er - 5.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                              3.0 * dM_FeO_er + 3.0 * dM_FeSiO3_er + 3.0 * dM_MgO_er + 3.0 * dM_MgSiO3_er)) + M_FeO * (
                      M_MgO * (
                      M_SiO2 * (-dM_FeO_er - 2 * dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                      -2.0 * dM_FeSiO3_er - 2.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                      3.0 * dM_MgO_er + 3.0 * dM_MgSiO3_er)) + M_FeSiO3 * (M_FeO * (M_MgO * (
                      -3.0 * dM_FeO_er - 3.0 * dM_FeSiO3_er - 3.0 * dM_MgO_er - 3.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                    -3.0 * dM_MgO_er - 3.0 * dM_MgSiO3_er) + M_m * (
                                                                                    3.0 * dM_MgO_er + 3.0 * dM_MgSiO3_er)) + M_MgO * (
                                                                           M_SiO2 * (
                                                                           -dM_FeSiO3_er - dM_MgO_er - 2 * dM_MgSiO3_er - dM_SiO2_er) + M_m * (
                                                                           2.0 * dM_FeO_er - 2.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er)) + M_SiO2 * M_m * (
                                                                           3.0 * dM_MgO_er + 3.0 * dM_MgSiO3_er)) + M_MgSiO3 * (
                      M_FeO * (M_MgO * (
                      -3.0 * dM_FeO_er - 3.0 * dM_FeSiO3_er - 3.0 * dM_MgO_er - 3.0 * dM_MgSiO3_er) + M_SiO2 * (
                               2.0 * dM_FeO_er + 4.0 * dM_FeSiO3_er + 2.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_m * (
                               -5.0 * dM_FeSiO3_er + 5.0 * dM_MgO_er - 5.0 * dM_SiO2_er)) + M_FeSiO3 * (
                      M_FeO * (-3.0 * dM_FeO_er - 3.0 * dM_FeSiO3_er - 3.0 * dM_MgO_er - 3.0 * dM_MgSiO3_er) + M_MgO * (
                      -3.0 * dM_FeO_er - 3.0 * dM_FeSiO3_er - 3.0 * dM_MgO_er - 3.0 * dM_MgSiO3_er) + M_SiO2 * (
                      2.0 * dM_FeSiO3_er + 2.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er) + M_m * (
                      5.0 * dM_FeO_er + 5.0 * dM_MgO_er - 5.0 * dM_SiO2_er)))) + dKSiO2_KSiO2 * (M_Fe * (M_MgO * (
                      2.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                      M_FeO * (-2.0 * M_SiO2 + 2.0 * M_m) + 2.0 * M_SiO2 * M_m)) + M_MgSiO3 * (M_FeO * (
                      M_MgO * (-2.0 * M_SiO2 + 2.0 * M_m) + 2.0 * M_SiO2 * M_m) + M_FeSiO3 * (M_FeO * (
                      -2.0 * M_SiO2 + 2.0 * M_m) + M_MgO * (
                                                                                              -2.0 * M_SiO2 + 2.0 * M_m) + 2.0 * M_SiO2 * M_m))) + M_Mg * (
                                                                                                 M_Fe * (M_FeSiO3 * (
                                                                                                 M_FeO * (
                                                                                                 -2.0 * M_SiO2 + 2.0 * M_m) + M_MgO * (
                                                                                                 -2.0 * M_SiO2 + 2.0 * M_m) + 2.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                                         M_FeO * (
                                                                                                         -2.0 * M_SiO2 + 2.0 * M_m) + M_MgO * (
                                                                                                         -2.0 * M_SiO2 + 2.0 * M_m) + 2.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                                                                                                         2.0 * M_FeO + 2.0 * M_MgO)) + M_MgO * (
                                                                                                 2.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                 M_FeO * (
                                                                                                 -2.0 * M_SiO2 + 2.0 * M_m) + 2.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                 M_FeO * (M_MgO * (
                                                                                                 -2.0 * M_SiO2 + 2.0 * M_m) + 2.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                 M_FeO * (
                                                                                                 -2.0 * M_SiO2 + 2.0 * M_m) + M_MgO * (
                                                                                                 -2.0 * M_SiO2 + 2.0 * M_m) + 2.0 * M_SiO2 * M_m))) + M_c * (
                                                                                                 M_Fe * (M_MgO * (
                                                                                                 M_FeSiO3 * (
                                                                                                 2.0 * M_FeO - 3.0 * M_m) + M_SiO2 * (
                                                                                                 M_FeO - 2.0 * M_m)) + M_MgSiO3 * (
                                                                                                         M_FeO * (
                                                                                                         2.0 * M_MgO + M_SiO2) + M_FeSiO3 * (
                                                                                                         2.0 * M_FeO + 2.0 * M_MgO - 3.0 * M_m) + M_MgO * (
                                                                                                         2.0 * M_SiO2 - 2.0 * M_m) - 2.0 * M_SiO2 * M_m)) + M_Mg * (
                                                                                                 M_Fe * (M_FeSiO3 * (
                                                                                                 2.0 * M_FeO + 2.0 * M_MgO - 3.0 * M_m) + M_MgSiO3 * (
                                                                                                         2.0 * M_FeO + 2.0 * M_MgO - 3.0 * M_m) + M_SiO2 * (
                                                                                                         M_FeO + M_MgO - 2.0 * M_m)) + M_FeO * M_SiO2 * (
                                                                                                 M_MgO - 2.0 * M_m) + M_FeSiO3 * (
                                                                                                 M_FeO * (
                                                                                                 2.0 * M_MgO + 2.0 * M_SiO2 - 2.0 * M_m) + M_SiO2 * (
                                                                                                 M_MgO - 2.0 * M_m)) + M_MgSiO3 * (
                                                                                                 M_FeO * (
                                                                                                 2.0 * M_MgO - 3.0 * M_m) + M_FeSiO3 * (
                                                                                                 2.0 * M_FeO + 2.0 * M_MgO - 3.0 * M_m))) + M_MgO * (
                                                                                                 -2.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                                 M_FeO * (
                                                                                                 2.0 * M_SiO2 - 2.0 * M_m) - 2.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                                 M_FeO * (M_MgO * (
                                                                                                 2.0 * M_SiO2 - 2.0 * M_m) - 2.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                                 M_FeO * (
                                                                                                 2.0 * M_SiO2 - 2.0 * M_m) + M_MgO * (
                                                                                                 2.0 * M_SiO2 - 2.0 * M_m) - 2.0 * M_SiO2 * M_m))))) + M_c * (
                      M_Fe * (M_MgO * (M_FeSiO3 * (
                      M_FeO * (-dM_FeO_er - 1.0 * dM_FeSiO3_er - dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_m * (
                      dM_FeO_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er)) + M_SiO2 * (M_FeO * (
                      -dM_FeO_er - 2.0 * dM_FeSiO3_er - dM_MgO_er - 2.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
                                                                                      dM_FeO_er + 1.0 * dM_FeSiO3_er))) + M_MgSiO3 * (
                              M_FeO * (
                              M_MgO * (-dM_FeO_er - 1.0 * dM_FeSiO3_er - dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_SiO2 * (
                              -dM_FeO_er - 2.0 * dM_FeSiO3_er - 1.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er)) + M_FeSiO3 * (
                              M_FeO * (-dM_FeO_er - 1.0 * dM_FeSiO3_er - dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_MgO * (
                              -dM_FeO_er - 1.0 * dM_FeSiO3_er - dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_m * (
                              dM_FeO_er + dM_MgO_er - 1.0 * dM_SiO2_er)) + M_MgO * (
                              M_SiO2 * (-dM_FeO_er - 1.0 * dM_FeSiO3_er) + M_m * (
                              dM_FeO_er + 1.0 * dM_FeSiO3_er)) + M_SiO2 * M_m * (
                              dM_FeO_er + 1.0 * dM_FeSiO3_er))) + M_Mg * (M_Fe * (M_FeSiO3 * (
                      M_FeO * (-dM_FeO_er - 1.0 * dM_FeSiO3_er - dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_MgO * (
                      -dM_FeO_er - 1.0 * dM_FeSiO3_er - dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_m * (
                      dM_FeO_er + dM_MgO_er - 1.0 * dM_SiO2_er)) + M_MgSiO3 * (M_FeO * (
                      -dM_FeO_er - 1.0 * dM_FeSiO3_er - dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_MgO * (
                                                                               -dM_FeO_er - 1.0 * dM_FeSiO3_er - dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_m * (
                                                                               dM_FeO_er + dM_MgO_er - 1.0 * dM_SiO2_er)) + M_SiO2 * (
                                                                                  M_FeO * (
                                                                                  -dM_FeO_er - 2.0 * dM_FeSiO3_er - dM_MgO_er - 2.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_MgO * (
                                                                                  -dM_FeO_er - 2.0 * dM_FeSiO3_er - dM_MgO_er - 2.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
                                                                                  dM_FeO_er + 1.0 * dM_FeSiO3_er + dM_MgO_er + 1.0 * dM_MgSiO3_er))) + M_FeO * M_SiO2 * (
                                                                          M_MgO * (
                                                                          -dM_FeO_er - 2.0 * dM_FeSiO3_er - dM_MgO_er - 2.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
                                                                          dM_MgO_er + 1.0 * dM_MgSiO3_er)) + M_FeSiO3 * (
                                                                          M_FeO * (M_MgO * (
                                                                          -dM_FeO_er - 1.0 * dM_FeSiO3_er - dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_SiO2 * (
                                                                                   -dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_m * (
                                                                                   dM_MgO_er + 1.0 * dM_MgSiO3_er)) + M_SiO2 * (
                                                                          M_MgO * (
                                                                          -1.0 * dM_FeSiO3_er - dM_MgO_er - 2.0 * dM_MgSiO3_er - 1.0 * dM_SiO2_er) + M_m * (
                                                                          dM_MgO_er + 1.0 * dM_MgSiO3_er))) + M_MgSiO3 * (
                                                                          M_FeO * (M_MgO * (
                                                                          -dM_FeO_er - 1.0 * dM_FeSiO3_er - dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_m * (
                                                                                   -1.0 * dM_FeSiO3_er + dM_MgO_er - 1.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                          M_FeO * (
                                                                          -dM_FeO_er - 1.0 * dM_FeSiO3_er - dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_MgO * (
                                                                          -dM_FeO_er - 1.0 * dM_FeSiO3_er - dM_MgO_er - 1.0 * dM_MgSiO3_er) + M_m * (
                                                                          dM_FeO_er + dM_MgO_er - 1.0 * dM_SiO2_er)))) + M_Si * (
                      M_Fe * (M_MgO * (M_FeO * (
                      -dM_FeO_er - 3.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_FeSiO3 * (
                                       -dM_FeO_er - 3.0 * dM_FeSiO3_er - 3.0 * dM_MgO_er - 5.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_SiO2 * (
                                       -dM_FeO_er - 3.0 * dM_FeSiO3_er - 2 * dM_MgO_er - 4.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                                       dM_FeO_er + 3.0 * dM_FeSiO3_er + 2.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er)) + M_MgSiO3 * (
                              M_FeO * (
                              -dM_FeO_er - 3.0 * dM_FeSiO3_er + 1.0 * dM_MgO_er - 1.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_FeSiO3 * (
                              -dM_FeO_er - 3.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_SiO2 * (
                              -dM_FeO_er - 3.0 * dM_FeSiO3_er - 2.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                              dM_FeO_er + 3.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er + 2.0 * dM_SiO2_er))) + M_Mg * (M_Fe * (
                      M_FeO * (
                      -dM_FeO_er - 3.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_FeSiO3 * (
                      -dM_FeO_er - 3.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_MgO * (
                      -dM_FeO_er - 3.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_MgSiO3 * (
                      -dM_FeO_er - 3.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_SiO2 * (
                      -dM_FeO_er - 3.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                      dM_FeO_er + 3.0 * dM_FeSiO3_er + dM_MgO_er + 3.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er)) + M_FeO * (
                                                                                                               M_MgO * (
                                                                                                               -dM_FeO_er - 3.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                                               -2 * dM_FeO_er - 4.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                                                                                                               2.0 * dM_FeSiO3_er + dM_MgO_er + 3.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er)) + M_FeSiO3 * (
                                                                                                               M_MgO * (
                                                                                                               1.0 * dM_FeO_er - 1.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_SiO2 * (
                                                                                                               -2.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                                                                                                               -2.0 * dM_FeO_er + dM_MgO_er + 3.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er)) + M_MgSiO3 * (
                                                                                                               M_FeO * (
                                                                                                               -3.0 * dM_FeO_er - 5.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_FeSiO3 * (
                                                                                                               -dM_FeO_er - 3.0 * dM_FeSiO3_er - dM_MgO_er - 3.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er))) + M_MgO * (
                      M_FeO * (M_SiO2 * (
                      -2 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 2 * dM_MgO_er - 4.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                               2.0 * dM_FeSiO3_er + 2.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er)) + M_FeSiO3 * (
                      M_SiO2 * (-2.0 * dM_FeSiO3_er - 2 * dM_MgO_er - 4.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                      -2.0 * dM_FeO_er + 2.0 * dM_MgSiO3_er + 2.0 * dM_SiO2_er))) + M_MgSiO3 * (M_FeO * (
                      M_SiO2 * (-2 * dM_FeO_er - 4.0 * dM_FeSiO3_er - 2.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                      2.0 * dM_FeSiO3_er - 2.0 * dM_MgO_er + 2.0 * dM_SiO2_er)) + M_FeSiO3 * (M_SiO2 * (
                      -2.0 * dM_FeSiO3_er - 2.0 * dM_MgSiO3_er - 2.0 * dM_SiO2_er) + M_m * (
                                                                                              -2.0 * dM_FeO_er - 2.0 * dM_MgO_er + 2.0 * dM_SiO2_er)))))) / (
               M_O * (M_Fe * (M_MgO * (4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                              M_FeO * (M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                              M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                              -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m))) + M_Mg * (M_Fe * (M_FeSiO3 * (
               M_FeO * (-4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
               -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (-4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO)) + M_MgO * (
                                                                                           4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (M_MgO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + M_MgO * (
                                                                                           -4.0 * M_SiO2 + 4.0 * M_m) + 4.0 * M_SiO2 * M_m)))) + M_Si * (
               M_Fe * (
               M_MgO * (M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (M_MgO * (
               M_FeO * (-M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
               -9.0 * M_FeO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               -9.0 * M_MgO - M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                                                                                             -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                                                                                             -9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m))) + M_Mg * (
               M_Fe * (M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_O * (
                       M_FeO * (-M_SiO2 + 4.0 * M_m) + M_FeSiO3 * (
                       -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_MgO * (
                       -M_SiO2 + 4.0 * M_m) + M_MgSiO3 * (
                       -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_SiO2 * M_m * (
                       M_FeO + M_MgO)) + M_MgO * (
               M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-M_SiO2 + M_m) + M_MgO * (-M_SiO2 + M_m) + M_SiO2 * M_m)) + M_O * (
               M_FeO * (M_MgO * (-M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (-9.0 * M_MgO - 9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (
               -M_SiO2 + 4.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (-9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m) + M_FeSiO3 * (
               -9.0 * M_FeO - 9.0 * M_MgO - 4.0 * M_SiO2 + 25.0 * M_m)))) + M_O * (M_MgO * (
               9.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)) + M_MgSiO3 * (M_FeO * (
               M_MgO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m) + M_FeSiO3 * (M_FeO * (
               -9.0 * M_SiO2 + 9.0 * M_m) + M_MgO * (-9.0 * M_SiO2 + 9.0 * M_m) + 9.0 * M_SiO2 * M_m)))) + M_c * (
               M_Fe * (
               M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
               M_MgO * (M_FeSiO3 * (M_FeO - M_m) + M_SiO2 * (M_FeO - M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO + M_SiO2) + M_FeSiO3 * (M_FeO + M_MgO - M_m) + M_MgO * (
               M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Mg * (M_Fe * (
               M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_MgSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_O * (
               M_FeSiO3 * (M_FeO + M_MgO - M_m) + M_MgSiO3 * (M_FeO + M_MgO - M_m) + M_SiO2 * (
               M_FeO + M_MgO - M_m)) + M_SiO2 * M_m * (-M_FeO - M_MgO)) + M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
                                                         M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
                                                         M_FeO * (M_SiO2 - M_m) + M_MgO * (
                                                         M_SiO2 - M_m) - M_SiO2 * M_m)) + M_O * (
                                                         M_FeO * M_SiO2 * (M_MgO - M_m) + M_FeSiO3 * (
                                                         M_FeO * (M_MgO + M_SiO2 - M_m) + M_SiO2 * (
                                                         M_MgO - M_m)) + M_MgSiO3 * (
                                                         M_FeO * (M_MgO - M_m) + M_FeSiO3 * (
                                                         M_FeO + M_MgO - M_m)))) + M_O * (
               M_MgO * (-M_FeO * M_SiO2 * M_m + M_FeSiO3 * (M_FeO * (M_SiO2 - M_m) - M_SiO2 * M_m)) + M_MgSiO3 * (
               M_FeO * (M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m) + M_FeSiO3 * (
               M_FeO * (M_SiO2 - M_m) + M_MgO * (M_SiO2 - M_m) - M_SiO2 * M_m))) + M_Si * (M_Fe * (M_MgO * (
               M_FeO * (M_SiO2 - M_m) + M_FeSiO3 * (
               4.0 * M_FeO + M_SiO2 - 9.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (M_FeO * (
               4.0 * M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                     4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_O * (
                                                                                                   M_MgO * (
                                                                                                   M_FeO + M_FeSiO3 + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                   M_FeO + M_FeSiO3 + M_SiO2 - M_m))) + M_Mg * (
                                                                                           M_Fe * (M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_MgO * (
                                                                                                   M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                                   4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_O * (
                                                                                                   M_FeO + M_FeSiO3 + M_MgO + M_MgSiO3 + M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeO * (
                                                                                           M_MgO * (
                                                                                           M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_MgO + 4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                           M_SiO2 - M_m) - 4.0 * M_SiO2 * M_m) + M_MgSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_MgO + M_SiO2 - 9.0 * M_m) + M_FeSiO3 * (
                                                                                           4.0 * M_FeO + 4.0 * M_MgO + M_SiO2 - 9.0 * M_m)) + M_O * (
                                                                                           M_FeO * (
                                                                                           M_MgO + M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                           M_MgO + M_SiO2 - M_m) + M_MgSiO3 * (
                                                                                           M_FeO + M_FeSiO3))) + M_MgO * (
                                                                                           -4.0 * M_FeO * M_SiO2 * M_m + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (M_MgO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m) + M_FeSiO3 * (
                                                                                           M_FeO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) + M_MgO * (
                                                                                           4.0 * M_SiO2 - 4.0 * M_m) - 4.0 * M_SiO2 * M_m)) + M_O * (
                                                                                           M_MgO * (M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                                    M_SiO2 - M_m)) + M_MgSiO3 * (
                                                                                           M_FeO * (
                                                                                           M_SiO2 - M_m) + M_FeSiO3 * (
                                                                                           M_SiO2 - M_m))))))

    def dM_m_dTc_dMi(self, dM_im_dTc):
        '''alternate way to compute dM_m given the results of the mole changes of all mantle components
        inputs:
        dMi_m: [dM_MgO, dM_SiO2, dM_FeO, dM_MgSiO3, dM_FeSiO3]
        '''
        return np.sum(dM_im_dTc)

    def smooth_nonnegative(self, x):
        if x < 0:
            return 0
        elif x < 10:
            return x / (1 + 2 ** (-(x - 0.3) * 10))
        else:
            return x

    def smooth_nonpositive(self, x):
        if x > 0:
            return 0
        elif x > -10:
            return x / (1 + 2 ** ((x + 0.3) * 10))
        else:
            return x

    def logit(self, x, mu=None, sig=1):
        if mu is None:
            mu = sig * 10
        if x > mu + 10 * sig:
            return 1
        elif x < mu - 10 * sig:
            return 0
        else:
            return 1 / (1 + 2 ** -((x - mu) / sig))


