class Radiogenics(object):
    def __init__(self):
        self.CGyr2s = 1e9 * 365.25 * 24 * 3600
        self.l_238U = 4.9160e-18  # [1/Gyr > 1/s]
        self.l_235U = 3.1209e-17  # [1/Gyr > 1/s]
        self.l_232Th = 1.5633e-18  # [1/Gyr > 1/s]
        self.l_40K = 1.7558e-17  # [1/Gyr > 1/s]

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
        return Hp * (self._heat_hp(h_238U, self.l_238U, t, tp=tp)
                     + self._heat_hp(h_235U, self.l_235U, t, tp=tp)
                     + self._heat_hp(h_232Th, self.l_232Th, t, tp=tp)
                     + self._heat_hp(h_40K, self.l_40K, t, tp=tp))

    def heat_production_core(self, Hp, t, tp=1.44155e17):
        '''
        Calculates heat production through time given present-day heat production

        :param Hp: heat production at present day in [W, W/kg, W/m^3]
        :param t: time [s]
        :param tp: present-day time [s] default 1.44155e17 s (4.568 Byr)
        :return: Hc(t) heat production in core same units as provided in Hp
        '''
        h_40K = 1.  # normalized W/kg
        return Hp * self._heat_hp(h_40K, self.l_40K, t, tp=tp)

    def _heat_h0(self, h0, lam, t):
        return h0 * np.exp(-lam * t)

    def _heat_hp(self, hp, lam, t, tp=1.44155e17):
        return hp * np.exp(lam * (tp - t))

    def ppm2Hp_40K_Driscoll(self, ppm):
        '''
        Converts 40K ppm to present-day heat-production in W, using numbers from Nimmo 2004

        :param ppm:
        :return:
        '''
        return 2e12 / 255 * ppm

    def ppm2Hp_40K_Nimmo(self, ppm):
        '''
        Converts 40K ppm to present-day heat-production in W, using numbers from Driscoll 2014

        :param ppm:
        :return:
        '''
        return 2.1e12 / 300 * ppm

    def ppm2Hp_40K_Driscoll(self, Hp):
        '''
        Converts present-day heat-production in W to 40K ppm, using numbers from Driscoll 2014

        :param Hp:
        :return:
        '''
        return 255/2e12* Hp

    def Hp2ppm_Nimmo(self, Hp):
        '''
        Converts present-day heat-production in W to 40K ppm, using numbers from Nimmo 2004

        :param Hp:
        :return:
        '''
        return 300/2.1e12* Hp
