def unwrap_dKs(dKs):
    ''' helper function to unwrap dK values from dKs'''
    dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs
    return dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3

def unwrap_moles(Moles):
    ''' helper function to unwrap mol values from Moles'''
    M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = Moles
    return M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m

def dM_Mg(Moles, dKs):
    M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = Moles
    dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs
    return M_Mg*(M_Fe*dKFeO_KFeO*(M_MgO*(M_FeO*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m +
    M_c*(M_O + M_m))) + M_FeSiO3*M_SiO2*M_m*(-4*M_O + M_c) + M_MgSiO3*(M_FeO*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m))) + M_Si*(M_FeO*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c -
    M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(2*M_SiO2 + 10*M_m) - M_SiO2*M_m
    + M_c*(-M_O - M_SiO2 - 3*M_m)) + M_MgSiO3*(M_FeO*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c -
    M_m)))) + M_MgSiO3*(M_FeO*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*M_m*(-4*M_O + M_c)) +
    M_FeSiO3*M_m*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2)) + M_Si*(M_FeO*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) +
    M_O*(2*M_SiO2 + 10*M_m) - M_SiO2*M_m + M_c*(-M_O - M_SiO2 - 3*M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m +
    M_c*(-M_O + M_SiO2 - 9*M_m))))) + M_FeSiO3*dKFeSiO3_KFeSiO3*(M_MgO*(M_Fe*(M_MgSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O -
    M_SiO2 + M_m)) + M_SiO2*M_m*(-4*M_O + M_c)) + M_FeO*M_O*M_SiO2*M_c + M_Si*(M_Fe*(M_MgSiO3*(-9*M_O + M_SiO2 + 4*M_c -
    M_m) + M_O*(2*M_SiO2 + 10*M_m) - M_SiO2*M_m + M_c*(-M_O - M_SiO2 - 3*M_m)) + M_FeO*(M_O*(-3*M_SiO2 - 6*M_m) +
    M_c*(2*M_SiO2 + 2*M_m)))) + M_MgSiO3*(M_Fe*(M_FeO*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) +
    M_m*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2))) + M_FeO*M_O*M_c*M_m + M_Si*(M_Fe*(M_FeO*(-9*M_O + M_SiO2 + 4*M_c - M_m) +
    M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - 9*M_m)) + M_FeO*(M_O*(6*M_SiO2 - 15*M_m) + M_c*(-2*M_SiO2 +
    6*M_m))))) + M_MgSiO3*dKMgSiO3_KMgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) +
    M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_MgO*(M_Fe*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) - M_O*M_SiO2*M_c) + M_O*M_SiO2*M_c*(M_FeO +
    M_FeSiO3) + M_Si*(M_Fe*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(3*M_SiO2 + 6*M_m) + M_c*(-2*M_SiO2 - 2*M_m)) +
    M_FeO*(M_O*(-3*M_SiO2 - 6*M_m) + M_c*(2*M_SiO2 + 2*M_m)) + M_FeSiO3*(M_O*(-3*M_SiO2 - 6*M_m) + M_c*(2*M_SiO2 + 2*M_m))))
    + M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m) + M_Si*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O -
    M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 -
    25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) +
    M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + 9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)))) +
    M_Si*dKSiO2_KSiO2*(M_MgO*(M_Fe*(M_FeSiO3*(M_O*(-4*M_SiO2 + 10*M_m) + M_c*(-M_O + M_SiO2 - 3*M_m)) + M_SiO2*(6*M_O*M_m +
    M_c*(-M_O - 2*M_m))) + M_FeO*(M_FeSiO3*(M_O*(6*M_SiO2 - 6*M_m) + M_c*(-2*M_SiO2 + 2*M_m)) + M_SiO2*(-6*M_O*M_m +
    M_c*(M_O + 2*M_m))) + M_FeSiO3*M_SiO2*(-6*M_O*M_m + M_c*(M_O + 2*M_m)) + M_MgSiO3*(M_Fe*(M_O*(-6*M_SiO2 + 6*M_m) +
    M_c*(2*M_SiO2 - 2*M_m)) + M_FeO*(M_O*(6*M_SiO2 - 6*M_m) + M_c*(-2*M_SiO2 + 2*M_m)) + M_FeSiO3*(M_O*(6*M_SiO2 - 6*M_m) +
    M_c*(-2*M_SiO2 + 2*M_m)))) + M_MgSiO3*(M_Fe*(M_FeO*(M_O*(-2*M_SiO2 - 4*M_m) + M_c*(M_O + M_SiO2 + M_m)) +
    M_m*(6*M_O*M_SiO2 + M_c*(-M_O - 2*M_SiO2))) + M_FeO*(M_FeSiO3*(M_O*(6*M_SiO2 - 6*M_m) + M_c*(-2*M_SiO2 + 2*M_m)) +
    M_m*(-6*M_O*M_SiO2 + M_c*(M_O + 2*M_SiO2))) + M_FeSiO3*M_m*(-6*M_O*M_SiO2 + M_c*(M_O + 2*M_SiO2)))) +
    dKMgO_KMgO*(M_MgO*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_SiO2*(4*M_O*M_m +
    M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) +
    M_c*(-M_O + M_SiO2 - M_m)) + M_O*M_c*(-M_SiO2 + M_m)) + M_O*M_c*(M_FeO*(M_SiO2 - M_m) + M_FeSiO3*(M_SiO2 - M_m))) +
    M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m) + M_Si*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O -
    M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 -
    25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) +
    M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + 9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_MgSiO3*(M_Fe*(M_FeO*(9*M_O -
    M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(9*M_SiO2 - 9*M_m) + M_c*(-4*M_SiO2 + 4*M_m)) +
    M_FeO*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 -
    4*M_m))))) + M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_SiO2*(4*M_O*M_m +
    M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m) + M_Si*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O -
    M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 -
    25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) +
    M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + 9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m))))))/(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m
    + M_c*(M_O + M_m))) + M_m*(M_FeSiO3*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2)) - M_O*M_SiO2*M_c)) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_m*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2))) +
    M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) +
    M_FeSiO3*M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) + M_O*M_c*(M_FeO*(M_FeSiO3*(-M_SiO2 + M_m) + M_SiO2*M_m) +
    M_FeSiO3*M_SiO2*M_m)) + M_MgO*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) +
    M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_Mg*(M_Fe*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m + M_c*(M_O + M_m))) +
    M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) +
    M_FeSiO3*M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m)) + M_MgSiO3*(M_Fe*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) +
    M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)))) + M_MgSiO3*(M_Fe*(M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 +
    4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_O*M_c*(-M_SiO2 + M_m)) + M_O*M_c*(M_FeO*(M_SiO2 - M_m) + M_FeSiO3*(M_SiO2 -
    M_m))) + M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m) +
    M_Si*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)) + M_Mg*(M_Fe*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O +
    M_SiO2 - M_m)) + M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2
    + M_m)) + M_FeSiO3*(M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) + M_MgSiO3*(M_Fe*(-9*M_O + M_SiO2 +
    4*M_c - M_m) + M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m))) +
    M_MgSiO3*(M_Fe*(M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(9*M_SiO2 - 9*M_m)
    + M_c*(-4*M_SiO2 + 4*M_m)) + M_FeO*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + M_FeSiO3*(M_O*(-9*M_SiO2 +
    9*M_m) + M_c*(4*M_SiO2 - 4*M_m))))) + M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m)) +
    M_Si*(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O +
    M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - 9*M_m)) + 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(9*M_SiO2 - 9*M_m) + M_c*(-4*M_SiO2 + 4*M_m)) -
    9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeSiO3*(-9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) +
    4*M_SiO2*M_m)) + M_MgSiO3*(M_Fe*(M_FeO*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m +
    M_c*(-M_O + M_SiO2 - 9*M_m)) + M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m +
    M_c*(M_O - M_SiO2 + 9*M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)))) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)))))

def dM_Si(Moles, dKs):
    M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = Moles
    dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs
    return M_Si*(M_Fe*dKFeO_KFeO*(M_Mg*(M_FeO*(M_FeSiO3*(M_O*(-6*M_SiO2 + 6*M_m) + M_c*(2*M_SiO2 - 2*M_m)) +
    M_MgSiO3*(M_O*(-4*M_SiO2 + 10*M_m) + M_c*(-M_O + M_SiO2 - 3*M_m)) + M_SiO2*(6*M_O*M_m + M_c*(-M_O - 2*M_m))) +
    M_FeSiO3*M_m*(6*M_O*M_SiO2 + M_c*(-M_O - 2*M_SiO2))) + M_MgO*(M_FeO*(M_FeSiO3*(M_O*(6*M_SiO2 - 6*M_m) + M_c*(-2*M_SiO2 +
    2*M_m)) + M_SiO2*(-6*M_O*M_m + M_c*(M_O + 2*M_m))) + M_FeSiO3*(M_Mg*(M_O*(-2*M_SiO2 - 4*M_m) + M_c*(M_O + M_SiO2 + M_m))
    + M_m*(-6*M_O*M_SiO2 + M_c*(M_O + 2*M_SiO2))) + M_MgSiO3*(M_FeO*(M_O*(6*M_SiO2 - 6*M_m) + M_c*(-2*M_SiO2 + 2*M_m)) +
    M_FeSiO3*(M_O*(6*M_SiO2 - 6*M_m) + M_c*(-2*M_SiO2 + 2*M_m)))) + M_MgSiO3*(M_FeO*(M_FeSiO3*(M_O*(6*M_SiO2 - 6*M_m) +
    M_c*(-2*M_SiO2 + 2*M_m)) + M_SiO2*(-6*M_O*M_m + M_c*(M_O + 2*M_m))) + M_FeSiO3*M_m*(-6*M_O*M_SiO2 + M_c*(M_O +
    2*M_SiO2)))) + M_FeSiO3*dKFeSiO3_KFeSiO3*(M_Mg*(M_Fe*(M_FeO*(M_O*(-2*M_SiO2 - 4*M_m) + M_c*(M_O + M_SiO2 + M_m)) +
    M_m*(6*M_O*M_SiO2 + M_c*(-M_O - 2*M_SiO2))) + M_FeO*(M_MgSiO3*(M_O*(-4*M_SiO2 + 10*M_m) + M_c*(-M_O + M_SiO2 - 3*M_m)) +
    M_O*M_c*(-M_SiO2 + M_m))) + M_MgO*(M_Fe*(M_FeO*(M_O*(2*M_SiO2 + 4*M_m) + M_c*(-M_O - M_SiO2 - M_m)) +
    M_MgSiO3*(M_O*(6*M_SiO2 - 6*M_m) + M_c*(-2*M_SiO2 + 2*M_m)) + M_m*(-6*M_O*M_SiO2 + M_c*(M_O + 2*M_SiO2))) +
    M_FeO*M_O*M_c*(M_SiO2 - M_m) + M_Mg*(M_Fe*(M_O*(-2*M_SiO2 - 4*M_m) + M_c*(M_O + M_SiO2 + M_m)) + M_FeO*(M_O*(2*M_SiO2 +
    4*M_m) + M_c*(-M_O - M_SiO2 - M_m)))) + M_MgSiO3*(M_Fe*(M_FeO*(M_O*(2*M_SiO2 + 4*M_m) + M_c*(-M_O - M_SiO2 - M_m)) +
    M_m*(-6*M_O*M_SiO2 + M_c*(M_O + 2*M_SiO2))) + M_FeO*M_O*M_c*(M_SiO2 - M_m))) +
    M_Mg*dKMgO_KMgO*(M_MgO*(M_Fe*(M_FeSiO3*(M_O*(-4*M_SiO2 + 10*M_m) + M_c*(-M_O + M_SiO2 - 3*M_m)) + M_SiO2*(6*M_O*M_m +
    M_c*(-M_O - 2*M_m))) + M_FeO*(M_FeSiO3*(M_O*(6*M_SiO2 - 6*M_m) + M_c*(-2*M_SiO2 + 2*M_m)) + M_SiO2*(-6*M_O*M_m +
    M_c*(M_O + 2*M_m))) + M_FeSiO3*M_SiO2*(-6*M_O*M_m + M_c*(M_O + 2*M_m)) + M_MgSiO3*(M_Fe*(M_O*(-6*M_SiO2 + 6*M_m) +
    M_c*(2*M_SiO2 - 2*M_m)) + M_FeO*(M_O*(6*M_SiO2 - 6*M_m) + M_c*(-2*M_SiO2 + 2*M_m)) + M_FeSiO3*(M_O*(6*M_SiO2 - 6*M_m) +
    M_c*(-2*M_SiO2 + 2*M_m)))) + M_MgSiO3*(M_Fe*(M_FeO*(M_O*(-2*M_SiO2 - 4*M_m) + M_c*(M_O + M_SiO2 + M_m)) +
    M_m*(6*M_O*M_SiO2 + M_c*(-M_O - 2*M_SiO2))) + M_FeO*(M_FeSiO3*(M_O*(6*M_SiO2 - 6*M_m) + M_c*(-2*M_SiO2 + 2*M_m)) +
    M_m*(-6*M_O*M_SiO2 + M_c*(M_O + 2*M_SiO2))) + M_FeSiO3*M_m*(-6*M_O*M_SiO2 + M_c*(M_O + 2*M_SiO2)))) +
    M_MgSiO3*dKMgSiO3_KMgSiO3*(M_Mg*(M_Fe*(M_FeO*(M_O*(-2*M_SiO2 - 4*M_m) + M_c*(M_O + M_SiO2 + M_m)) + M_m*(6*M_O*M_SiO2 +
    M_c*(-M_O - 2*M_SiO2))) + M_FeO*(M_FeSiO3*(M_O*(6*M_SiO2 - 6*M_m) + M_c*(-2*M_SiO2 + 2*M_m)) + M_m*(-6*M_O*M_SiO2 +
    M_c*(M_O + 2*M_SiO2))) + M_FeSiO3*M_m*(-6*M_O*M_SiO2 + M_c*(M_O + 2*M_SiO2))) + M_MgO*(M_Fe*(M_FeO*(M_O*(2*M_SiO2 +
    4*M_m) + M_c*(-M_O - M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 10*M_m) + M_c*(-M_O + M_SiO2 - 3*M_m)) +
    M_O*M_c*(-M_SiO2 + M_m)) + M_Mg*(M_Fe*(M_O*(-2*M_SiO2 - 4*M_m) + M_c*(M_O + M_SiO2 + M_m)) + M_FeO*(M_O*(2*M_SiO2 +
    4*M_m) + M_c*(-M_O - M_SiO2 - M_m)) + M_FeSiO3*(M_O*(2*M_SiO2 + 4*M_m) + M_c*(-M_O - M_SiO2 - M_m))) +
    M_O*M_c*(M_FeO*(M_SiO2 - M_m) + M_FeSiO3*(M_SiO2 - M_m)))) + dKSiO2_KSiO2*(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(4*M_SiO2 -
    4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m + M_c*(M_O + M_m))) + M_m*(M_FeSiO3*(-4*M_O*M_SiO2 + M_c*(-M_O +
    M_SiO2)) - M_O*M_SiO2*M_c)) + M_MgSiO3*(M_Fe*(M_FeO*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) +
    M_m*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2))) + M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) +
    M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) + M_FeSiO3*M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) +
    M_O*M_c*(M_FeO*(M_FeSiO3*(-M_SiO2 + M_m) + M_SiO2*M_m) + M_FeSiO3*M_SiO2*M_m)) +
    M_MgO*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O -
    M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) + M_Mg*(M_Fe*(M_FeSiO3*(M_O*(4*M_SiO2 -
    4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m + M_c*(M_O + M_m))) + M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) +
    M_c*(-M_O + M_SiO2 - M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_FeSiO3*M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m)) +
    M_MgSiO3*(M_Fe*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2
    - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)))) + M_MgSiO3*(M_Fe*(M_FeO*(M_O*(-4*M_SiO2 +
    4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_O*M_c*(-M_SiO2
    + M_m)) + M_O*M_c*(M_FeO*(M_SiO2 - M_m) + M_FeSiO3*(M_SiO2 - M_m))) + M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) -
    M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m)) + M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) -
    M_FeSiO3*M_SiO2*M_m))))/(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) +
    M_SiO2*(-4*M_O*M_m + M_c*(M_O + M_m))) + M_m*(M_FeSiO3*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2)) - M_O*M_SiO2*M_c)) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_m*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2))) +
    M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) +
    M_FeSiO3*M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) + M_O*M_c*(M_FeO*(M_FeSiO3*(-M_SiO2 + M_m) + M_SiO2*M_m) +
    M_FeSiO3*M_SiO2*M_m)) + M_MgO*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) +
    M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_Mg*(M_Fe*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m + M_c*(M_O + M_m))) +
    M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) +
    M_FeSiO3*M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m)) + M_MgSiO3*(M_Fe*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) +
    M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)))) + M_MgSiO3*(M_Fe*(M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 +
    4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_O*M_c*(-M_SiO2 + M_m)) + M_O*M_c*(M_FeO*(M_SiO2 - M_m) + M_FeSiO3*(M_SiO2 -
    M_m))) + M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m) +
    M_Si*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)) + M_Mg*(M_Fe*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O +
    M_SiO2 - M_m)) + M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2
    + M_m)) + M_FeSiO3*(M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) + M_MgSiO3*(M_Fe*(-9*M_O + M_SiO2 +
    4*M_c - M_m) + M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m))) +
    M_MgSiO3*(M_Fe*(M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(9*M_SiO2 - 9*M_m)
    + M_c*(-4*M_SiO2 + 4*M_m)) + M_FeO*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + M_FeSiO3*(M_O*(-9*M_SiO2 +
    9*M_m) + M_c*(4*M_SiO2 - 4*M_m))))) + M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m)) +
    M_Si*(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O +
    M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - 9*M_m)) + 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(9*M_SiO2 - 9*M_m) + M_c*(-4*M_SiO2 + 4*M_m)) -
    9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeSiO3*(-9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) +
    4*M_SiO2*M_m)) + M_MgSiO3*(M_Fe*(M_FeO*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m +
    M_c*(-M_O + M_SiO2 - 9*M_m)) + M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m +
    M_c*(M_O - M_SiO2 + 9*M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)))) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)))))

def dM_Fe(Moles, dKs):
    M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = Moles
    dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs
    return M_Fe*(M_FeSiO3*dKFeSiO3_KFeSiO3*(M_Mg*(M_MgSiO3*(M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) +
    M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) + M_O*M_SiO2*M_c*(-M_FeO + M_m)) + M_MgO*(M_Mg*(M_MgSiO3*(M_O*(-4*M_SiO2 +
    4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_O*M_c*(M_MgSiO3*(M_SiO2 - M_m) +
    M_SiO2*(M_FeO - M_m)) + M_Si*(M_FeO*(M_O*(-3*M_SiO2 - 6*M_m) + M_c*(2*M_SiO2 + 2*M_m)) + M_Mg*(M_MgSiO3*(9*M_O - M_SiO2
    - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) + M_MgSiO3*(M_O*(-9*M_SiO2 + 9*M_m) +
    M_c*(4*M_SiO2 - 4*M_m)) + 9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m))) + M_MgSiO3*M_O*M_SiO2*M_c*(M_FeO
    - M_m) + M_Si*(M_Mg*(M_FeO*(M_O*(3*M_SiO2 + 6*M_m) + M_c*(-2*M_SiO2 - 2*M_m)) + M_MgSiO3*(M_FeO*(9*M_O - M_SiO2 - 4*M_c
    + M_m) + M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m)
    + 4*M_SiO2*M_m)) + M_MgSiO3*(M_FeO*(M_O*(-3*M_SiO2 - 6*M_m) + M_c*(2*M_SiO2 + 2*M_m)) + 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)))) + M_Mg*dKMgO_KMgO*(M_MgO*(M_FeO*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O
    - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m + M_c*(M_O + M_m))) + M_FeSiO3*M_SiO2*M_m*(-4*M_O + M_c) +
    M_MgSiO3*(M_FeO*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O -
    M_SiO2 + M_m))) + M_Si*(M_FeO*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O
    + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(2*M_SiO2 + 10*M_m) - M_SiO2*M_m + M_c*(-M_O - M_SiO2 - 3*M_m)) +
    M_MgSiO3*(M_FeO*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m)))) +
    M_MgSiO3*(M_FeO*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*M_m*(-4*M_O + M_c)) +
    M_FeSiO3*M_m*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2)) + M_Si*(M_FeO*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) +
    M_O*(2*M_SiO2 + 10*M_m) - M_SiO2*M_m + M_c*(-M_O - M_SiO2 - 3*M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m +
    M_c*(-M_O + M_SiO2 - 9*M_m))))) + M_MgSiO3*dKMgSiO3_KMgSiO3*(M_Mg*(M_FeO*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O -
    M_SiO2 + M_m)) + M_SiO2*M_m*(-4*M_O + M_c)) + M_FeSiO3*M_m*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2)) +
    M_Si*(M_FeO*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(2*M_SiO2 + 10*M_m) - M_SiO2*M_m + M_c*(-M_O - M_SiO2 -
    3*M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - 9*M_m)))) + M_MgO*(M_FeO*M_O*M_SiO2*M_c
    + M_FeSiO3*(M_Mg*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_O*M_c*M_m) + M_Si*(M_FeO*(M_O*(-3*M_SiO2 -
    6*M_m) + M_c*(2*M_SiO2 + 2*M_m)) + M_FeSiO3*(M_Mg*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(6*M_SiO2 - 15*M_m) +
    M_c*(-2*M_SiO2 + 6*M_m))))) + M_Si*dKSiO2_KSiO2*(M_Mg*(M_FeO*(M_FeSiO3*(M_O*(-6*M_SiO2 + 6*M_m) + M_c*(2*M_SiO2 -
    2*M_m)) + M_MgSiO3*(M_O*(-4*M_SiO2 + 10*M_m) + M_c*(-M_O + M_SiO2 - 3*M_m)) + M_SiO2*(6*M_O*M_m + M_c*(-M_O - 2*M_m))) +
    M_FeSiO3*M_m*(6*M_O*M_SiO2 + M_c*(-M_O - 2*M_SiO2))) + M_MgO*(M_FeO*(M_FeSiO3*(M_O*(6*M_SiO2 - 6*M_m) + M_c*(-2*M_SiO2 +
    2*M_m)) + M_SiO2*(-6*M_O*M_m + M_c*(M_O + 2*M_m))) + M_FeSiO3*(M_Mg*(M_O*(-2*M_SiO2 - 4*M_m) + M_c*(M_O + M_SiO2 + M_m))
    + M_m*(-6*M_O*M_SiO2 + M_c*(M_O + 2*M_SiO2))) + M_MgSiO3*(M_FeO*(M_O*(6*M_SiO2 - 6*M_m) + M_c*(-2*M_SiO2 + 2*M_m)) +
    M_FeSiO3*(M_O*(6*M_SiO2 - 6*M_m) + M_c*(-2*M_SiO2 + 2*M_m)))) + M_MgSiO3*(M_FeO*(M_FeSiO3*(M_O*(6*M_SiO2 - 6*M_m) +
    M_c*(-2*M_SiO2 + 2*M_m)) + M_SiO2*(-6*M_O*M_m + M_c*(M_O + 2*M_m))) + M_FeSiO3*M_m*(-6*M_O*M_SiO2 + M_c*(M_O +
    2*M_SiO2)))) + dKFeO_KFeO*(M_Mg*(M_MgSiO3*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) +
    M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) + M_FeSiO3*M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) +
    M_O*M_c*(M_FeO*(M_FeSiO3*(-M_SiO2 + M_m) + M_SiO2*M_m) + M_FeSiO3*M_SiO2*M_m)) +
    M_MgO*(M_Mg*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O -
    M_m))) + M_FeSiO3*M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m)) + M_MgSiO3*(M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2
    - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)))) + M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) -
    M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m + M_MgSiO3*(M_FeO*(M_SiO2 - M_m) + M_FeSiO3*(M_SiO2 - M_m))) +
    M_Si*(M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + 9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_Mg*(M_FeO*(M_FeSiO3*(9*M_O -
    M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) + M_FeSiO3*(M_O*(M_SiO2 - 4*M_m) +
    M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) + M_MgSiO3*(M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 -
    4*M_c + M_m))) + M_MgSiO3*(M_FeO*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m)
    + M_c*(4*M_SiO2 - 4*M_m))))) + M_MgSiO3*M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m) +
    M_Si*(M_Mg*(M_FeO*(M_FeSiO3*(M_O*(9*M_SiO2 - 9*M_m) + M_c*(-4*M_SiO2 + 4*M_m)) - 9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 -
    M_m) + 4*M_SiO2*M_m)) + M_FeSiO3*(-9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) +
    M_MgSiO3*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    9*M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)))) +
    M_MgSiO3*(M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + 9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 +
    M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m))))))/(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m
    + M_c*(M_O + M_m))) + M_m*(M_FeSiO3*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2)) - M_O*M_SiO2*M_c)) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_m*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2))) +
    M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) +
    M_FeSiO3*M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) + M_O*M_c*(M_FeO*(M_FeSiO3*(-M_SiO2 + M_m) + M_SiO2*M_m) +
    M_FeSiO3*M_SiO2*M_m)) + M_MgO*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) +
    M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_Mg*(M_Fe*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m + M_c*(M_O + M_m))) +
    M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) +
    M_FeSiO3*M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m)) + M_MgSiO3*(M_Fe*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) +
    M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)))) + M_MgSiO3*(M_Fe*(M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 +
    4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_O*M_c*(-M_SiO2 + M_m)) + M_O*M_c*(M_FeO*(M_SiO2 - M_m) + M_FeSiO3*(M_SiO2 -
    M_m))) + M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m) +
    M_Si*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)) + M_Mg*(M_Fe*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O +
    M_SiO2 - M_m)) + M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2
    + M_m)) + M_FeSiO3*(M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) + M_MgSiO3*(M_Fe*(-9*M_O + M_SiO2 +
    4*M_c - M_m) + M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m))) +
    M_MgSiO3*(M_Fe*(M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(9*M_SiO2 - 9*M_m)
    + M_c*(-4*M_SiO2 + 4*M_m)) + M_FeO*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + M_FeSiO3*(M_O*(-9*M_SiO2 +
    9*M_m) + M_c*(4*M_SiO2 - 4*M_m))))) + M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m)) +
    M_Si*(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O +
    M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - 9*M_m)) + 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(9*M_SiO2 - 9*M_m) + M_c*(-4*M_SiO2 + 4*M_m)) -
    9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeSiO3*(-9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) +
    4*M_SiO2*M_m)) + M_MgSiO3*(M_Fe*(M_FeO*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m +
    M_c*(-M_O + M_SiO2 - 9*M_m)) + M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m +
    M_c*(M_O - M_SiO2 + 9*M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)))) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)))))

def dM_O(Moles, dKs):
    M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = Moles
    dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs
    return M_O*(M_Fe*dKFeO_KFeO*(M_MgO*(M_Si*(M_FeO*(M_FeSiO3*(3*M_SiO2 - 3*M_m) - 3*M_SiO2*M_m + M_c*(M_SiO2 + M_m)) +
    M_FeSiO3*(M_Mg*(-M_SiO2 + 2*M_c - 2*M_m) - 3*M_SiO2*M_m + M_c*(-M_SiO2 + 3*M_m)) + M_MgSiO3*(M_FeO*(3*M_SiO2 - 3*M_m) +
    M_FeSiO3*(3*M_SiO2 - 3*M_m))) + M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) + M_FeSiO3*M_SiO2*(-M_Mg - M_m) +
    M_MgSiO3*(M_FeO*(M_SiO2 - M_m) + M_FeSiO3*(M_SiO2 - M_m)))) + M_Si*(M_Mg*(M_FeO*(M_FeSiO3*(-3*M_SiO2 + 3*M_m) +
    M_MgSiO3*(-2*M_SiO2 - 2*M_c + 5*M_m) + 3*M_SiO2*M_m + M_c*(-M_SiO2 - M_m)) + M_FeSiO3*(3*M_SiO2*M_m + M_c*(M_SiO2 -
    3*M_m))) + M_MgSiO3*(M_FeO*(M_FeSiO3*(3*M_SiO2 - 3*M_m) - 3*M_SiO2*M_m + M_c*(M_SiO2 + M_m)) + M_FeSiO3*(-3*M_SiO2*M_m +
    M_c*(-M_SiO2 + 3*M_m)))) + M_c*(M_Mg*(M_FeO*(M_FeSiO3*(-M_SiO2 + M_m) + M_m*(M_MgSiO3 + M_SiO2)) + M_FeSiO3*M_SiO2*M_m)
    + M_MgSiO3*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m))) +
    M_FeSiO3*dKFeSiO3_KFeSiO3*(M_MgO*(M_Si*(M_Fe*(M_FeO*(M_SiO2 - 2*M_c + 2*M_m) + M_MgSiO3*(3*M_SiO2 - 3*M_m) -
    3*M_SiO2*M_m + M_c*(-M_SiO2 + 3*M_m)) + M_FeO*M_c*(2*M_SiO2 - 2*M_m) + M_Mg*(M_Fe*(-M_SiO2 + 2*M_c - 2*M_m) +
    M_FeO*(M_SiO2 - 2*M_c + 2*M_m))) + M_c*(M_Fe*(M_MgSiO3*(M_SiO2 - M_m) + M_SiO2*(M_FeO - M_m)) + M_Mg*M_SiO2*(-M_Fe +
    M_FeO))) + M_Si*(M_Mg*(M_Fe*(M_FeO*(-M_SiO2 + 2*M_c - 2*M_m) + 3*M_SiO2*M_m + M_c*(M_SiO2 - 3*M_m)) +
    M_FeO*(M_MgSiO3*(-2*M_SiO2 - 2*M_c + 5*M_m) + M_c*(-2*M_SiO2 + 2*M_m))) + M_MgSiO3*(M_Fe*(M_FeO*(M_SiO2 - 2*M_c + 2*M_m)
    - 3*M_SiO2*M_m + M_c*(-M_SiO2 + 3*M_m)) + M_FeO*M_c*(2*M_SiO2 - 2*M_m))) + M_c*(M_Fe*M_MgSiO3*M_SiO2*(M_FeO - M_m) +
    M_Mg*(M_Fe*M_SiO2*(-M_FeO + M_m) + M_FeO*M_MgSiO3*M_m))) + M_Mg*dKMgO_KMgO*(M_MgO*(M_Si*(M_Fe*(M_FeSiO3*(-2*M_SiO2 -
    2*M_c + 5*M_m) + 3*M_SiO2*M_m + M_c*(-M_SiO2 - M_m)) + M_FeO*(M_FeSiO3*(3*M_SiO2 - 3*M_m) - 3*M_SiO2*M_m + M_c*(M_SiO2 +
    M_m)) + M_FeSiO3*(-3*M_SiO2*M_m + M_c*(M_SiO2 + M_m)) + M_MgSiO3*(M_Fe*(-3*M_SiO2 + 3*M_m) + M_FeO*(3*M_SiO2 - 3*M_m) +
    M_FeSiO3*(3*M_SiO2 - 3*M_m))) + M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) + M_MgSiO3*(M_Fe*(-M_SiO2 + M_m) +
    M_FeO*(M_SiO2 - M_m) + M_FeSiO3*(M_SiO2 - M_m)) + M_m*(M_Fe*(M_FeSiO3 + M_SiO2) - M_FeSiO3*M_SiO2))) +
    M_MgSiO3*(M_Si*(M_Fe*(M_FeO*(-M_SiO2 + 2*M_c - 2*M_m) + 3*M_SiO2*M_m + M_c*(M_SiO2 - 3*M_m)) + M_FeO*(M_FeSiO3*(3*M_SiO2
    - 3*M_m) - 3*M_SiO2*M_m + M_c*(-M_SiO2 + 3*M_m)) + M_FeSiO3*(-3*M_SiO2*M_m + M_c*(-M_SiO2 + 3*M_m))) +
    M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) + M_SiO2*(M_Fe*(-M_FeO + M_m) - M_FeSiO3*M_m)))) +
    M_MgSiO3*dKMgSiO3_KMgSiO3*(M_Mg*(M_Si*(M_Fe*(M_FeO*(-M_SiO2 + 2*M_c - 2*M_m) + 3*M_SiO2*M_m + M_c*(M_SiO2 - 3*M_m)) +
    M_FeO*(M_FeSiO3*(3*M_SiO2 - 3*M_m) - 3*M_SiO2*M_m + M_c*(-M_SiO2 + 3*M_m)) + M_FeSiO3*(-3*M_SiO2*M_m + M_c*(-M_SiO2 +
    3*M_m))) + M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) + M_SiO2*(M_Fe*(-M_FeO + M_m) - M_FeSiO3*M_m))) +
    M_MgO*(M_Si*(M_Fe*(M_FeO*(M_SiO2 - 2*M_c + 2*M_m) + M_FeSiO3*(-2*M_SiO2 - 2*M_c + 5*M_m) + M_c*(-2*M_SiO2 + 2*M_m)) +
    M_Mg*(M_Fe*(-M_SiO2 + 2*M_c - 2*M_m) + M_FeO*(M_SiO2 - 2*M_c + 2*M_m) + M_FeSiO3*(M_SiO2 - 2*M_c + 2*M_m)) +
    M_c*(M_FeO*(2*M_SiO2 - 2*M_m) + M_FeSiO3*(2*M_SiO2 - 2*M_m))) + M_c*(M_Fe*(M_FeO*M_SiO2 + M_FeSiO3*M_m) +
    M_Mg*M_SiO2*(-M_Fe + M_FeO + M_FeSiO3)))) + M_Si*dKSiO2_KSiO2*(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(2*M_SiO2 + 2*M_c - 2*M_m) +
    M_SiO2*(M_c - 2*M_m)) + M_m*(M_FeSiO3*(-2*M_SiO2 - 3*M_c) - 2*M_SiO2*M_c)) + M_MgSiO3*(M_Fe*(M_FeO*(2*M_SiO2 + 2*M_c -
    2*M_m) + M_m*(-2*M_SiO2 - 3*M_c)) + M_FeO*(M_FeSiO3*(-2*M_SiO2 - 2*M_c + 2*M_m) + M_m*(2*M_SiO2 + 3*M_c)) +
    M_FeSiO3*M_m*(2*M_SiO2 + 3*M_c)) + M_c*(M_FeO*(M_FeSiO3*(-2*M_SiO2 + 2*M_m) + 2*M_SiO2*M_m) + 2*M_FeSiO3*M_SiO2*M_m)) +
    M_MgO*(M_Fe*(M_FeO*(M_FeSiO3*(-2*M_SiO2 - 2*M_c + 2*M_m) + M_SiO2*(-M_c + 2*M_m)) + M_m*(M_FeSiO3*(2*M_SiO2 + 3*M_c) +
    2*M_SiO2*M_c)) + M_Mg*(M_Fe*(M_FeSiO3*(2*M_SiO2 + 2*M_c - 2*M_m) + M_SiO2*(M_c - 2*M_m)) + M_FeO*(M_FeSiO3*(-2*M_SiO2 -
    2*M_c + 2*M_m) + M_SiO2*(-M_c + 2*M_m)) + M_FeSiO3*M_SiO2*(-M_c + 2*M_m) + M_MgSiO3*(M_Fe*(2*M_SiO2 + 2*M_c - 2*M_m) +
    M_FeO*(-2*M_SiO2 - 2*M_c + 2*M_m) + M_FeSiO3*(-2*M_SiO2 - 2*M_c + 2*M_m))) + M_MgSiO3*(M_Fe*(M_FeO*(-2*M_SiO2 - 2*M_c +
    2*M_m) + M_FeSiO3*(-2*M_SiO2 - 2*M_c + 2*M_m) + M_c*(-2*M_SiO2 + 2*M_m)) + M_c*(M_FeO*(2*M_SiO2 - 2*M_m) +
    M_FeSiO3*(2*M_SiO2 - 2*M_m))) + M_c*(M_FeO*(M_FeSiO3*(2*M_SiO2 - 2*M_m) - 2*M_SiO2*M_m) - 2*M_FeSiO3*M_SiO2*M_m)) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(-2*M_SiO2 - 2*M_c + 2*M_m) + M_SiO2*(-M_c + 2*M_m)) + M_m*(M_FeSiO3*(2*M_SiO2 + 3*M_c)
    + 2*M_SiO2*M_c)) + M_c*(M_FeO*(M_FeSiO3*(2*M_SiO2 - 2*M_m) - 2*M_SiO2*M_m) -
    2*M_FeSiO3*M_SiO2*M_m))))/(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) +
    M_SiO2*(-4*M_O*M_m + M_c*(M_O + M_m))) + M_m*(M_FeSiO3*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2)) - M_O*M_SiO2*M_c)) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_m*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2))) +
    M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) +
    M_FeSiO3*M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) + M_O*M_c*(M_FeO*(M_FeSiO3*(-M_SiO2 + M_m) + M_SiO2*M_m) +
    M_FeSiO3*M_SiO2*M_m)) + M_MgO*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) +
    M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_Mg*(M_Fe*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m + M_c*(M_O + M_m))) +
    M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) +
    M_FeSiO3*M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m)) + M_MgSiO3*(M_Fe*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) +
    M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)))) + M_MgSiO3*(M_Fe*(M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 +
    4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_O*M_c*(-M_SiO2 + M_m)) + M_O*M_c*(M_FeO*(M_SiO2 - M_m) + M_FeSiO3*(M_SiO2 -
    M_m))) + M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m) +
    M_Si*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)) + M_Mg*(M_Fe*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O +
    M_SiO2 - M_m)) + M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2
    + M_m)) + M_FeSiO3*(M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) + M_MgSiO3*(M_Fe*(-9*M_O + M_SiO2 +
    4*M_c - M_m) + M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m))) +
    M_MgSiO3*(M_Fe*(M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(9*M_SiO2 - 9*M_m)
    + M_c*(-4*M_SiO2 + 4*M_m)) + M_FeO*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + M_FeSiO3*(M_O*(-9*M_SiO2 +
    9*M_m) + M_c*(4*M_SiO2 - 4*M_m))))) + M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m)) +
    M_Si*(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O +
    M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - 9*M_m)) + 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(9*M_SiO2 - 9*M_m) + M_c*(-4*M_SiO2 + 4*M_m)) -
    9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeSiO3*(-9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) +
    4*M_SiO2*M_m)) + M_MgSiO3*(M_Fe*(M_FeO*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m +
    M_c*(-M_O + M_SiO2 - 9*M_m)) + M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m +
    M_c*(M_O - M_SiO2 + 9*M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)))) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)))))

def dM_c(Moles, dKs):
    M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = Moles
    dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs
    return M_c*(M_Fe*dKFeO_KFeO*(M_MgO*(M_O*(M_FeO*(M_FeSiO3*(2*M_SiO2 - 2*M_m) - 2*M_SiO2*M_m) + M_FeSiO3*M_SiO2*(-2*M_Mg - 2*M_m)
    + M_MgSiO3*(M_FeO*(2*M_SiO2 - 2*M_m) + M_FeSiO3*(2*M_SiO2 - 2*M_m))) + M_Si*(M_FeO*(M_FeSiO3*(2*M_SiO2 - 2*M_m) +
    M_O*(M_SiO2 + 2*M_m) - 2*M_SiO2*M_m) + M_FeSiO3*(M_Mg*(3*M_O - M_SiO2 - M_m) + M_O*(-2*M_SiO2 + 5*M_m) - 2*M_SiO2*M_m) +
    M_MgSiO3*(M_FeO*(2*M_SiO2 - 2*M_m) + M_FeSiO3*(2*M_SiO2 - 2*M_m)))) + M_O*(M_Mg*(M_FeO*(M_FeSiO3*(-2*M_SiO2 + 2*M_m) +
    M_m*(2*M_MgSiO3 + 2*M_SiO2)) + 2*M_FeSiO3*M_SiO2*M_m) + M_MgSiO3*(M_FeO*(M_FeSiO3*(2*M_SiO2 - 2*M_m) - 2*M_SiO2*M_m) -
    2*M_FeSiO3*M_SiO2*M_m)) + M_Si*(M_Mg*(M_FeO*(M_FeSiO3*(-2*M_SiO2 + 2*M_m) + M_MgSiO3*(-3*M_O - M_SiO2 + 3*M_m) +
    M_O*(-M_SiO2 - 2*M_m) + 2*M_SiO2*M_m) + M_FeSiO3*(M_O*(2*M_SiO2 - 5*M_m) + 2*M_SiO2*M_m)) +
    M_MgSiO3*(M_FeO*(M_FeSiO3*(2*M_SiO2 - 2*M_m) + M_O*(M_SiO2 + 2*M_m) - 2*M_SiO2*M_m) + M_FeSiO3*(M_O*(-2*M_SiO2 + 5*M_m)
    - 2*M_SiO2*M_m)))) + M_FeSiO3*dKFeSiO3_KFeSiO3*(M_MgO*(M_O*(M_Fe*(M_MgSiO3*(2*M_SiO2 - 2*M_m) + M_SiO2*(2*M_FeO -
    2*M_m)) + M_Mg*M_SiO2*(-2*M_Fe + 2*M_FeO)) + M_Si*(M_Fe*(M_FeO*(-3*M_O + M_SiO2 + M_m) + M_MgSiO3*(2*M_SiO2 - 2*M_m) +
    M_O*(-2*M_SiO2 + 5*M_m) - 2*M_SiO2*M_m) + M_FeO*M_O*(3*M_SiO2 - 3*M_m) + M_Mg*(M_Fe*(3*M_O - M_SiO2 - M_m) +
    M_FeO*(-3*M_O + M_SiO2 + M_m)))) + M_O*(M_Fe*M_MgSiO3*M_SiO2*(2*M_FeO - 2*M_m) + M_Mg*(M_Fe*M_SiO2*(-2*M_FeO + 2*M_m) +
    2*M_FeO*M_MgSiO3*M_m)) + M_Si*(M_Mg*(M_Fe*(M_FeO*(3*M_O - M_SiO2 - M_m) + M_O*(2*M_SiO2 - 5*M_m) + 2*M_SiO2*M_m) +
    M_FeO*(M_MgSiO3*(-3*M_O - M_SiO2 + 3*M_m) + M_O*(-3*M_SiO2 + 3*M_m))) + M_MgSiO3*(M_Fe*(M_FeO*(-3*M_O + M_SiO2 + M_m) +
    M_O*(-2*M_SiO2 + 5*M_m) - 2*M_SiO2*M_m) + M_FeO*M_O*(3*M_SiO2 - 3*M_m)))) +
    M_Mg*dKMgO_KMgO*(M_MgO*(M_O*(M_FeO*(M_FeSiO3*(2*M_SiO2 - 2*M_m) - 2*M_SiO2*M_m) + M_MgSiO3*(M_Fe*(-2*M_SiO2 + 2*M_m) +
    M_FeO*(2*M_SiO2 - 2*M_m) + M_FeSiO3*(2*M_SiO2 - 2*M_m)) + M_m*(M_Fe*(2*M_FeSiO3 + 2*M_SiO2) - 2*M_FeSiO3*M_SiO2)) +
    M_Si*(M_Fe*(M_FeSiO3*(-3*M_O - M_SiO2 + 3*M_m) + M_O*(-M_SiO2 - 2*M_m) + 2*M_SiO2*M_m) + M_FeO*(M_FeSiO3*(2*M_SiO2 -
    2*M_m) + M_O*(M_SiO2 + 2*M_m) - 2*M_SiO2*M_m) + M_FeSiO3*(M_O*(M_SiO2 + 2*M_m) - 2*M_SiO2*M_m) +
    M_MgSiO3*(M_Fe*(-2*M_SiO2 + 2*M_m) + M_FeO*(2*M_SiO2 - 2*M_m) + M_FeSiO3*(2*M_SiO2 - 2*M_m)))) +
    M_MgSiO3*(M_O*(M_FeO*(M_FeSiO3*(2*M_SiO2 - 2*M_m) - 2*M_SiO2*M_m) + M_SiO2*(M_Fe*(-2*M_FeO + 2*M_m) - 2*M_FeSiO3*M_m)) +
    M_Si*(M_Fe*(M_FeO*(3*M_O - M_SiO2 - M_m) + M_O*(2*M_SiO2 - 5*M_m) + 2*M_SiO2*M_m) + M_FeO*(M_FeSiO3*(2*M_SiO2 - 2*M_m) +
    M_O*(-2*M_SiO2 + 5*M_m) - 2*M_SiO2*M_m) + M_FeSiO3*(M_O*(-2*M_SiO2 + 5*M_m) - 2*M_SiO2*M_m)))) +
    M_MgSiO3*dKMgSiO3_KMgSiO3*(M_Mg*(M_O*(M_FeO*(M_FeSiO3*(2*M_SiO2 - 2*M_m) - 2*M_SiO2*M_m) + M_SiO2*(M_Fe*(-2*M_FeO +
    2*M_m) - 2*M_FeSiO3*M_m)) + M_Si*(M_Fe*(M_FeO*(3*M_O - M_SiO2 - M_m) + M_O*(2*M_SiO2 - 5*M_m) + 2*M_SiO2*M_m) +
    M_FeO*(M_FeSiO3*(2*M_SiO2 - 2*M_m) + M_O*(-2*M_SiO2 + 5*M_m) - 2*M_SiO2*M_m) + M_FeSiO3*(M_O*(-2*M_SiO2 + 5*M_m) -
    2*M_SiO2*M_m))) + M_MgO*(M_O*(M_Fe*(2*M_FeO*M_SiO2 + 2*M_FeSiO3*M_m) + M_Mg*M_SiO2*(-2*M_Fe + 2*M_FeO + 2*M_FeSiO3)) +
    M_Si*(M_Fe*(M_FeO*(-3*M_O + M_SiO2 + M_m) + M_FeSiO3*(-3*M_O - M_SiO2 + 3*M_m) + M_O*(-3*M_SiO2 + 3*M_m)) +
    M_Mg*(M_Fe*(3*M_O - M_SiO2 - M_m) + M_FeO*(-3*M_O + M_SiO2 + M_m) + M_FeSiO3*(-3*M_O + M_SiO2 + M_m)) +
    M_O*(M_FeO*(3*M_SiO2 - 3*M_m) + M_FeSiO3*(3*M_SiO2 - 3*M_m))))) + M_Si*dKSiO2_KSiO2*(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(3*M_O
    + M_SiO2 - M_m) + M_SiO2*(M_O - M_m)) + M_m*(M_FeSiO3*(-5*M_O - M_SiO2) - 3*M_O*M_SiO2)) + M_MgSiO3*(M_Fe*(M_FeO*(3*M_O
    + M_SiO2 - M_m) + M_m*(-5*M_O - M_SiO2)) + M_FeO*(M_FeSiO3*(-3*M_O - M_SiO2 + M_m) + M_m*(5*M_O + M_SiO2)) +
    M_FeSiO3*M_m*(5*M_O + M_SiO2)) + M_O*(M_FeO*(M_FeSiO3*(-3*M_SiO2 + 3*M_m) + 3*M_SiO2*M_m) + 3*M_FeSiO3*M_SiO2*M_m)) +
    M_MgO*(M_Fe*(M_FeO*(M_FeSiO3*(-3*M_O - M_SiO2 + M_m) + M_SiO2*(-M_O + M_m)) + M_m*(M_FeSiO3*(5*M_O + M_SiO2) +
    3*M_O*M_SiO2)) + M_Mg*(M_Fe*(M_FeSiO3*(3*M_O + M_SiO2 - M_m) + M_SiO2*(M_O - M_m)) + M_FeO*(M_FeSiO3*(-3*M_O - M_SiO2 +
    M_m) + M_SiO2*(-M_O + M_m)) + M_FeSiO3*M_SiO2*(-M_O + M_m) + M_MgSiO3*(M_Fe*(3*M_O + M_SiO2 - M_m) + M_FeO*(-3*M_O -
    M_SiO2 + M_m) + M_FeSiO3*(-3*M_O - M_SiO2 + M_m))) + M_MgSiO3*(M_Fe*(M_FeO*(-3*M_O - M_SiO2 + M_m) + M_FeSiO3*(-3*M_O -
    M_SiO2 + M_m) + M_O*(-3*M_SiO2 + 3*M_m)) + M_O*(M_FeO*(3*M_SiO2 - 3*M_m) + M_FeSiO3*(3*M_SiO2 - 3*M_m))) +
    M_O*(M_FeO*(M_FeSiO3*(3*M_SiO2 - 3*M_m) - 3*M_SiO2*M_m) - 3*M_FeSiO3*M_SiO2*M_m)) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(-3*M_O - M_SiO2 + M_m) + M_SiO2*(-M_O + M_m)) + M_m*(M_FeSiO3*(5*M_O + M_SiO2) +
    3*M_O*M_SiO2)) + M_O*(M_FeO*(M_FeSiO3*(3*M_SiO2 - 3*M_m) - 3*M_SiO2*M_m) -
    3*M_FeSiO3*M_SiO2*M_m))))/(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) +
    M_SiO2*(-4*M_O*M_m + M_c*(M_O + M_m))) + M_m*(M_FeSiO3*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2)) - M_O*M_SiO2*M_c)) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_m*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2))) +
    M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) +
    M_FeSiO3*M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) + M_O*M_c*(M_FeO*(M_FeSiO3*(-M_SiO2 + M_m) + M_SiO2*M_m) +
    M_FeSiO3*M_SiO2*M_m)) + M_MgO*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) +
    M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_Mg*(M_Fe*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m + M_c*(M_O + M_m))) +
    M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) +
    M_FeSiO3*M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m)) + M_MgSiO3*(M_Fe*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) +
    M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)))) + M_MgSiO3*(M_Fe*(M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 +
    4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_O*M_c*(-M_SiO2 + M_m)) + M_O*M_c*(M_FeO*(M_SiO2 - M_m) + M_FeSiO3*(M_SiO2 -
    M_m))) + M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m) +
    M_Si*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)) + M_Mg*(M_Fe*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O +
    M_SiO2 - M_m)) + M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2
    + M_m)) + M_FeSiO3*(M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) + M_MgSiO3*(M_Fe*(-9*M_O + M_SiO2 +
    4*M_c - M_m) + M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m))) +
    M_MgSiO3*(M_Fe*(M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(9*M_SiO2 - 9*M_m)
    + M_c*(-4*M_SiO2 + 4*M_m)) + M_FeO*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + M_FeSiO3*(M_O*(-9*M_SiO2 +
    9*M_m) + M_c*(4*M_SiO2 - 4*M_m))))) + M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m)) +
    M_Si*(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O +
    M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - 9*M_m)) + 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(9*M_SiO2 - 9*M_m) + M_c*(-4*M_SiO2 + 4*M_m)) -
    9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeSiO3*(-9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) +
    4*M_SiO2*M_m)) + M_MgSiO3*(M_Fe*(M_FeO*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m +
    M_c*(-M_O + M_SiO2 - 9*M_m)) + M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m +
    M_c*(M_O - M_SiO2 + 9*M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)))) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)))))

def dM_c_dMi(dMi_c):
    '''alternate way to compute dM_c given the results of the mole changes of all core components 
    inputs:
    dMi_c: [dM_Mg, dM_Si, dM_Fe, dM_O]
    '''
    return np.sum(dMi_c)

def dM_MgO(Moles, dKs):
    M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = Moles
    dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs
    return M_MgO*(M_Fe*dKFeO_KFeO*(M_Mg*(M_FeO*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) +
    M_MgSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m + M_c*(M_O + M_m))) +
    M_FeSiO3*M_SiO2*M_m*(-4*M_O + M_c)) + M_MgSiO3*M_O*M_c*(-M_FeO*M_SiO2 - M_FeSiO3*M_m) +
    M_Si*(M_Mg*(M_FeO*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_MgSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 +
    4*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(2*M_SiO2 + 10*M_m) - M_SiO2*M_m + M_c*(-M_O - M_SiO2 -
    3*M_m))) + M_MgSiO3*(M_FeO*(M_O*(3*M_SiO2 + 6*M_m) + M_c*(-2*M_SiO2 - 2*M_m)) + M_FeSiO3*(M_O*(-6*M_SiO2 + 15*M_m) +
    M_c*(2*M_SiO2 - 6*M_m))))) + M_FeSiO3*dKFeSiO3_KFeSiO3*(M_Mg*(M_Fe*M_SiO2*M_m*(-4*M_O + M_c) +
    M_FeO*(M_MgSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_O*M_SiO2*M_c)) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) - M_O*M_c*M_m) + M_FeO*M_O*M_c*(-M_SiO2 +
    M_m)) + M_Si*(M_Mg*(M_Fe*(M_O*(2*M_SiO2 + 10*M_m) - M_SiO2*M_m + M_c*(-M_O - M_SiO2 - 3*M_m)) + M_FeO*(M_MgSiO3*(-9*M_O
    + M_SiO2 + 4*M_c - M_m) + M_O*(-3*M_SiO2 - 6*M_m) + M_c*(2*M_SiO2 + 2*M_m))) + M_MgSiO3*(M_Fe*(M_FeO*(-9*M_O + M_SiO2 +
    4*M_c - M_m) + M_O*(-6*M_SiO2 + 15*M_m) + M_c*(2*M_SiO2 - 6*M_m)) + M_FeO*(M_O*(9*M_SiO2 - 9*M_m) + M_c*(-4*M_SiO2 +
    4*M_m))))) + M_Mg*dKMgO_KMgO*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) +
    M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_O*M_c*M_m) + M_O*M_c*M_m*(-M_FeO -
    M_FeSiO3)) + M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m) +
    M_Si*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)) + M_MgSiO3*(M_Fe*(M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(6*M_SiO2 - 15*M_m) + M_c*(-2*M_SiO2 +
    6*M_m)) + M_FeO*(M_O*(-6*M_SiO2 + 15*M_m) + M_c*(2*M_SiO2 - 6*M_m)) + M_FeSiO3*(M_O*(-6*M_SiO2 + 15*M_m) + M_c*(2*M_SiO2
    - 6*M_m))))) + M_MgSiO3*dKMgSiO3_KMgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) +
    M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_Mg*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_SiO2*M_m*(4*M_O - M_c)) +
    M_SiO2*M_m*(M_Fe*(-4*M_O + M_c) + M_FeSiO3*(4*M_O - M_c))) + M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) -
    M_FeSiO3*M_SiO2*M_m) + M_Si*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m +
    M_c*(M_O - M_SiO2 + M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) -
    9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 -
    4*M_m)) + 9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2
    + M_m) - 4*M_SiO2*M_m)) + M_Mg*(M_Fe*(M_O*(2*M_SiO2 + 10*M_m) - M_SiO2*M_m + M_c*(-M_O - M_SiO2 - 3*M_m)) +
    M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(-2*M_SiO2 - 10*M_m) + M_SiO2*M_m + M_c*(M_O + M_SiO2 + 3*M_m)) +
    M_FeSiO3*(M_O*(-2*M_SiO2 - 10*M_m) + M_SiO2*M_m + M_c*(M_O + M_SiO2 + 3*M_m))))) +
    M_Si*dKSiO2_KSiO2*(M_Mg*(M_Fe*(M_FeSiO3*(M_O*(-4*M_SiO2 + 10*M_m) + M_c*(-M_O + M_SiO2 - 3*M_m)) + M_SiO2*(6*M_O*M_m +
    M_c*(-M_O - 2*M_m))) + M_FeO*(M_FeSiO3*(M_O*(6*M_SiO2 - 6*M_m) + M_c*(-2*M_SiO2 + 2*M_m)) + M_SiO2*(-6*M_O*M_m +
    M_c*(M_O + 2*M_m))) + M_FeSiO3*M_SiO2*(-6*M_O*M_m + M_c*(M_O + 2*M_m)) + M_MgSiO3*(M_Fe*(M_O*(-4*M_SiO2 + 10*M_m) +
    M_c*(-M_O + M_SiO2 - 3*M_m)) + M_FeO*(M_O*(4*M_SiO2 - 10*M_m) + M_c*(M_O - M_SiO2 + 3*M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 -
    10*M_m) + M_c*(M_O - M_SiO2 + 3*M_m)))) + M_MgSiO3*(M_Fe*(M_FeO*(M_O*(-2*M_SiO2 - 4*M_m) + M_c*(M_O + M_SiO2 + M_m)) +
    M_FeSiO3*(M_O*(4*M_SiO2 - 10*M_m) + M_c*(M_O - M_SiO2 + 3*M_m)) + M_O*M_c*(M_SiO2 - M_m)) + M_O*M_c*(M_FeO*(-M_SiO2 +
    M_m) + M_FeSiO3*(-M_SiO2 + M_m)))))/(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) +
    M_SiO2*(-4*M_O*M_m + M_c*(M_O + M_m))) + M_m*(M_FeSiO3*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2)) - M_O*M_SiO2*M_c)) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_m*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2))) +
    M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) +
    M_FeSiO3*M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) + M_O*M_c*(M_FeO*(M_FeSiO3*(-M_SiO2 + M_m) + M_SiO2*M_m) +
    M_FeSiO3*M_SiO2*M_m)) + M_MgO*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) +
    M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_Mg*(M_Fe*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m + M_c*(M_O + M_m))) +
    M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) +
    M_FeSiO3*M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m)) + M_MgSiO3*(M_Fe*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) +
    M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)))) + M_MgSiO3*(M_Fe*(M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 +
    4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_O*M_c*(-M_SiO2 + M_m)) + M_O*M_c*(M_FeO*(M_SiO2 - M_m) + M_FeSiO3*(M_SiO2 -
    M_m))) + M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m) +
    M_Si*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)) + M_Mg*(M_Fe*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O +
    M_SiO2 - M_m)) + M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2
    + M_m)) + M_FeSiO3*(M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) + M_MgSiO3*(M_Fe*(-9*M_O + M_SiO2 +
    4*M_c - M_m) + M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m))) +
    M_MgSiO3*(M_Fe*(M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(9*M_SiO2 - 9*M_m)
    + M_c*(-4*M_SiO2 + 4*M_m)) + M_FeO*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + M_FeSiO3*(M_O*(-9*M_SiO2 +
    9*M_m) + M_c*(4*M_SiO2 - 4*M_m))))) + M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m)) +
    M_Si*(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O +
    M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - 9*M_m)) + 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(9*M_SiO2 - 9*M_m) + M_c*(-4*M_SiO2 + 4*M_m)) -
    9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeSiO3*(-9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) +
    4*M_SiO2*M_m)) + M_MgSiO3*(M_Fe*(M_FeO*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m +
    M_c*(-M_O + M_SiO2 - 9*M_m)) + M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m +
    M_c*(M_O - M_SiO2 + 9*M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)))) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)))))

def dM_SiO2(Moles, dKs):
    M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = Moles
    dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs
    return M_SiO2*(M_Fe*dKFeO_KFeO*(M_FeSiO3*M_MgSiO3*M_O*M_c*(-M_FeO + M_m) + M_Mg*(M_FeO*(M_FeSiO3*M_O*M_c + M_MgSiO3*M_m*(4*M_O
    - M_c)) - M_FeSiO3*M_O*M_c*M_m) + M_MgO*(M_FeSiO3*(M_Mg*(-4*M_O*M_m + M_c*(M_O + M_m)) + M_O*M_c*(-M_FeO + M_m)) +
    M_MgSiO3*M_O*M_c*(-M_FeO - M_FeSiO3) + M_Si*(M_FeO*(M_FeSiO3*(9*M_O - 4*M_c) - 6*M_O*M_m + M_c*(M_O + 2*M_m)) +
    M_FeSiO3*(M_Mg*(-3*M_O + 2*M_c - M_m) - 15*M_O*M_m + M_c*(M_O + 6*M_m)) + M_MgSiO3*(M_FeO*(9*M_O - 4*M_c) +
    M_FeSiO3*(9*M_O - 4*M_c)))) + M_Si*(M_Mg*(M_FeO*(M_FeSiO3*(-9*M_O + 4*M_c) + M_MgSiO3*(-6*M_O + 2*M_c + M_m) + 6*M_O*M_m
    + M_c*(-M_O - 2*M_m)) + M_FeSiO3*(15*M_O*M_m + M_c*(-M_O - 6*M_m))) + M_MgSiO3*(M_FeO*(M_FeSiO3*(9*M_O - 4*M_c) -
    6*M_O*M_m + M_c*(M_O + 2*M_m)) + M_FeSiO3*(-15*M_O*M_m + M_c*(M_O + 6*M_m))))) +
    M_FeSiO3*dKFeSiO3_KFeSiO3*(M_Mg*(M_Fe*(M_FeO*(-4*M_O*M_m + M_c*(M_O + M_m)) - M_O*M_c*M_m) + M_FeO*M_m*(M_MgSiO3*(4*M_O
    - M_c) + M_O*M_c)) + M_MgO*(M_Fe*(M_FeO*(4*M_O*M_m + M_c*(-M_O - M_m)) + M_O*M_c*(-M_MgSiO3 + M_m)) - M_FeO*M_O*M_c*M_m
    + M_Mg*(M_Fe*(-4*M_O*M_m + M_c*(M_O + M_m)) + M_FeO*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_Si*(M_Fe*(M_FeO*(3*M_O - 2*M_c +
    M_m) + M_MgSiO3*(9*M_O - 4*M_c) - 15*M_O*M_m + M_c*(M_O + 6*M_m)) + M_FeO*M_m*(9*M_O - 4*M_c) + M_Mg*(M_Fe*(-3*M_O +
    2*M_c - M_m) + M_FeO*(3*M_O - 2*M_c + M_m)))) + M_MgSiO3*(M_Fe*(M_FeO*(4*M_O*M_m + M_c*(-M_O - M_m)) + M_O*M_c*M_m) -
    M_FeO*M_O*M_c*M_m) + M_Si*(M_Mg*(M_Fe*(M_FeO*(-3*M_O + 2*M_c - M_m) + 15*M_O*M_m + M_c*(-M_O - 6*M_m)) +
    M_FeO*(M_MgSiO3*(-6*M_O + 2*M_c + M_m) + M_m*(-9*M_O + 4*M_c))) + M_MgSiO3*(M_Fe*(M_FeO*(3*M_O - 2*M_c + M_m) -
    15*M_O*M_m + M_c*(M_O + 6*M_m)) + M_FeO*M_m*(9*M_O - 4*M_c)))) + M_Mg*dKMgO_KMgO*(M_MgO*(M_FeSiO3*(M_Fe*M_m*(4*M_O -
    M_c) - M_FeO*M_O*M_c) + M_MgSiO3*M_O*M_c*(M_Fe - M_FeO - M_FeSiO3) + M_Si*(M_Fe*(M_FeSiO3*(-6*M_O + 2*M_c + M_m) +
    6*M_O*M_m + M_c*(-M_O - 2*M_m)) + M_FeO*(M_FeSiO3*(9*M_O - 4*M_c) - 6*M_O*M_m + M_c*(M_O + 2*M_m)) +
    M_FeSiO3*(-6*M_O*M_m + M_c*(M_O + 2*M_m)) + M_MgSiO3*(M_Fe*(-9*M_O + 4*M_c) + M_FeO*(9*M_O - 4*M_c) + M_FeSiO3*(9*M_O -
    4*M_c)))) + M_MgSiO3*(M_Fe*(M_FeO*(-4*M_O*M_m + M_c*(M_O + M_m)) - M_O*M_c*M_m) + M_O*M_c*(M_FeO*(-M_FeSiO3 + M_m) +
    M_FeSiO3*M_m) + M_Si*(M_Fe*(M_FeO*(-3*M_O + 2*M_c - M_m) + 15*M_O*M_m + M_c*(-M_O - 6*M_m)) + M_FeO*(M_FeSiO3*(9*M_O -
    4*M_c) - 15*M_O*M_m + M_c*(M_O + 6*M_m)) + M_FeSiO3*(-15*M_O*M_m + M_c*(M_O + 6*M_m))))) +
    M_MgSiO3*dKMgSiO3_KMgSiO3*(M_Mg*(M_Fe*(M_FeO*(-4*M_O*M_m + M_c*(M_O + M_m)) - M_O*M_c*M_m) + M_O*M_c*(M_FeO*(-M_FeSiO3 +
    M_m) + M_FeSiO3*M_m) + M_Si*(M_Fe*(M_FeO*(-3*M_O + 2*M_c - M_m) + 15*M_O*M_m + M_c*(-M_O - 6*M_m)) +
    M_FeO*(M_FeSiO3*(9*M_O - 4*M_c) - 15*M_O*M_m + M_c*(M_O + 6*M_m)) + M_FeSiO3*(-15*M_O*M_m + M_c*(M_O + 6*M_m)))) +
    M_MgO*(M_Fe*(M_FeO*(4*M_O*M_m + M_c*(-M_O - M_m)) + M_m*(M_FeSiO3*(4*M_O - M_c) + M_O*M_c)) + M_Mg*(M_Fe*(-4*M_O*M_m +
    M_c*(M_O + M_m)) + M_FeO*(4*M_O*M_m + M_c*(-M_O - M_m)) + M_FeSiO3*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_O*M_c*M_m*(-M_FeO
    - M_FeSiO3) + M_Si*(M_Fe*(M_FeO*(3*M_O - 2*M_c + M_m) + M_FeSiO3*(-6*M_O + 2*M_c + M_m) + M_m*(-9*M_O + 4*M_c)) +
    M_Mg*(M_Fe*(-3*M_O + 2*M_c - M_m) + M_FeO*(3*M_O - 2*M_c + M_m) + M_FeSiO3*(3*M_O - 2*M_c + M_m)) + M_m*(M_FeO*(9*M_O -
    4*M_c) + M_FeSiO3*(9*M_O - 4*M_c))))) + M_Si*dKSiO2_KSiO2*(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(6*M_O - 2*M_c) - 4*M_O*M_m +
    M_c*(M_O + M_m)) + M_m*(M_FeSiO3*(-10*M_O + 3*M_c) - M_O*M_c)) + M_MgSiO3*(M_Fe*(M_FeO*(6*M_O - 2*M_c) + M_m*(-10*M_O +
    3*M_c)) + M_FeO*(M_FeSiO3*(-6*M_O + 2*M_c) + M_m*(10*M_O - 3*M_c)) + M_FeSiO3*M_m*(10*M_O - 3*M_c)) + M_O*M_c*M_m*(M_FeO
    + M_FeSiO3)) + M_MgO*(M_Fe*(M_FeO*(M_FeSiO3*(-6*M_O + 2*M_c) + 4*M_O*M_m + M_c*(-M_O - M_m)) + M_MgSiO3*(M_FeO*(-6*M_O +
    2*M_c) + M_FeSiO3*(-6*M_O + 2*M_c)) + M_m*(M_FeSiO3*(10*M_O - 3*M_c) + M_O*M_c)) + M_Mg*(M_Fe*(M_FeSiO3*(6*M_O - 2*M_c)
    - 4*M_O*M_m + M_c*(M_O + M_m)) + M_FeO*(M_FeSiO3*(-6*M_O + 2*M_c) + 4*M_O*M_m + M_c*(-M_O - M_m)) + M_FeSiO3*(4*M_O*M_m
    + M_c*(-M_O - M_m)) + M_MgSiO3*(M_Fe*(6*M_O - 2*M_c) + M_FeO*(-6*M_O + 2*M_c) + M_FeSiO3*(-6*M_O + 2*M_c))) +
    M_O*M_c*M_m*(-M_FeO - M_FeSiO3)) + M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(-6*M_O + 2*M_c) + 4*M_O*M_m + M_c*(-M_O - M_m)) +
    M_m*(M_FeSiO3*(10*M_O - 3*M_c) + M_O*M_c)) + M_O*M_c*M_m*(-M_FeO -
    M_FeSiO3))))/(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m +
    M_c*(M_O + M_m))) + M_m*(M_FeSiO3*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2)) - M_O*M_SiO2*M_c)) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_m*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2))) +
    M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) +
    M_FeSiO3*M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) + M_O*M_c*(M_FeO*(M_FeSiO3*(-M_SiO2 + M_m) + M_SiO2*M_m) +
    M_FeSiO3*M_SiO2*M_m)) + M_MgO*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) +
    M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_Mg*(M_Fe*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m + M_c*(M_O + M_m))) +
    M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) +
    M_FeSiO3*M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m)) + M_MgSiO3*(M_Fe*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) +
    M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)))) + M_MgSiO3*(M_Fe*(M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 +
    4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_O*M_c*(-M_SiO2 + M_m)) + M_O*M_c*(M_FeO*(M_SiO2 - M_m) + M_FeSiO3*(M_SiO2 -
    M_m))) + M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m) +
    M_Si*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)) + M_Mg*(M_Fe*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O +
    M_SiO2 - M_m)) + M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2
    + M_m)) + M_FeSiO3*(M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) + M_MgSiO3*(M_Fe*(-9*M_O + M_SiO2 +
    4*M_c - M_m) + M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m))) +
    M_MgSiO3*(M_Fe*(M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(9*M_SiO2 - 9*M_m)
    + M_c*(-4*M_SiO2 + 4*M_m)) + M_FeO*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + M_FeSiO3*(M_O*(-9*M_SiO2 +
    9*M_m) + M_c*(4*M_SiO2 - 4*M_m))))) + M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m)) +
    M_Si*(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O +
    M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - 9*M_m)) + 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(9*M_SiO2 - 9*M_m) + M_c*(-4*M_SiO2 + 4*M_m)) -
    9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeSiO3*(-9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) +
    4*M_SiO2*M_m)) + M_MgSiO3*(M_Fe*(M_FeO*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m +
    M_c*(-M_O + M_SiO2 - 9*M_m)) + M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m +
    M_c*(M_O - M_SiO2 + 9*M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)))) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)))))

def dM_FeO(Moles, dKs):
    M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = Moles
    dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs
    return M_FeO*(M_Fe*dKFeO_KFeO*(M_MgO*(M_Mg*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) +
    M_MgSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) +
    M_O*M_c*(M_MgSiO3*(M_SiO2 - M_m) + M_m*(-M_FeSiO3 - M_SiO2)) + M_Si*(M_FeSiO3*(M_O*(-6*M_SiO2 + 15*M_m) + M_c*(2*M_SiO2
    - 6*M_m)) + M_Mg*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_MgSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 -
    4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) + M_MgSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m))) + M_Si*(M_Mg*(M_FeSiO3*(M_O*(6*M_SiO2 - 15*M_m) +
    M_c*(-2*M_SiO2 + 6*M_m)) + M_MgSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) -
    9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_MgSiO3*(M_FeSiO3*(M_O*(-6*M_SiO2 + 15*M_m) +
    M_c*(2*M_SiO2 - 6*M_m)) + 9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m))) +
    M_m*(M_Mg*(M_MgSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_c*(M_FeSiO3 + M_SiO2)) + M_MgSiO3*M_O*M_c*(-M_FeSiO3 -
    M_SiO2))) + M_FeSiO3*dKFeSiO3_KFeSiO3*(M_MgO*(M_Mg*(M_MgSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) +
    M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_MgSiO3*(M_Fe*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) +
    M_O*M_c*(M_SiO2 - M_m)) + M_Si*(M_Fe*(M_O*(-2*M_SiO2 - 10*M_m) + M_SiO2*M_m + M_c*(M_O + M_SiO2 + 3*M_m)) +
    M_Mg*(M_MgSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) +
    M_MgSiO3*(M_Fe*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_SiO2*M_m*(M_Fe*(4*M_O - M_c) - M_O*M_c)) + M_Si*(M_Mg*(M_Fe*(M_O*(2*M_SiO2
    + 10*M_m) - M_SiO2*M_m + M_c*(-M_O - M_SiO2 - 3*M_m)) + M_MgSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O -
    M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_MgSiO3*(M_Fe*(M_O*(-2*M_SiO2 -
    10*M_m) + M_SiO2*M_m + M_c*(M_O + M_SiO2 + 3*M_m)) + 9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m))) +
    M_m*(M_Mg*(M_MgSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_SiO2*(M_Fe*(-4*M_O + M_c) + M_O*M_c)) +
    M_MgSiO3*M_SiO2*(M_Fe*(4*M_O - M_c) - M_O*M_c))) + M_Mg*dKMgO_KMgO*(M_MgO*(M_Fe*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) +
    M_c*(M_O - M_SiO2 + M_m)) + M_MgSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m + M_c*(M_O
    + M_m))) - M_FeSiO3*M_O*M_SiO2*M_c + M_Si*(M_Fe*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_MgSiO3*(-9*M_O + M_SiO2 +
    4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(3*M_SiO2 + 6*M_m) +
    M_c*(-2*M_SiO2 - 2*M_m)))) + M_MgSiO3*(M_Si*(M_Fe*(M_O*(2*M_SiO2 + 10*M_m) - M_SiO2*M_m + M_c*(-M_O - M_SiO2 - 3*M_m)) +
    M_FeSiO3*(M_O*(-6*M_SiO2 + 15*M_m) + M_c*(2*M_SiO2 - 6*M_m))) + M_m*(M_Fe*M_SiO2*(-4*M_O + M_c) - M_FeSiO3*M_O*M_c))) +
    M_MgSiO3*dKMgSiO3_KMgSiO3*(M_Mg*(M_Si*(M_Fe*(M_O*(2*M_SiO2 + 10*M_m) - M_SiO2*M_m + M_c*(-M_O - M_SiO2 - 3*M_m)) +
    M_FeSiO3*(M_O*(-6*M_SiO2 + 15*M_m) + M_c*(2*M_SiO2 - 6*M_m))) + M_m*(M_Fe*M_SiO2*(-4*M_O + M_c) - M_FeSiO3*M_O*M_c)) +
    M_MgO*(M_Fe*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_O*M_SiO2*M_c) +
    M_FeSiO3*(M_Mg*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_O*M_c*(-M_SiO2 + M_m)) +
    M_Si*(M_Fe*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-3*M_SiO2 - 6*M_m) + M_c*(2*M_SiO2 + 2*M_m)) +
    M_FeSiO3*(M_Mg*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(9*M_SiO2 - 9*M_m) + M_c*(-4*M_SiO2 + 4*M_m))))) +
    M_Si*dKSiO2_KSiO2*(M_Mg*(M_Fe*(M_FeSiO3*(M_O*(-4*M_SiO2 + 10*M_m) + M_c*(-M_O + M_SiO2 - 3*M_m)) + M_SiO2*(6*M_O*M_m +
    M_c*(-M_O - 2*M_m))) + M_FeSiO3*M_O*M_c*(M_SiO2 - M_m) + M_MgSiO3*(M_Fe*(M_O*(-4*M_SiO2 + 10*M_m) + M_c*(-M_O + M_SiO2 -
    3*M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 10*M_m) + M_c*(M_O - M_SiO2 + 3*M_m)))) + M_MgO*(M_Fe*(M_FeSiO3*(M_O*(4*M_SiO2 -
    10*M_m) + M_c*(M_O - M_SiO2 + 3*M_m)) + M_MgSiO3*(M_O*(6*M_SiO2 - 6*M_m) + M_c*(-2*M_SiO2 + 2*M_m)) + M_SiO2*(-6*M_O*M_m
    + M_c*(M_O + 2*M_m))) + M_FeSiO3*(M_Mg*(M_O*(-2*M_SiO2 - 4*M_m) + M_c*(M_O + M_SiO2 + M_m)) + M_O*M_c*(-M_SiO2 + M_m)))
    + M_MgSiO3*(M_Fe*(M_FeSiO3*(M_O*(4*M_SiO2 - 10*M_m) + M_c*(M_O - M_SiO2 + 3*M_m)) + M_SiO2*(-6*M_O*M_m + M_c*(M_O +
    2*M_m))) + M_FeSiO3*M_O*M_c*(-M_SiO2 + M_m))))/(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2
    + M_m)) + M_SiO2*(-4*M_O*M_m + M_c*(M_O + M_m))) + M_m*(M_FeSiO3*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2)) -
    M_O*M_SiO2*M_c)) + M_MgSiO3*(M_Fe*(M_FeO*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_m*(-4*M_O*M_SiO2 +
    M_c*(-M_O + M_SiO2))) + M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_m*(4*M_O*M_SiO2 +
    M_c*(M_O - M_SiO2))) + M_FeSiO3*M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) + M_O*M_c*(M_FeO*(M_FeSiO3*(-M_SiO2 + M_m) +
    M_SiO2*M_m) + M_FeSiO3*M_SiO2*M_m)) + M_MgO*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m))
    + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_Mg*(M_Fe*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m + M_c*(M_O + M_m))) +
    M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) +
    M_FeSiO3*M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m)) + M_MgSiO3*(M_Fe*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) +
    M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)))) + M_MgSiO3*(M_Fe*(M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 +
    4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_O*M_c*(-M_SiO2 + M_m)) + M_O*M_c*(M_FeO*(M_SiO2 - M_m) + M_FeSiO3*(M_SiO2 -
    M_m))) + M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m) +
    M_Si*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)) + M_Mg*(M_Fe*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O +
    M_SiO2 - M_m)) + M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2
    + M_m)) + M_FeSiO3*(M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) + M_MgSiO3*(M_Fe*(-9*M_O + M_SiO2 +
    4*M_c - M_m) + M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m))) +
    M_MgSiO3*(M_Fe*(M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(9*M_SiO2 - 9*M_m)
    + M_c*(-4*M_SiO2 + 4*M_m)) + M_FeO*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + M_FeSiO3*(M_O*(-9*M_SiO2 +
    9*M_m) + M_c*(4*M_SiO2 - 4*M_m))))) + M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m)) +
    M_Si*(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O +
    M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - 9*M_m)) + 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(9*M_SiO2 - 9*M_m) + M_c*(-4*M_SiO2 + 4*M_m)) -
    9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeSiO3*(-9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) +
    4*M_SiO2*M_m)) + M_MgSiO3*(M_Fe*(M_FeO*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m +
    M_c*(-M_O + M_SiO2 - 9*M_m)) + M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m +
    M_c*(M_O - M_SiO2 + 9*M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)))) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)))))

def dM_MgSiO3(Moles, dKs):
    M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = Moles
    dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs
    return M_MgSiO3*(M_Fe*dKFeO_KFeO*(M_Mg*(M_FeO*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) +
    M_SiO2*M_m*(-4*M_O + M_c)) + M_FeSiO3*M_m*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2)) + M_Si*(M_FeO*(M_FeSiO3*(-9*M_O + M_SiO2
    + 4*M_c - M_m) + M_O*(2*M_SiO2 + 10*M_m) - M_SiO2*M_m + M_c*(-M_O - M_SiO2 - 3*M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 +
    25*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - 9*M_m)))) + M_MgO*(M_FeO*M_O*M_SiO2*M_c + M_FeSiO3*(M_Mg*(M_O*(4*M_SiO2 -
    4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_O*M_c*M_m) + M_Si*(M_FeO*(M_O*(-3*M_SiO2 - 6*M_m) + M_c*(2*M_SiO2 + 2*M_m)) +
    M_FeSiO3*(M_Mg*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(6*M_SiO2 - 15*M_m) + M_c*(-2*M_SiO2 + 6*M_m))))) +
    M_FeSiO3*dKFeSiO3_KFeSiO3*(M_Mg*(M_Fe*(M_FeO*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_m*(-4*M_O*M_SiO2 +
    M_c*(-M_O + M_SiO2))) + M_FeO*M_O*M_c*M_m + M_Si*(M_Fe*(M_FeO*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-4*M_SiO2 + 25*M_m)
    - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - 9*M_m)) + M_FeO*(M_O*(6*M_SiO2 - 15*M_m) + M_c*(-2*M_SiO2 + 6*M_m)))) +
    M_MgO*(M_Fe*(M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_O*M_c*M_m) + M_FeO*M_O*M_c*(M_SiO2 - M_m) +
    M_Mg*(M_Fe*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m))) + M_Si*(M_Fe*(M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(6*M_SiO2 - 15*M_m) + M_c*(-2*M_SiO2 + 6*M_m)) +
    M_FeO*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + M_Mg*(M_Fe*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_FeO*(9*M_O -
    M_SiO2 - 4*M_c + M_m))))) + M_Mg*dKMgO_KMgO*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m))
    + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_MgO*(M_Fe*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) - M_O*M_SiO2*M_c) + M_O*M_SiO2*M_c*(M_FeO +
    M_FeSiO3) + M_Si*(M_Fe*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(3*M_SiO2 + 6*M_m) + M_c*(-2*M_SiO2 - 2*M_m)) +
    M_FeO*(M_O*(-3*M_SiO2 - 6*M_m) + M_c*(2*M_SiO2 + 2*M_m)) + M_FeSiO3*(M_O*(-3*M_SiO2 - 6*M_m) + M_c*(2*M_SiO2 + 2*M_m))))
    + M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m) + M_Si*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O -
    M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 -
    25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) +
    M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + 9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)))) +
    M_Si*dKSiO2_KSiO2*(M_Mg*(M_Fe*(M_FeO*(M_O*(-2*M_SiO2 - 4*M_m) + M_c*(M_O + M_SiO2 + M_m)) + M_m*(6*M_O*M_SiO2 +
    M_c*(-M_O - 2*M_SiO2))) + M_FeO*(M_FeSiO3*(M_O*(6*M_SiO2 - 6*M_m) + M_c*(-2*M_SiO2 + 2*M_m)) + M_m*(-6*M_O*M_SiO2 +
    M_c*(M_O + 2*M_SiO2))) + M_FeSiO3*M_m*(-6*M_O*M_SiO2 + M_c*(M_O + 2*M_SiO2))) + M_MgO*(M_Fe*(M_FeO*(M_O*(2*M_SiO2 +
    4*M_m) + M_c*(-M_O - M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 10*M_m) + M_c*(-M_O + M_SiO2 - 3*M_m)) +
    M_O*M_c*(-M_SiO2 + M_m)) + M_Mg*(M_Fe*(M_O*(-2*M_SiO2 - 4*M_m) + M_c*(M_O + M_SiO2 + M_m)) + M_FeO*(M_O*(2*M_SiO2 +
    4*M_m) + M_c*(-M_O - M_SiO2 - M_m)) + M_FeSiO3*(M_O*(2*M_SiO2 + 4*M_m) + M_c*(-M_O - M_SiO2 - M_m))) +
    M_O*M_c*(M_FeO*(M_SiO2 - M_m) + M_FeSiO3*(M_SiO2 - M_m)))) +
    dKMgSiO3_KMgSiO3*(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_SiO2*(4*M_O*M_m
    + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m) + M_Si*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O -
    M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 -
    25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) +
    M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + 9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)))) +
    M_MgO*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m + M_c*(M_O +
    M_m))) + M_m*(M_FeSiO3*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2)) - M_O*M_SiO2*M_c)) + M_Mg*(M_Fe*(M_FeSiO3*(M_O*(-4*M_SiO2 +
    4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_FeO*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) +
    M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m + M_c*(M_O + M_m))) + M_FeSiO3*M_SiO2*(-4*M_O*M_m + M_c*(M_O + M_m))) +
    M_O*M_c*(M_FeO*(M_FeSiO3*(-M_SiO2 + M_m) + M_SiO2*M_m) + M_FeSiO3*M_SiO2*M_m) + M_Si*(M_Fe*(M_FeO*(M_FeSiO3*(-9*M_O +
    M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 +
    25*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - 9*M_m)) + 9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) +
    M_FeO*(M_FeSiO3*(M_O*(9*M_SiO2 - 9*M_m) + M_c*(-4*M_SiO2 + 4*M_m)) - 9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) +
    4*M_SiO2*M_m)) + M_FeSiO3*(-9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_Mg*(M_Fe*(M_FeSiO3*(9*M_O -
    M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) + M_FeO*(M_FeSiO3*(-9*M_O + M_SiO2
    + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-M_SiO2 + 4*M_m) -
    M_SiO2*M_m + M_c*(-M_O + M_SiO2 - M_m)))))))/(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 +
    M_m)) + M_SiO2*(-4*M_O*M_m + M_c*(M_O + M_m))) + M_m*(M_FeSiO3*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2)) - M_O*M_SiO2*M_c))
    + M_MgSiO3*(M_Fe*(M_FeO*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_m*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2)))
    + M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) +
    M_FeSiO3*M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) + M_O*M_c*(M_FeO*(M_FeSiO3*(-M_SiO2 + M_m) + M_SiO2*M_m) +
    M_FeSiO3*M_SiO2*M_m)) + M_MgO*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) +
    M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_Mg*(M_Fe*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m + M_c*(M_O + M_m))) +
    M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) +
    M_FeSiO3*M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m)) + M_MgSiO3*(M_Fe*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) +
    M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)))) + M_MgSiO3*(M_Fe*(M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 +
    4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_O*M_c*(-M_SiO2 + M_m)) + M_O*M_c*(M_FeO*(M_SiO2 - M_m) + M_FeSiO3*(M_SiO2 -
    M_m))) + M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m) +
    M_Si*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)) + M_Mg*(M_Fe*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O +
    M_SiO2 - M_m)) + M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2
    + M_m)) + M_FeSiO3*(M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) + M_MgSiO3*(M_Fe*(-9*M_O + M_SiO2 +
    4*M_c - M_m) + M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m))) +
    M_MgSiO3*(M_Fe*(M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(9*M_SiO2 - 9*M_m)
    + M_c*(-4*M_SiO2 + 4*M_m)) + M_FeO*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + M_FeSiO3*(M_O*(-9*M_SiO2 +
    9*M_m) + M_c*(4*M_SiO2 - 4*M_m))))) + M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m)) +
    M_Si*(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O +
    M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - 9*M_m)) + 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(9*M_SiO2 - 9*M_m) + M_c*(-4*M_SiO2 + 4*M_m)) -
    9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeSiO3*(-9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) +
    4*M_SiO2*M_m)) + M_MgSiO3*(M_Fe*(M_FeO*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m +
    M_c*(-M_O + M_SiO2 - 9*M_m)) + M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m +
    M_c*(M_O - M_SiO2 + 9*M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)))) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)))))

def dM_FeSiO3(Moles, dKs):
    M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = Moles
    dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs
    return M_FeSiO3*(M_Fe*dKFeO_KFeO*(M_Mg*(M_MgSiO3*(M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) +
    M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) + M_O*M_SiO2*M_c*(-M_FeO + M_m)) + M_MgO*(M_Mg*(M_MgSiO3*(M_O*(-4*M_SiO2 +
    4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_O*M_c*(M_MgSiO3*(M_SiO2 - M_m) +
    M_SiO2*(M_FeO - M_m)) + M_Si*(M_FeO*(M_O*(-3*M_SiO2 - 6*M_m) + M_c*(2*M_SiO2 + 2*M_m)) + M_Mg*(M_MgSiO3*(9*M_O - M_SiO2
    - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) + M_MgSiO3*(M_O*(-9*M_SiO2 + 9*M_m) +
    M_c*(4*M_SiO2 - 4*M_m)) + 9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m))) + M_MgSiO3*M_O*M_SiO2*M_c*(M_FeO
    - M_m) + M_Si*(M_Mg*(M_FeO*(M_O*(3*M_SiO2 + 6*M_m) + M_c*(-2*M_SiO2 - 2*M_m)) + M_MgSiO3*(M_FeO*(9*M_O - M_SiO2 - 4*M_c
    + M_m) + M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m)
    + 4*M_SiO2*M_m)) + M_MgSiO3*(M_FeO*(M_O*(-3*M_SiO2 - 6*M_m) + M_c*(2*M_SiO2 + 2*M_m)) + 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)))) + M_Mg*dKMgO_KMgO*(M_MgO*(M_Fe*(M_MgSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O
    - M_SiO2 + M_m)) + M_SiO2*M_m*(-4*M_O + M_c)) + M_FeO*M_O*M_SiO2*M_c + M_Si*(M_Fe*(M_MgSiO3*(-9*M_O + M_SiO2 + 4*M_c -
    M_m) + M_O*(2*M_SiO2 + 10*M_m) - M_SiO2*M_m + M_c*(-M_O - M_SiO2 - 3*M_m)) + M_FeO*(M_O*(-3*M_SiO2 - 6*M_m) +
    M_c*(2*M_SiO2 + 2*M_m)))) + M_MgSiO3*(M_Fe*(M_FeO*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) +
    M_m*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2))) + M_FeO*M_O*M_c*M_m + M_Si*(M_Fe*(M_FeO*(-9*M_O + M_SiO2 + 4*M_c - M_m) +
    M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - 9*M_m)) + M_FeO*(M_O*(6*M_SiO2 - 15*M_m) + M_c*(-2*M_SiO2 +
    6*M_m))))) + M_MgSiO3*dKMgSiO3_KMgSiO3*(M_Mg*(M_Fe*(M_FeO*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) +
    M_m*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2))) + M_FeO*M_O*M_c*M_m + M_Si*(M_Fe*(M_FeO*(-9*M_O + M_SiO2 + 4*M_c - M_m) +
    M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - 9*M_m)) + M_FeO*(M_O*(6*M_SiO2 - 15*M_m) + M_c*(-2*M_SiO2 +
    6*M_m)))) + M_MgO*(M_Fe*(M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_O*M_c*M_m) +
    M_FeO*M_O*M_c*(M_SiO2 - M_m) + M_Mg*(M_Fe*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_FeO*(M_O*(-4*M_SiO2 +
    4*M_m) + M_c*(-M_O + M_SiO2 - M_m))) + M_Si*(M_Fe*(M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(6*M_SiO2 - 15*M_m) +
    M_c*(-2*M_SiO2 + 6*M_m)) + M_FeO*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + M_Mg*(M_Fe*(-9*M_O + M_SiO2 +
    4*M_c - M_m) + M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m))))) + M_Si*dKSiO2_KSiO2*(M_Mg*(M_Fe*(M_FeO*(M_O*(-2*M_SiO2 - 4*M_m)
    + M_c*(M_O + M_SiO2 + M_m)) + M_m*(6*M_O*M_SiO2 + M_c*(-M_O - 2*M_SiO2))) + M_FeO*(M_MgSiO3*(M_O*(-4*M_SiO2 + 10*M_m) +
    M_c*(-M_O + M_SiO2 - 3*M_m)) + M_O*M_c*(-M_SiO2 + M_m))) + M_MgO*(M_Fe*(M_FeO*(M_O*(2*M_SiO2 + 4*M_m) + M_c*(-M_O -
    M_SiO2 - M_m)) + M_MgSiO3*(M_O*(6*M_SiO2 - 6*M_m) + M_c*(-2*M_SiO2 + 2*M_m)) + M_m*(-6*M_O*M_SiO2 + M_c*(M_O +
    2*M_SiO2))) + M_FeO*M_O*M_c*(M_SiO2 - M_m) + M_Mg*(M_Fe*(M_O*(-2*M_SiO2 - 4*M_m) + M_c*(M_O + M_SiO2 + M_m)) +
    M_FeO*(M_O*(2*M_SiO2 + 4*M_m) + M_c*(-M_O - M_SiO2 - M_m)))) + M_MgSiO3*(M_Fe*(M_FeO*(M_O*(2*M_SiO2 + 4*M_m) + M_c*(-M_O
    - M_SiO2 - M_m)) + M_m*(-6*M_O*M_SiO2 + M_c*(M_O + 2*M_SiO2))) + M_FeO*M_O*M_c*(M_SiO2 - M_m))) +
    dKFeSiO3_KFeSiO3*(M_Mg*(M_MgSiO3*(M_Fe*(M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_m*(4*M_O*M_SiO2
    + M_c*(M_O - M_SiO2))) + M_FeO*M_m*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2))) + M_SiO2*(M_Fe*(M_FeO*(4*M_O*M_m + M_c*(-M_O -
    M_m)) + M_O*M_c*M_m) - M_FeO*M_O*M_c*M_m)) + M_MgO*(M_Mg*(M_MgSiO3*(M_Fe*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)) + M_FeO*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m))) + M_SiO2*(M_Fe*(4*M_O*M_m + M_c*(-M_O - M_m)) +
    M_FeO*(-4*M_O*M_m + M_c*(M_O + M_m)))) + M_MgSiO3*(M_Fe*(M_FeO*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) +
    M_O*M_c*(M_SiO2 - M_m)) + M_FeO*M_O*M_c*(-M_SiO2 + M_m)) + M_Si*(M_Fe*(M_FeO*(M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m +
    M_c*(-M_O + M_SiO2 - M_m)) + 9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeO*(-9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_Mg*(M_Fe*(M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) +
    M_FeO*(M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - M_m)) + M_MgSiO3*(M_Fe*(9*M_O - M_SiO2 - 4*M_c + M_m) +
    M_FeO*(-9*M_O + M_SiO2 + 4*M_c - M_m))) + M_MgSiO3*(M_Fe*(M_FeO*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-9*M_SiO2 +
    9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + M_FeO*(M_O*(9*M_SiO2 - 9*M_m) + M_c*(-4*M_SiO2 + 4*M_m)))) +
    M_SiO2*(M_Fe*(M_FeO*(-4*M_O*M_m + M_c*(M_O + M_m)) - M_O*M_c*M_m) + M_FeO*M_O*M_c*M_m)) +
    M_MgSiO3*M_SiO2*(M_Fe*(M_FeO*(-4*M_O*M_m + M_c*(M_O + M_m)) - M_O*M_c*M_m) + M_FeO*M_O*M_c*M_m) +
    M_Si*(M_Mg*(M_Fe*(M_FeO*(M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) +
    M_MgSiO3*(M_Fe*(M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    9*M_m)) + M_FeO*(M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - 9*M_m)))) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - M_m)) + 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeO*(-9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) +
    4*M_SiO2*M_m))))))/(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m
    + M_c*(M_O + M_m))) + M_m*(M_FeSiO3*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2)) - M_O*M_SiO2*M_c)) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_m*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2))) +
    M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) +
    M_FeSiO3*M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) + M_O*M_c*(M_FeO*(M_FeSiO3*(-M_SiO2 + M_m) + M_SiO2*M_m) +
    M_FeSiO3*M_SiO2*M_m)) + M_MgO*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) +
    M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_Mg*(M_Fe*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m + M_c*(M_O + M_m))) +
    M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) +
    M_FeSiO3*M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m)) + M_MgSiO3*(M_Fe*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) +
    M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)))) + M_MgSiO3*(M_Fe*(M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 +
    4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_O*M_c*(-M_SiO2 + M_m)) + M_O*M_c*(M_FeO*(M_SiO2 - M_m) + M_FeSiO3*(M_SiO2 -
    M_m))) + M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m) +
    M_Si*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)) + M_Mg*(M_Fe*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O +
    M_SiO2 - M_m)) + M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2
    + M_m)) + M_FeSiO3*(M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) + M_MgSiO3*(M_Fe*(-9*M_O + M_SiO2 +
    4*M_c - M_m) + M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m))) +
    M_MgSiO3*(M_Fe*(M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(9*M_SiO2 - 9*M_m)
    + M_c*(-4*M_SiO2 + 4*M_m)) + M_FeO*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + M_FeSiO3*(M_O*(-9*M_SiO2 +
    9*M_m) + M_c*(4*M_SiO2 - 4*M_m))))) + M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m)) +
    M_Si*(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O +
    M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - 9*M_m)) + 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(9*M_SiO2 - 9*M_m) + M_c*(-4*M_SiO2 + 4*M_m)) -
    9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeSiO3*(-9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) +
    4*M_SiO2*M_m)) + M_MgSiO3*(M_Fe*(M_FeO*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m +
    M_c*(-M_O + M_SiO2 - 9*M_m)) + M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m +
    M_c*(M_O - M_SiO2 + 9*M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)))) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)))))

def dM_m(Moles, dKs):
    M_Mg, M_Si, M_Fe, M_O, M_c, M_MgO, M_SiO2, M_FeO, M_MgSiO3, M_FeSiO3, M_m = Moles
    dKMgO_KMgO, dKSiO2_KSiO2, dKFeO_KFeO, dKMgSiO3_KMgSiO3, dKFeSiO3_KFeSiO3 = dKs
    return M_m*(M_Fe*dKFeO_KFeO*(M_FeO*(M_Mg*(M_MgSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_c*(M_FeSiO3 + M_SiO2)) +
    M_MgSiO3*M_O*M_c*(-M_FeSiO3 - M_SiO2)) + M_MgO*(M_FeSiO3*M_Mg*M_SiO2*(-4*M_O + M_c) + M_O*M_c*(M_FeO*(-M_FeSiO3 -
    M_SiO2) + M_MgSiO3*(-M_FeO - M_FeSiO3)) + M_Si*(M_FeO*(M_FeSiO3*(9*M_O - 4*M_c) + 3*M_O*M_SiO2 + M_c*(M_O - 2*M_SiO2)) +
    M_FeSiO3*(M_Mg*(6*M_O - M_SiO2 - 2*M_c) - 6*M_O*M_SiO2 + M_c*(M_O + 2*M_SiO2)) + M_MgSiO3*(M_FeO*(9*M_O - 4*M_c) +
    M_FeSiO3*(9*M_O - 4*M_c)))) + M_Si*(M_Mg*(M_FeO*(M_FeSiO3*(-9*M_O + 4*M_c) + M_MgSiO3*(-15*M_O + M_SiO2 + 6*M_c) -
    3*M_O*M_SiO2 + M_c*(-M_O + 2*M_SiO2)) + M_FeSiO3*(6*M_O*M_SiO2 + M_c*(-M_O - 2*M_SiO2))) +
    M_MgSiO3*(M_FeO*(M_FeSiO3*(9*M_O - 4*M_c) + 3*M_O*M_SiO2 + M_c*(M_O - 2*M_SiO2)) + M_FeSiO3*(-6*M_O*M_SiO2 + M_c*(M_O +
    2*M_SiO2))))) + M_FeSiO3*dKFeSiO3_KFeSiO3*(M_FeO*(M_Mg*(M_MgSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) +
    M_SiO2*(M_Fe*(-4*M_O + M_c) + M_O*M_c)) + M_MgSiO3*M_SiO2*(M_Fe*(4*M_O - M_c) - M_O*M_c)) +
    M_MgO*(M_Fe*(M_FeO*M_SiO2*(4*M_O - M_c) - M_MgSiO3*M_O*M_c) + M_Si*(M_Fe*(M_FeO*(-6*M_O + M_SiO2 + 2*M_c) +
    M_MgSiO3*(9*M_O - 4*M_c) - 6*M_O*M_SiO2 + M_c*(M_O + 2*M_SiO2)) + M_FeO*M_SiO2*(9*M_O - 4*M_c) + M_Mg*(M_Fe*(6*M_O -
    M_SiO2 - 2*M_c) + M_FeO*(-6*M_O + M_SiO2 + 2*M_c))) + M_SiO2*(-M_FeO*M_O*M_c + M_Mg*(M_Fe*(-4*M_O + M_c) + M_FeO*(4*M_O
    - M_c)))) + M_Si*(M_Mg*(M_Fe*(M_FeO*(6*M_O - M_SiO2 - 2*M_c) + 6*M_O*M_SiO2 + M_c*(-M_O - 2*M_SiO2)) +
    M_FeO*(M_MgSiO3*(-15*M_O + M_SiO2 + 6*M_c) + M_SiO2*(-9*M_O + 4*M_c))) + M_MgSiO3*(M_Fe*(M_FeO*(-6*M_O + M_SiO2 + 2*M_c)
    - 6*M_O*M_SiO2 + M_c*(M_O + 2*M_SiO2)) + M_FeO*M_SiO2*(9*M_O - 4*M_c)))) +
    M_Mg*dKMgO_KMgO*(M_MgO*(M_Fe*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c) + M_O*M_c*(M_FeO*(-M_FeSiO3
    - M_SiO2) - M_FeSiO3*M_SiO2 + M_MgSiO3*(M_Fe - M_FeO - M_FeSiO3)) + M_Si*(M_Fe*(M_FeSiO3*(-15*M_O + M_SiO2 + 6*M_c) -
    3*M_O*M_SiO2 + M_c*(-M_O + 2*M_SiO2)) + M_FeO*(M_FeSiO3*(9*M_O - 4*M_c) + 3*M_O*M_SiO2 + M_c*(M_O - 2*M_SiO2)) +
    M_FeSiO3*(3*M_O*M_SiO2 + M_c*(M_O - 2*M_SiO2)) + M_MgSiO3*(M_Fe*(-9*M_O + 4*M_c) + M_FeO*(9*M_O - 4*M_c) +
    M_FeSiO3*(9*M_O - 4*M_c)))) + M_MgSiO3*(M_FeO*(M_Fe*M_SiO2*(-4*M_O + M_c) - M_FeSiO3*M_O*M_c) + M_Si*(M_Fe*(M_FeO*(6*M_O
    - M_SiO2 - 2*M_c) + 6*M_O*M_SiO2 + M_c*(-M_O - 2*M_SiO2)) + M_FeO*(M_FeSiO3*(9*M_O - 4*M_c) - 6*M_O*M_SiO2 + M_c*(M_O +
    2*M_SiO2)) + M_FeSiO3*(-6*M_O*M_SiO2 + M_c*(M_O + 2*M_SiO2))))) +
    M_MgSiO3*dKMgSiO3_KMgSiO3*(M_Mg*(M_FeO*(M_Fe*M_SiO2*(-4*M_O + M_c) - M_FeSiO3*M_O*M_c) + M_Si*(M_Fe*(M_FeO*(6*M_O -
    M_SiO2 - 2*M_c) + 6*M_O*M_SiO2 + M_c*(-M_O - 2*M_SiO2)) + M_FeO*(M_FeSiO3*(9*M_O - 4*M_c) - 6*M_O*M_SiO2 + M_c*(M_O +
    2*M_SiO2)) + M_FeSiO3*(-6*M_O*M_SiO2 + M_c*(M_O + 2*M_SiO2)))) + M_MgO*(M_Fe*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O -
    M_SiO2)) + M_SiO2*(M_FeO*(4*M_O - M_c) + M_O*M_c)) + M_Si*(M_Fe*(M_FeO*(-6*M_O + M_SiO2 + 2*M_c) + M_FeSiO3*(-15*M_O +
    M_SiO2 + 6*M_c) + M_SiO2*(-9*M_O + 4*M_c)) + M_Mg*(M_Fe*(6*M_O - M_SiO2 - 2*M_c) + M_FeO*(-6*M_O + M_SiO2 + 2*M_c) +
    M_FeSiO3*(-6*M_O + M_SiO2 + 2*M_c)) + M_SiO2*(M_FeO*(9*M_O - 4*M_c) + M_FeSiO3*(9*M_O - 4*M_c))) +
    M_SiO2*(M_Mg*(M_Fe*(-4*M_O + M_c) + M_FeO*(4*M_O - M_c) + M_FeSiO3*(4*M_O - M_c)) + M_O*M_c*(-M_FeO - M_FeSiO3)))) +
    M_Si*dKSiO2_KSiO2*(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(6*M_O - 2*M_c) + M_SiO2*(2*M_O - M_c)) + M_FeSiO3*(-4*M_O*M_SiO2 +
    M_c*(-M_O + M_SiO2)) - M_O*M_SiO2*M_c) + M_MgSiO3*(M_Fe*(M_FeO*(6*M_O - 2*M_c) - 4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2)) +
    M_FeO*(M_FeSiO3*(-6*M_O + 2*M_c) + 4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) +
    M_O*M_SiO2*M_c*(M_FeO + M_FeSiO3)) + M_MgO*(M_Fe*(M_FeO*(M_FeSiO3*(-6*M_O + 2*M_c) + M_SiO2*(-2*M_O + M_c)) +
    M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_MgSiO3*(M_FeO*(-6*M_O + 2*M_c) + M_FeSiO3*(-6*M_O + 2*M_c)) +
    M_O*M_SiO2*M_c) + M_Mg*(M_Fe*(M_FeSiO3*(6*M_O - 2*M_c) + M_SiO2*(2*M_O - M_c)) + M_FeO*(M_FeSiO3*(-6*M_O + 2*M_c) +
    M_SiO2*(-2*M_O + M_c)) + M_FeSiO3*M_SiO2*(-2*M_O + M_c) + M_MgSiO3*(M_Fe*(6*M_O - 2*M_c) + M_FeO*(-6*M_O + 2*M_c) +
    M_FeSiO3*(-6*M_O + 2*M_c))) + M_O*M_SiO2*M_c*(-M_FeO - M_FeSiO3)) + M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(-6*M_O + 2*M_c) +
    M_SiO2*(-2*M_O + M_c)) + M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c) + M_O*M_SiO2*M_c*(-M_FeO -
    M_FeSiO3))))/(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m +
    M_c*(M_O + M_m))) + M_m*(M_FeSiO3*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2)) - M_O*M_SiO2*M_c)) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_m*(-4*M_O*M_SiO2 + M_c*(-M_O + M_SiO2))) +
    M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) +
    M_FeSiO3*M_m*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2))) + M_O*M_c*(M_FeO*(M_FeSiO3*(-M_SiO2 + M_m) + M_SiO2*M_m) +
    M_FeSiO3*M_SiO2*M_m)) + M_MgO*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) +
    M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_Mg*(M_Fe*(M_FeSiO3*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) + M_SiO2*(-4*M_O*M_m + M_c*(M_O + M_m))) +
    M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) +
    M_FeSiO3*M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m)) + M_MgSiO3*(M_Fe*(M_O*(4*M_SiO2 - 4*M_m) + M_c*(M_O - M_SiO2 + M_m)) +
    M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)))) + M_MgSiO3*(M_Fe*(M_FeO*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 +
    4*M_m) + M_c*(-M_O + M_SiO2 - M_m)) + M_O*M_c*(-M_SiO2 + M_m)) + M_O*M_c*(M_FeO*(M_SiO2 - M_m) + M_FeSiO3*(M_SiO2 -
    M_m))) + M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m) +
    M_Si*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)) + M_Mg*(M_Fe*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O +
    M_SiO2 - M_m)) + M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2
    + M_m)) + M_FeSiO3*(M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + M_m)) + M_MgSiO3*(M_Fe*(-9*M_O + M_SiO2 +
    4*M_c - M_m) + M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m))) +
    M_MgSiO3*(M_Fe*(M_FeO*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(9*M_SiO2 - 9*M_m)
    + M_c*(-4*M_SiO2 + 4*M_m)) + M_FeO*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) + M_FeSiO3*(M_O*(-9*M_SiO2 +
    9*M_m) + M_c*(4*M_SiO2 - 4*M_m))))) + M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(M_O*(-4*M_SiO2 + 4*M_m) + M_c*(-M_O + M_SiO2 -
    M_m)) + M_SiO2*(4*M_O*M_m + M_c*(-M_O - M_m))) + M_m*(M_FeSiO3*(4*M_O*M_SiO2 + M_c*(M_O - M_SiO2)) + M_O*M_SiO2*M_c)) +
    M_O*M_c*(M_FeO*(M_FeSiO3*(M_SiO2 - M_m) - M_SiO2*M_m) - M_FeSiO3*M_SiO2*M_m)) +
    M_Si*(M_Mg*(M_Fe*(M_FeO*(M_FeSiO3*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-M_SiO2 + 4*M_m) - M_SiO2*M_m + M_c*(-M_O +
    M_SiO2 - M_m)) + M_FeSiO3*(M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m + M_c*(-M_O + M_SiO2 - 9*M_m)) + 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(9*M_SiO2 - 9*M_m) + M_c*(-4*M_SiO2 + 4*M_m)) -
    9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeSiO3*(-9*M_O*M_SiO2*M_m + M_c*(M_O*(M_SiO2 - M_m) +
    4*M_SiO2*M_m)) + M_MgSiO3*(M_Fe*(M_FeO*(-9*M_O + M_SiO2 + 4*M_c - M_m) + M_O*(-4*M_SiO2 + 25*M_m) - M_SiO2*M_m +
    M_c*(-M_O + M_SiO2 - 9*M_m)) + M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m +
    M_c*(M_O - M_SiO2 + 9*M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)))) +
    M_MgSiO3*(M_Fe*(M_FeO*(M_FeSiO3*(9*M_O - M_SiO2 - 4*M_c + M_m) + M_O*(M_SiO2 - 4*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 +
    M_m)) + M_FeSiO3*(M_O*(4*M_SiO2 - 25*M_m) + M_SiO2*M_m + M_c*(M_O - M_SiO2 + 9*M_m)) - 9*M_O*M_SiO2*M_m +
    M_c*(M_O*(M_SiO2 - M_m) + 4*M_SiO2*M_m)) + M_FeO*(M_FeSiO3*(M_O*(-9*M_SiO2 + 9*M_m) + M_c*(4*M_SiO2 - 4*M_m)) +
    9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) - 4*M_SiO2*M_m)) + M_FeSiO3*(9*M_O*M_SiO2*M_m + M_c*(M_O*(-M_SiO2 + M_m) -
    4*M_SiO2*M_m)))))

def dM_m_dMi(dMi_m):
    '''alternate way to compute dM_m given the results of the mole changes of all mantle components 
    inputs:
    dMi_m: [dM_MgO, dM_SiO2, dM_FeO, dM_MgSiO3, dM_FeSiO3]
    '''
    return np.sum(dMi_m)
