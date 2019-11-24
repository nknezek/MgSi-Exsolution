# Mg Si Parameters:

### Exsolution parameters:

alpha_FeO = 0.6 % (Jackie's estimate from FeO vs Fe density)

alpha_SiO2 = 1.17 % (Hirose 2017)

alpha_MgO = 0.89 % (Hard-shell density estimate with Van Der Waals radius)

dE_FeO = 1.01 MJ/kg (Van't Hoff equation given b values from Hirose)

dE_SiO2 = 4.307 MJ/kg (Van't Hoff equation given b values from Fischer)

dE_MgO = 9.005 MJ/kg  (Van't Hoff equation given b values from Badro)

### Mantle parameters

viscosity_now = 10^21 Pa-s

timescale of mantle overturn = 600 Myrs at present

layer thickness = 100 m 

### Cottaar mantle composition: 

Mole frac Pv = .674, Fp = 0.311, Stishovite = 0.015

Mg# Fp = 0.8

Mg# Pv = 0.93

### Parameters to vary

Temperature at CMB initial = linspace(4800K, 6500K, 100K) = 18 values

wt% Mg, Si = linspace(1e-5 %, 5 %, 0.5%) = 11 x 2 values

wt% O = linspace(1e-5 %, 15 %, 0.5%) = 31 values

 = 18 x 11 x 11 x 31 = 67,518 runs to compute