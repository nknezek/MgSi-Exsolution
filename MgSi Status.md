# MgSi Status

## 2/21/18

* Tushar - run lot's o shit
  * T_cmb = [ 5000,5250, 5500, 5750, 6000]
  * dT = [2400, 2600, 2800]
  * nu = 10**[19, 19.6, 20.2, 20.8]
  * X_Mg = [1e-5, 0.01, 0.02, 0.04]
  * X_Si = [1e-5, 0.01, 0.02, 0.04]
  * X_O = [0.10, 0.16, 0.22] 
    * today - low wt%? m
  * overturn = [800, 400, 200]
  * layer thickness = [30, 300,1000]

* Nick 
  * email collabs
    - parameter space of runs
    - background mantle composition
  * optimize interpolation plotting memory use
  * look at AWS to run many code? 
    * Google cloud 
    * Open Science Grid? osgconnect.net
  * pull Olson dipole scaling to get dipole moment 
    * Olson and Christensen 2013 Dipole Moment scaling
  * powering dynamo from below vs both boundaries
    * driven from below vs driven both boundaries

* Plots to produce
  * cartoon diagram of system
    - mantle dynamics / box diagram
  * thermal history / Ric
    * include dynamo power in subplot?
  * parameter tradeoff grid
    * 5x5
    * T_cmb, X_Mg, X_Si, X_O, layer thick
    * histogram KDE along axis
  * ensemble plot of dynamo power
    * all runs within 5% of ic size
    * lower axis is histogram of ic nucleation date?
  * inner core nucleation signal histogram ? 
    * d E_phi / (E_phi0 * dt)
    * dt ~ 10Myr
    * X_Si, X_Mg, means of change in power?

* supplement Plots
  * 7x7 parameter tradeoff
  * X_Mg vs X_Si grid for 
    * T_cmb vs nu?
    * overturn vs layter thick?
    * ​

  ### 3/8/18 Status

  Todo: Implement change in Kd over time to push things into Bridgemanite

  maybe: link overturn time to mantle temperature

  ### 3/8/18 Status

  * Notable problem: Core is significantly sub-adiabatic throughout most of it's history, but our core assumes adiabatic profile
  * lots of Mg and Si probably indicates significantly sub-adiabatic heatflow at the present day — strong thermal boundary layer
  * ​