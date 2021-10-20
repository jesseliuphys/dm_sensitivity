#!/usr/bin/env python3
'''
Collect common functions for calculating BREAD axion/dark photon observables/sensitivity
'''
import math

#__________________________________________
def calc_coupling_nep(nep, mirrorArea, Bfield, minMass, maxMass, snr=5., effic=0.5, time=1., relicDensity=0.45, do_axion=True):
  '''
  Convert instrument parameters and detected signal power to axion coupling
   - nep          = noise equivalent power is units of Watts per sqrt(Hz)
   - mirrorArea   = dish area in m^2
   - Bfield       = magnetic field strength in units of Tesla
   - min/maxMass  = min and max mass to plot in meV
   - snr          = required signal to noise ratio 
   - effic        = overall signal power efficiency 
   - time         = integration time in hours
   - relicDensity = dark matter relic density in GeV/cm^3
  Returns [minMass, maxMass] in eV and 
  [minCoupling, maxCoupling] coupling values in units of 1/GeV 
  '''

  snratio  = ( snr / 5. )
  noise    = ( nep / 1.0e-21 )
  rho      = ( 0.45 / relicDensity )
  area     = ( 10.  / mirrorArea  ) 
  epsilon  = ( 0.5  / effic )
  magnet   = ( 10.  / Bfield )**2 
  dt       = math.sqrt( 1. / time )

  # Common experimental parameters
  common_exp_param = snratio * noise * rho * area * epsilon  * dt
  
  # Axion
  minCoupling_axion = math.sqrt( 1.9 * ( minMass ** 2 ) * magnet * common_exp_param ) * 1e-11
  maxCoupling_axion = math.sqrt( 1.9 * ( maxMass ** 2 ) * magnet * common_exp_param ) * 1e-11

  # Dark photon
  minCoupling_dp    = math.sqrt( 7.6 * common_exp_param ) * 1e-14
  maxCoupling_dp    = math.sqrt( 7.6 * common_exp_param ) * 1e-14

  # For the input min and max masses, calculate corresponding couplings
  if do_axion:
    minCoupling = minCoupling_axion
    maxCoupling = maxCoupling_axion
  else:
    minCoupling = minCoupling_dp
    maxCoupling = maxCoupling_dp

  return [minMass*1e-3, maxMass*1e-3], [minCoupling, maxCoupling]

#__________________________________________
def calc_coupling_dcr(dcr, mirrorArea, Bfield, minMass, maxMass, Zsignif=5., effic=0.5, time=1., relicDensity=0.45, do_axion=True):
  '''
  Convert instrument parameters and detected signal power to axion coupling
   - dcr          = dark count rate in Hz
   - mirrorArea   = dish area in m^2
   - Bfield       = magnetic field in Tesla
   - min/maxMass  = min and max mass to plot in meV
   - Zsignif      = required significance 
   - effic        = overall signal power efficiency 
   - time         = integration time in hours
   - relicDensity = dark matter relic density in GeV/cm^3
  Returns [minMass, maxMass] in eV and 
  [minCoupling, maxCoupling] coupling values in units of 1/GeV 
  '''
  zsig     = ( Zsignif / 5. )
  area     = ( 10.  / mirrorArea  ) 
  rho      = ( 0.45 / relicDensity )
  epsilon  = ( 0.5  / effic )
  magnet   = ( 10.  / Bfield )**2 
  noise    = math.sqrt( dcr / 0.01 )
  dt       = math.sqrt( 1. / time )

  # Common experimental parameters
  common_exp_param = zsig * noise * rho * area * epsilon  * dt
  
  # Axion
  minCoupling_axion = math.sqrt( 3.0 * ( minMass ** 3 ) * magnet * common_exp_param ) * 1e-12
  maxCoupling_axion = math.sqrt( 3.0 * ( maxMass ** 3 ) * magnet * common_exp_param ) * 1e-12

  # Dark photon
  minCoupling_dp    = math.sqrt( 11.9 * minMass * common_exp_param ) * 1e-15
  maxCoupling_dp    = math.sqrt( 11.9 * maxMass * common_exp_param ) * 1e-15

  # For the input min and max masses, calculate corresponding couplings
  if do_axion:
    minCoupling = minCoupling_axion
    maxCoupling = maxCoupling_axion
  else:
    minCoupling = minCoupling_dp
    maxCoupling = maxCoupling_dp

  return [minMass*1e-3, maxMass*1e-3], [minCoupling, maxCoupling]

#__________________________________________
if __name__ == "__main__":
    main()