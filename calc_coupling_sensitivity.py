#!/usr/bin/env python3
'''
Printout sensitivity to various couplings
Given various experimental parameters of BREAD
'''

import math

#__________________________________________
def main():
  '''
  Display sensitivity to couplings given experimental parameters
  '''

  # ------------------------------------------------------- 
  # Sensitivity lines and regions
  # ------------------------------------------------------- 
  Adish  = 10.0 # m^2
  Bfield = 10.0 # Tesla
  hours  = 24000. # hours
  nep    = 1e-22 # W/sqrt(Hz)

  # Constant rate
  couplingKSVZ, couplingDFSZ = calc_QCDaxion_coupling(nep, Adish, Bfield, snr=5., effic=0.5, time=hours, relicDensity = 0.45)
  kineticMix = calc_darkPhoton_coupling(nep, Adish, snr=5., effic=0.5, time=hours, relicDensity = 0.45)

  printout = 'axion/KSVZ: {0}\naxion/DFSZ: {1}\ndark-photon/1e-14: {2}'.format(couplingKSVZ, couplingDFSZ, kineticMix)
  print(printout)

#__________________________________________
def calc_darkPhoton_coupling(nep, mirrorArea, snr=5., effic=0.5, time=1., relicDensity = 0.45):
  '''
  Convert instrument parameters and detected signal power to dark photon coupling
   - nep          = overall noise equivalent power of photonsensor in units of Watts/sqrt(Hz)
   - mirrorArea   = area of dish antenna units is m^2
   - min/maxMass  = min and max mass to plot
   - snr          = required signal to noise ratio 
   - effic        = overall signal power efficiency 
   - time         = integration time in hours
   - relicDensity = dark matter relic density in GeV/cm^3
  '''

  snratio  = ( snr / 5. )
  noise    = ( nep / 1.0e-21 )
  rho      = ( 0.45 / relicDensity )
  area     = ( 10.  / mirrorArea  ) 
  epsilon  = ( 0.5  / effic )
  dt       = math.sqrt( 1. / time )

  # Common experimental parameters
  common_exp_param = snratio * noise * rho * area * epsilon  * dt

  kineticMixing    = math.sqrt( 7.6 * common_exp_param ) * 1e-14

  return kineticMixing / 1e-14

#__________________________________________
def calc_QCDaxion_coupling(nep, mirrorArea, Bfield, snr=5., effic=0.5, time=1., relicDensity = 0.45):
  '''
  Convert instrument parameters and detected signal power to axion coupling
   - power is units of Watts
   - mirrorArea units is m^2
   - Bfield units is Tesla
   - nep          = NEP considered
   - snr          = required signal to noise ratio 
   - effic        = overall signal power efficiency 
   - time         = integration time in hours
   - relicDensity = dark matter relic density in GeV/cm^3
  Returns lowest [minCoupling, maxCoupling] coupling values probed 
  in units of 10^{-11}/GeV
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
  coupling = math.sqrt( 1.9 * magnet * common_exp_param ) * 1e-11

  gKSVZ = 3.9e-13 #(1/GeV) (mass/meV )
  gDFSZ = 1.5e-13 #(1/GeV) (mass/meV )

  # For the input min and max masses in units of 1 meV
  # Calculate corresponding couplings normalized to DFSZ coupling
  couplingKSVZ = ( coupling ) * (1. / gKSVZ)
  couplingDFSZ = ( coupling ) * (1. / gDFSZ)

  return couplingKSVZ, couplingDFSZ

#__________________________________________
if __name__ == "__main__":
    main()