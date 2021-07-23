#!/usr/bin/env python
'''
Plot sensitivity to axion-photon coupling relative to DFSZ
as a function of photosensor noise equivalent power (NEP)
'''

import matplotlib as mplt
mplt.use('Agg') # So we can use without X forwarding

import numpy as np
import os, json, math, csv, argparse, datetime
import matplotlib.pyplot  as plt
import matplotlib.lines   as mlines
import matplotlib.patches as mpatches
import matplotlib.ticker  as ticker

# So we can produce PDFs
from matplotlib.backends.backend_pdf import PdfPages

doLogY = False # Draw the y-axis on log scale

mplt.rc("text", usetex=True)

#__________________________________________
def main():

  #mkdir('figs') # For the figures

  mk_plot()
  # ----------------------------------------------------------

#__________________________________________
def mk_plot():
  '''
  This plots the contours from the raw (x,y) values of each contour
  '''
  print('Plotting contours for' )
 
  # Figures
  fig, ax = plt.subplots()
  fig.set_size_inches(11, 9)
  text_size = 28

  # Define various colours
  myLightestBlue = '#deebf7'
  myLighterBlue  = '#c6dbef'
  myLightBlue    = '#9ecae1'
  myMediumBlue   = '#6baed6'
  myDarkBlue     = '#4292c6'
  myDarkerBlue   = '#08519c'

  myLightPurple  = '#bcbddc'
  myMediumPurple = '#9e9ac8'
  myPurple       = '#807dba'
  myDarkPurple   = '#54278f'
  myDarkerPurple = '#3f007d'

  myLightestOrange = '#ffffcc'
  myLighterOrange  = '#ffeda0'
  myLightOrange    = '#fed976'
  myMediumOrange   = '#fc4e2a'
  myDarkOrange     = '#993404'

  myLightPink     = '#fcc5c0'
  myDarkPink      = '#ce1256'
  myMediumGreen   = '#41ab5d'
  myDarkGreen     = '#006d2c'

  myLightGray  = '#d9d9d9'
  myMediumGray = '#969696'
  myDarkGray   = '#525252'
  # ------------------------------------------------------- 
  # Sensitivity lines and regions
  # ------------------------------------------------------- 
  Adish = 10.0 # m^2
  Bfield = 10.0

  # Constant rate
  x, y = calc_QCDaxion_coupling(1e-14, 1e-23, Adish, Bfield, snr=5., effic=0.5, time=24., relicDensity = 0.45)
  plt.plot(x, y, lw=2, ls='-', c=myMediumBlue, label=r'$1~\mathrm{day}$') 

  x, y = calc_QCDaxion_coupling(1e-14, 1e-23, Adish, Bfield, snr=5., effic=0.5, time=24000., relicDensity = 0.45)
  plt.plot(x, y, lw=4, ls='-', c=myMediumOrange, label=r'$1000~\mathrm{days}$') 

  x, y = calc_QCDaxion_coupling(1e-18, 1e-21, Adish, Bfield, snr=5., effic=0.5, time=240., relicDensity = 0.45)
  #plt.plot(x, y, lw=4, ls='-', c=myMediumOrange, label=r'$1000~\mathrm{days}$') 

  plt.plot([1e-13, 1e-23], [2.6, 2.6], lw=2, ls='-.', c=myMediumGray) 
  plt.plot([1e-13, 1e-23], [1, 1],     lw=2, ls='--', c=myMediumGray) 
  plt.plot([1e-20, 1e-20], [0.1, 1e2], lw=3, ls='-',  c=myDarkPurple) 

  #plt.plot([1, 1], [1, 1], lw=3, ls=':', c=myDarkGray, label=r'$A_\mathrm{dish}/10$') 
  #plt.plot([1, 1], [1, 1], lw=3, ls='--', c=myDarkGray, label=r'$B_\mathrm{ext}/10$') 
  plt.legend(loc='upper right', prop={'size':text_size*0.95}, frameon=False, handlelength=2.8, borderpad=0.8)

  # ------------------------------------------------------- 
  # Axis properties
  # ------------------------------------------------------- 

  plt.xlim(1e-23, 1e-16)
  plt.ylim(0.1, 1e5)

  # axes labels
  x_txt = r'$\mathrm{Noise~equivalent~power~[W/\sqrt{Hz}}]$'
  y_txt = r'$\left|g_{a\gamma\gamma} / g_{a\gamma\gamma}^\mathrm{DFSZ}\right|~\mathrm{sensitivity}$'
    
  # Axis label properties
  plt.xlabel(x_txt, labelpad=20, size=35) 
  plt.ylabel(y_txt, labelpad=20, size=35)

  ax.set_xscale('log')
  ax.set_yscale('log')
  # Draw the upper THz axis
  #ax2 = ax.twiny()
  #ax2.set_xlim(2.42e-5, 242)
  #ax2.set_xlabel(r'$\mathrm{Frequency}~[\mathrm{THz}]$', labelpad=20, size=35) 
  #ax2.set_xscale('log')

  # ------------------------------------------------------- 
  # Add plot text
  # ------------------------------------------------------- 


  #fig.text(0.89, 0.263, r'$10^{-25}~\mathrm{W}$', color=myMediumGray, size=0.7*text_size)

  #fig.text(0.87, 0.50, r'Photon/second', color=myMediumGray, size=0.5*text_size, rotation=15)
  #fig.text(0.87, 0.375, r'Photon/hour',    color=myMediumGray, size=0.5*text_size, rotation=15)

  #fig.text(0.77, 0.18, r'$A_\mathrm{antenna} = 10~\mathrm{m}^2$', color=myDarkGray, size=text_size)
  fig.text(0.21, 0.77, r'\textbf{BREAD}', color=myDarkGray, size=text_size)
  fig.text(0.21, 0.71, r'$A_\mathrm{dish} = 10~\mathrm{m}^2, B_\mathrm{ext} = 10~\mathrm{T}$', color=myDarkGray, size=text_size)
  fig.text(0.21, 0.65, r'$\mathrm{SNR} = 5, \epsilon_s = 0.5$', color=myDarkGray, size=text_size)

  fig.text(0.71, 0.33, r'$\mathrm{KSVZ}$', color=myDarkGray, size=text_size)
  fig.text(0.71, 0.23, r'$\mathrm{DFSZ}$', color=myDarkGray, size=text_size)

  fig.text(0.28, 0.55, r'$\mathrm{Quantum~capacitance~detector}$', color=myDarkPurple, size=text_size)
  fig.text(0.28, 0.50, r'$\mathrm{(State~of~the~art)}$', color=myDarkPurple, size=text_size*0.8)


  # Adjust axis ticks
  ax.minorticks_on()
  #ax2.minorticks_on()

  ax.tick_params('x', length=12, width=1, which='major', labelsize='28', pad=10)
  ax.tick_params('x', length=6,  width=1, which='minor') 
  #ax2.tick_params('x', length=12, width=1, which='major', labelsize='28', pad=10)
  #ax2.tick_params('x', length=6,  width=1, which='minor') 
  ax.tick_params('y', length=12, width=1, which='major', labelsize='28', pad=10)
  ax.tick_params('y', length=6,  width=1, which='minor') 
   
  plt.tight_layout(pad=0.3)
  plt.subplots_adjust( top=0.85,left=0.17 )
  
  save_name = 'coupling_vs_nep'
  print('Saving as {0}'.format(save_name))
  plt.savefig(save_name + '.pdf', format='pdf', dpi=50)
  #plt.savefig(save_name + '.png', format='png', dpi=400)
  #plt.savefig(save_name + '.eps', format='eps', dpi=150)
  #plt.savefig(save_name, format='png', dpi=150)

#__________________________________________
def csv_to_lists(csv_file):
  '''
  converts csv to dictionary of lists containing columns
  the dictionary keys is the header
  '''
  with open(csv_file) as input_file:
      reader = csv.reader(input_file)
      col_names = next(reader)
      data = {name: [] for name in col_names}
      for line in reader:
        for pos, name in enumerate(col_names):
          data[name].append(line[pos])
  
  return data

#__________________________________________
def calc_QCDaxion_coupling(minNEP, maxNEP, mirrorArea, Bfield, snr=5., effic=0.5, time=1., relicDensity = 0.45):
  '''
  Convert instrument parameters and detected signal power to axion coupling
   - power is units of Watts
   - mirrorArea units is m^2
   - Bfield units is Tesla
   - min/maxNEP   = min and max NEP to plot
   - snr          = required signal to noise ratio 
   - effic        = overall signal power efficiency 
   - time         = integration time in hours
   - relicDensity = dark matter relic density in GeV/cm^3
  Returns lowest [minCoupling, maxCoupling] coupling values probed 
  in units of 10^{-11}/GeV
  '''
  ratio    = ( snr / 5. )
  noiseMin = ( minNEP / 1.0e-21 )
  noiseMax = ( maxNEP / 1.0e-21 )
  area     = ( 10. / mirrorArea  ) 
  dt       = math.sqrt( 1. / time )
  epsilon  = ( 0.5 / effic )
  rho      = ( 0.45 / relicDensity )
  magnet   = ( 10. / Bfield )**2 

  couplingSqMin = 3.6e-24 * ratio * noiseMin * area * dt * epsilon * rho * magnet 
  couplingMin   = math.sqrt( couplingSqMin )
  couplingSqMax = 3.6e-24 * ratio * noiseMax * area * dt * epsilon * rho * magnet 
  couplingMax  = math.sqrt( couplingSqMax )

  gKSVZ = 3.9e-13 #(1/GeV) (mass/meV )
  gDFSZ = 1.5e-13 #(1/GeV) (mass/meV )

  # For the input min and max masses in units of 1 meV
  # Calculate corresponding couplings normalized to DFSZ coupling
  minCoupling = ( couplingMin ) * (1. / gDFSZ)
  maxCoupling = ( couplingMax ) * (1. / gDFSZ)

  return [minNEP, maxNEP], [minCoupling, maxCoupling]

#__________________________________________
def calc_axion_rate(mirrorArea, Bfield, rhoDM, minMass, maxMass, DFSZ=False):
  '''
  Calculate QCD axion photon emission rate with BREAD given input parameters:
  - mirrorArea in m^2
  - Bfield in Tesla
  - rhoDM in GeV/cm^3
  - minMass, maxMass in eV
  - DFSZ if True, else KSVZ axion
  '''

  B = Bfield     / 10.
  A = mirrorArea / 10.
  rho = rhoDM    / 0.45
  
  prefactor = 3.9
  if DFSZ:
    prefactor = 1.5

  minRate = prefactor * (1./minMass) * rho * A * (B**2)
  maxRate = prefactor * (1./maxMass) * rho * A * (B**2)

  return [minMass, maxMass], [minRate, maxRate]

#_________________________________________________________________________
def mkdir(dirPath):
  '''
  make directory for given input path
  '''
  try:
    os.makedirs(dirPath)
    print('Successfully made new directory ' + dirPath)
  except OSError:
    pass

#__________________________________________
if __name__ == "__main__":
    main()