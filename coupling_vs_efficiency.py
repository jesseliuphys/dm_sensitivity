#!/usr/bin/env python3
'''
Plot sensitivity to axion-photon coupling relative to DFSZ
as a function of overall emission-to-detection efficiency
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
mplt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

#__________________________________________
def main():
  '''
  This plots the contours from the raw (x,y) values of each contour
  '''
  print('Plotting contours for' )
 
  # Figures
  fig, ax = plt.subplots()
  fig.set_size_inches(11, 9)
  text_size = 38

  # Define various colours
  myMediumBlue   = '#6baed6'
  myDarkPurple   = '#54278f'
  myMediumOrange = '#fc4e2a'
  myMediumGray   = '#969696'
  myDarkGray     = '#525252'
  # ------------------------------------------------------- 
  # Sensitivity lines and regions
  # ------------------------------------------------------- 
  Adish = 10.0 # m^2
  Bfield = 10.0

  # Constant rate
  x, y = calc_QCDaxion_coupling(1e-21, Adish, Bfield, snr=5., minEffic=1e-3, maxEffic=1., time=24000., relicDensity = 0.45)
  plt.plot(x, y, lw=2, ls='-', c=myMediumBlue, label=r'$10^{-21}~\mathrm{W~Hz}^{-1/2}$') 

  x, y = calc_QCDaxion_coupling(1e-23, Adish, Bfield, snr=5., minEffic=1e-3, maxEffic=1., time=24000., relicDensity = 0.45)
  plt.plot(x, y, lw=4, ls='-', c=myMediumOrange, label=r'$10^{-23}~\mathrm{W~Hz}^{-1/2}$') 


  plt.plot([1e-3, 1], [2.6, 2.6], lw=2, ls='-.', c=myMediumGray) 
  plt.plot([1e-3, 1], [1, 1],     lw=2, ls='--', c=myMediumGray) 
  #plt.plot([1e-3, 1], [0.1, 1e2], lw=3, ls='-',  c=myDarkPurple) 

  plt.legend(loc='upper right', prop={'size':text_size*0.95}, frameon=False, handlelength=1.5, borderpad=0.6, handletextpad=0.4)

  # ------------------------------------------------------- 
  # Axis properties
  # ------------------------------------------------------- 

  plt.xlim(1e-3, 1)
  plt.ylim(0.1, 9e3)

  # axes labels
  x_txt = r'$\mathrm{Signal~efficiency}$'
  y_txt = r'$\left|g_{a\gamma\gamma} / g_{a\gamma\gamma}^\mathrm{DFSZ}\right|~\mathrm{sensitivity}$'
    
  # Axis label properties
  plt.xlabel(x_txt, labelpad=10, size=40) 
  plt.ylabel(y_txt, labelpad=5, size=40)

  ax.set_xscale('log')
  ax.set_yscale('log')
  ax.get_xaxis().set_major_formatter(mplt.ticker.ScalarFormatter())
  ax.get_xaxis().set_major_formatter(mplt.ticker.FormatStrFormatter(r'$%g$'))
  ax.get_yaxis().set_major_formatter(mplt.ticker.ScalarFormatter())
  ax.get_yaxis().set_major_formatter(mplt.ticker.FormatStrFormatter(r'$%g$'))

  # ------------------------------------------------------- 
  # Add plot text
  # ------------------------------------------------------- 


  #fig.text(0.77, 0.18, r'$A_\mathrm{antenna} = 10~\mathrm{m}^2$', color=myDarkGray, size=text_size)
  fig.text(0.23, 0.88, r'\textbf{BREAD}', color=myDarkGray, size=text_size)
  fig.text(0.23, 0.81, r'$A_\mathrm{dish} = 10~\mathrm{m}^2$', color=myDarkGray, size=text_size)
  fig.text(0.23, 0.74, r'$B_\mathrm{ext} = 10~\mathrm{T}$',    color=myDarkGray, size=text_size)
  fig.text(0.33, 0.66, r'$\mathrm{SNR} = 5, \Delta t = 1000~\mathrm{days}$', color=myDarkGray, size=text_size)

  fig.text(0.24, 0.42, r'$\mathrm{KSVZ}$', color=myDarkGray, size=text_size)
  fig.text(0.24, 0.27, r'$\mathrm{DFSZ}$', color=myDarkGray, size=text_size)

  # Adjust axis ticks
  ax.minorticks_on()
  ax.tick_params('x', length=12, width=1, which='major', labelsize='35', pad=20, direction="in", top="on")
  ax.tick_params('x', length=6,  width=1, which='minor', direction="in", top="on") 
  ax.tick_params('y', length=12, width=1, which='major', labelsize='35', pad=20, direction="in", right="on")
  ax.tick_params('y', length=6,  width=1, which='minor', direction="in", right="on") 
   
  plt.tight_layout(pad=0.3)
  plt.subplots_adjust( top=0.97, left=0.19, bottom=0.17 )
  
  save_name = 'coupling_vs_efficiency'
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
def calc_QCDaxion_coupling(NEP, mirrorArea, Bfield, snr=5., minEffic=1e-4, maxEffic=1.0, time=1., relicDensity = 0.45):
  '''
  Convert instrument parameters and detected signal power to axion coupling
   - power is units of Watts
   - mirrorArea units is m^2
   - Bfield units is Tesla
   - NEP          = NEP to plot
   - snr          = required signal to noise ratio 
   - min/maxEffic = overall signal power efficiency 
   - time         = integration time in hours
   - relicDensity = dark matter relic density in GeV/cm^3
  Returns lowest [minCoupling, maxCoupling] coupling values probed 
  in units of 10^{-11}/GeV
  '''
  ratio    = ( snr / 5. )
  noise    = ( NEP / 1.0e-21 )
  area     = ( 10. / mirrorArea  ) 
  dt       = math.sqrt( 1. / time )
  minEpsilon  = ( 0.5 / minEffic )
  maxEpsilon  = ( 0.5 / maxEffic )
  rho      = ( 0.45 / relicDensity )
  magnet   = ( 10. / Bfield )**2 

  couplingSqMin = 1.9 * ratio * noise * area * dt * minEpsilon * rho * magnet 
  couplingMin   = math.sqrt( couplingSqMin ) * 1e-11
  couplingSqMax = 1.9 * ratio * noise * area * dt * maxEpsilon * rho * magnet 
  couplingMax   = math.sqrt( couplingSqMax ) * 1e-11

  gDFSZ = 1.5e-13 #(1/GeV) (mass/meV )

  # For the input min and max masses in units of 1 meV
  # Calculate corresponding couplings normalized to DFSZ coupling
  minCoupling = ( couplingMin ) * (1. / gDFSZ)
  maxCoupling = ( couplingMax ) * (1. / gDFSZ)

  return [minEffic, maxEffic], [minCoupling, maxCoupling]

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