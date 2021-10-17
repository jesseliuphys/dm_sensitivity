#!/usr/bin/env python3
'''
Plot signal photon rate
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
  This plots the photon rate vs QCD axion mass
  '''
  print('Plotting contours for' )
 
  # Figures
  fig, ax = plt.subplots()
  fig.set_size_inches(11, 9)
  text_size = 32

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

  # Constant rate
  x, y = calc_axion_rate(10, 10, 0.45, 1e-6, 10, DFSZ=False)
  plt.plot(x, y, lw=4, ls='-', c=myMediumBlue, label=r'$\mathrm{KSVZ}$') 

  x, y = calc_axion_rate(10, 10, 0.45, 1e-6, 10, DFSZ=True)
  plt.plot(x, y, lw=2, ls='-', c=myMediumOrange, label=r'$\mathrm{DFSZ}$') 

  # A = 1 m^2
  x, y = calc_axion_rate(1, 10, 0.45, 1e-6, 10, DFSZ=False)
  plt.plot(x, y, lw=4, ls=':', c=myMediumBlue) 
  x, y = calc_axion_rate(1, 10, 0.45, 1e-6, 10, DFSZ=True)
  plt.plot(x, y, lw=2, ls=':', c=myMediumOrange) 

  # B = 1 T
  x, y = calc_axion_rate(10, 1, 0.45, 1e-6, 10, DFSZ=False)
  plt.plot(x, y, lw=4, ls='--', c=myMediumBlue) 
  x, y = calc_axion_rate(10, 1, 0.45, 1e-6, 10, DFSZ=True)
  plt.plot(x, y, lw=2, ls='--', c=myMediumOrange) 

  plt.plot([1, 1], [1, 1], lw=3, ls=':', c=myDarkGray, label=r'$A_\mathrm{dish}/10$') 
  plt.plot([1, 1], [1, 1], lw=3, ls='--', c=myDarkGray, label=r'$B_\mathrm{ext}/10$') 
  plt.legend(loc='upper right', prop={'size':text_size*0.95}, frameon=False, handlelength=2.8, borderpad=0.8)

  # ------------------------------------------------------- 
  # Axis properties
  # ------------------------------------------------------- 

  plt.xlim(1e-6, 10)
  plt.ylim(0.5, 5e6)

  # axes labels
  x_txt = r'$\mathrm{Axion~mass~[meV}]$'
  y_txt = r'$\mathrm{Emitted~photons~[day}^{-1}]$'
    
  # Axis label properties
  plt.xlabel(x_txt, labelpad=20, size=40) 
  plt.ylabel(y_txt, labelpad=20, size=40)

  ax.set_xscale('log')
  ax.set_yscale('log')

  # Draw the upper THz axis
  ax2 = ax.twiny()
  ax2.set_xlim(2.42e-4, 2417.98)
  ax2.set_xlabel(r'$\mathrm{Frequency}~[\mathrm{THz}]$', labelpad=20, size=35)
  ax2.set_xscale('log')

  # Rescale axis from eV to meV
  scale_x=1e3
  ticks_x = mplt.ticker.FuncFormatter(lambda x, pos: r'${0:g}$'.format(x*scale_x))
  ax.xaxis.set_major_formatter(ticks_x)
  ticks_x = mplt.ticker.FuncFormatter(lambda x, pos: r'${0:g}$'.format(x))
  ax2.xaxis.set_major_formatter(ticks_x)

  # ------------------------------------------------------- 
  # Add plot text
  # ------------------------------------------------------- 
  fig.text(0.21, 0.34, r'\textbf{BREAD}', color=myDarkGray, size=text_size)
  fig.text(0.21, 0.28, r'$A_\mathrm{dish} = 10~\mathrm{m}^2$', color=myDarkGray, size=text_size)
  fig.text(0.21, 0.22, r'$B_\mathrm{ext} = 10~\mathrm{T}$', color=myDarkGray, size=text_size)

  # Adjust axis ticks
  ax.minorticks_on()
  #ax2.minorticks_on()

  ax.tick_params('x', length=12, width=1, which='major', labelsize='33', pad=10, direction="in")
  ax.tick_params('x', length=6,  width=1, which='minor', direction="in") 
  ax2.tick_params('x', length=12, width=1, which='major', labelsize='33', pad=10, direction="in")
  ax2.tick_params('x', length=6,  width=1, which='minor', direction="in") 
  ax.tick_params('y', length=12, width=1, which='major', labelsize='33', pad=10, direction="in", right="on")
  ax.tick_params('y', length=6,  width=1, which='minor', direction="in", right="on") 
   
  # Force axis ticks to appear in log scale
  locmaj = mplt.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=100)
  locmin = mplt.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=100)
  ax.xaxis.set_major_locator(locmaj)
  ax.xaxis.set_minor_locator(locmin)
  locmaj = mplt.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=100)
  locmin = mplt.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=100)
  ax2.xaxis.set_major_locator(locmaj)
  ax2.xaxis.set_minor_locator(locmin)
  locmaj = mplt.ticker.LogLocator(base=10.0, subs=(1.0, ), numticks=100)
  locmin = mplt.ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*.1,numticks=100)
  ax.yaxis.set_major_locator(locmaj)
  ax.yaxis.set_minor_locator(locmin)
  ax2.xaxis.set_minor_formatter(mplt.ticker.NullFormatter())
  ax.xaxis.set_minor_formatter(mplt.ticker.NullFormatter())
  ax.yaxis.set_minor_formatter(mplt.ticker.NullFormatter())

  plt.tight_layout(pad=0.3)
  plt.subplots_adjust( top=0.83,left=0.17, bottom=0.18 )
  
  save_name = 'photon_rate_qcdaxion'
  print('Saving as {0}'.format(save_name))
  plt.savefig(save_name + '.pdf', format='pdf', dpi=50)
  #plt.savefig(save_name + '.png', format='png', dpi=400)
  #plt.savefig(save_name + '.eps', format='eps', dpi=150)
  #plt.savefig(save_name, format='png', dpi=150)

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