#!/usr/bin/env python3
'''
Dark photon DM sensitivity plotter
'''

import matplotlib as mplt
mplt.use('Agg') # So we can use without X forwarding

import numpy as np
import pandas as pd
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


  fig, ax = plt.subplots()
  fig.set_size_inches(11, 9)

  # Define various colours
  myLightestBlue   = '#deebf7'
  myLighterBlue    = '#c6dbef'
  myLightBlue      = '#9ecae1'
  myMediumBlue     = '#6baed6'
  myDarkBlue       = '#4292c6'
  myDarkerBlue     = '#08519c'

  myLightPurple    = '#bcbddc'
  myMediumPurple   = '#9e9ac8'
  myPurple         = '#807dba'
  myDarkPurple     = '#54278f'
  myDarkerPurple   = '#3f007d'

  myLightestOrange = '#ffffcc'
  myLighterOrange  = '#ffeda0'
  myLightOrange    = '#fed976'
  myMediumOrange   = '#fc4e2a'
  myDarkOrange     = '#993404'

  myLightPink      = '#fcc5c0'
  myDarkPink       = '#ce1256'
  myMediumGreen    = '#41ab5d'
  myDarkGreen      = '#006d2c'

  myLightGray      = '#d9d9d9'
  myMediumGray     = '#969696'
  myDarkGray       = '#525252'

  # Existing constraints from https://arxiv.org/abs/2105.04565 https://github.com/cajohare/AxionLimits/ 
  haloscope= pd.read_csv('limits/common/DP_Combined_Haloscope.txt',sep=' ', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(haloscope['x'], haloscope['y'], 1e-1, edgecolor='none', facecolor=myLightBlue)

  dmsearch = pd.read_csv('limits/common/DP_Combined_DarkMatterSearches.txt', sep=' ', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(dmsearch['x'], dmsearch['y'], 1e-1, edgecolor='none', facecolor=myLightBlue)  

  dmcosmo  = pd.read_csv('limits/common/DP_Combined.txt', sep=' ', lineterminator='\n', names=['x', 'y'])
  plt.plot(dmcosmo['x'],   dmcosmo['y'],  color=myLighterBlue, lw=1)
  plt.fill_between(dmcosmo['x'], dmcosmo['y'], 1e-1, edgecolor='none', facecolor=myLightestBlue)

  stellar  = pd.read_csv('limits/common/DP_Combined_Stellar.txt', sep=' ',lineterminator='\n', names=['x', 'y'])
  plt.plot(stellar['x'],   stellar['y'],  color=myLightBlue, lw=1, zorder=-1)
  plt.fill_between(stellar['x'], stellar['y'], 1e-1, edgecolor='none', facecolor=myLighterBlue)
  
  lab      = pd.read_csv('limits/common/DP_Combined_Laboratory.txt', sep=' ', lineterminator='\n', names=['x', 'y'])
  plt.plot(lab['x'], lab['y'], color=myLightBlue, lw=1)
  plt.fill_between(lab['x'], lab['y'], 1e-1, edgecolor='none', facecolor=myLightestBlue)

  print(stellar)

  # Default values
  Adish = 10. # dish area in m^2
  snr   = 5.  # signal to noise
  effic = 0.5 # signal detection efficiency
  time  = 240.  # integration time in hours

  #-----------------------------
  # 10x time of Upgrade 1
  #-----------------------------
  # Bolometer [1.65, 83] meV
  x, y = calc_darkPhoton_coupling(4e-17, Adish, 0.24, 248, snr, effic, time*100)
  plt.plot(x, y, alpha=0.7,lw=3, ls='--', c=myDarkGray, zorder=4) 
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myMediumGray, alpha=0.1)

  # KID [0.2, 5] meV
  x, y = calc_darkPhoton_coupling(3e-21, Adish, 0.2, 5, snr, effic, time*100)
  plt.plot(x, y, myLightOrange, alpha=0.7,lw=3, ls='--', c=myMediumOrange, zorder=4) 
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myLightOrange, alpha=0.2)

  # TES [0.2, 1.2] meV
  x, y = calc_darkPhoton_coupling(2e-21, Adish, 0.19, 1.2, snr, effic, time*100)
  plt.plot(x, y, alpha=0.8,lw=3, ls='--', c=myDarkPurple, zorder=4) 
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myDarkPurple, alpha=0.15)

  # SNSPD [207, 830] meV
  x, y = calc_darkPhoton_coupling(1e-21, Adish, 207, 830, snr, effic, time*100)
  plt.plot(x, y, alpha=0.3,lw=3, ls='--', c=myDarkPink, zorder=4) 
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myDarkPink, alpha=0.1)

  # QCD [6.2] meV
  x, y = calc_darkPhoton_coupling(1e-22, Adish, 6.18, 6.22, snr, effic, time*100)
  plt.fill_between(x, y, 1e-1, color=myDarkGray, alpha=0.3, zorder=4)

  #-----------------------------
  # Upgrade 1
  #-----------------------------
  # Bolometer [1.65, 83] meV
  x, y = calc_darkPhoton_coupling(4e-17, Adish, 0.24, 248, snr, effic, time)
  plt.plot(x, y, alpha=0.7,lw=3, ls='-.', c=myDarkGray, zorder=4) 
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myMediumGray, alpha=0.1, zorder=4)

  # KID [0.2, 5] meV
  x, y = calc_darkPhoton_coupling(3e-21, Adish, 0.2, 5, snr, effic, time)
  plt.plot(x, y, myLightOrange, alpha=0.7,lw=3, ls='-.', c=myMediumOrange, zorder=4) 
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myLightOrange, alpha=0.2, zorder=4)
  
  # TES [0.2, 1.2] meV
  x, y = calc_darkPhoton_coupling(2e-21, Adish, 0.19, 1.2, snr, effic, time)
  plt.plot(x, y, alpha=0.8, lw=3, ls='-.',c=myDarkPurple, zorder=4) 

  # SNSPD [207, 830] meV
  x, y = calc_darkPhoton_coupling(1e-21, Adish, 207, 830, snr, effic, time)
  plt.plot(x, y, alpha=0.3,lw=3, ls='-.', c=myDarkPink, zorder=4) 
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myDarkPink, alpha=0.1, zorder=4)

  #-----------------------------
  # Baseline
  #-----------------------------
  # Far-IR 1.6 K IR Labs bolometer
  x, y = calc_darkPhoton_coupling(4e-15, Adish, 0.24, 248, snr, effic, time)
  plt.plot(x, y, alpha=0.8,lw=3, c=myDarkGray, zorder=4) 
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myMediumGray, alpha=0.3, zorder=4)

  # KID [0.2, 5] meV
  x, y = calc_darkPhoton_coupling(3e-19, Adish, 0.2, 5, snr, effic, time)
  plt.plot(x, y, alpha=0.8,lw=3, c=myMediumOrange, zorder=4) 
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myLightOrange, alpha=0.6, zorder=4)

  # TES [0.2, 1.2] meV
  x, y = calc_darkPhoton_coupling(2e-19, Adish, 0.19, 1.2, snr, effic, time)
  plt.plot(x, y, alpha=0.8,lw=3, c=myDarkPurple, zorder=4) 
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myDarkPurple, alpha=0.15, zorder=4)

  # SNSPD [207, 830] meV
  x, y = calc_darkPhoton_coupling(1e-18, Adish, 207, 830, snr, effic, time)
  plt.plot(x, y, alpha=0.8,lw=3, c=myDarkPink, zorder=4) 
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myLightPink, alpha=0.6, zorder=4)

  # QCD [6.2] meV
  x, y = calc_darkPhoton_coupling(1e-20, Adish, 6.18, 6.22, snr, effic, time)
  #plt.plot(x, y, myMediumGray, alpha=0.8,lw=3, c=myDarkGray, zorder=4, label=r'$\mathrm{QCD}~10^{-20}~\mathrm{W}~\mathrm{Hz}^{-1/2}$') 
  plt.fill_between(x, y, 1e-1, color=myDarkGray, alpha=0.9, zorder=4)

  # Gentec [0.4, 120] meV
  x, y = calc_darkPhoton_coupling(1e-8, Adish, 0.4, 120, snr, effic, time)
  plt.plot(x, y, lw=3, ls='-', c=myDarkGreen) 
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myMediumGreen, alpha=0.3, zorder=4)

  #-----------------------------
  # Pilot dish size and 1 day 
  #-----------------------------

  # Far-IR 1.6 K IR Labs bolometer
  x, y = calc_darkPhoton_coupling(4e-15, Adish/10.0, 0.24, 248, snr, effic, time/10.0)
  plt.plot(x, y, alpha=0.8,lw=1, c=myDarkGray, zorder=4) 

  # SNSPD [207, 830] meV
  x, y = calc_darkPhoton_coupling(1e-18, Adish/10.0, 207, 830, snr, effic, time/10.0)
  plt.plot(x, y, alpha=0.8,lw=1, c=myDarkPink, zorder=4) 

  # Existing constraint labels
  text_size = 23
  fig.text(0.175, 0.45, r'Haloscope',      color=myDarkerBlue, size=text_size)
  fig.text(0.22, 0.67,  r'Cosmology',       color=myDarkBlue, size=text_size)
  fig.text(0.86, 0.63,  r'Stellar',         color=myDarkBlue, size=text_size)
  fig.text(0.23, 0.78,  r"$\gamma \to A'$", color=myDarkBlue, size=text_size)
  
  # Sensors labels
  fig.text(0.61, 0.78,r'Pyroelectric', color=myDarkGreen,   size=text_size)
  fig.text(0.61, 0.755,r'(293 K commercial)', color=myDarkGreen,   size=text_size*0.5)
  fig.text(0.63,0.62, r'IR Labs',      color=myDarkGray,    size=text_size)
  fig.text(0.63,0.60, r'(Commercial)', color=myDarkGray,    size=text_size*0.5)

  fig.text(0.78, 0.50, r'SNSPD',       color=myDarkPink,    size=text_size)
  fig.text(0.54, 0.43,r'KID',         color=myMediumOrange,size=text_size)
  fig.text(0.43, 0.39, r'TES',         color=myDarkPurple,  size=text_size)
  fig.text(0.61, 0.35, r'QCDet',       color=myDarkGray,    size=text_size, rotation=90)

  # Axis log scale
  ax.set_xscale('log')
  ax.set_yscale('log')
  # Axis limits
  ax.set_xlim(1e-6, 10)
  ax.set_ylim(1e-18, 1e-6)
  # Axis labels
  x_txt = r"$m_{A'}~[\mathrm{eV}]$"
  y_txt = r'$\kappa$'
  # Axis label properties
  plt.xlabel(x_txt, labelpad=20, size=35) 
  plt.ylabel(y_txt, labelpad=20, size=35)

  # Draw the upper THz axis
  ax2 = ax.twiny()
  ax2.set_xlim(241.79e-6, 2417.9)
  ax2.set_xscale('log')
  ax2.set_xlabel(r'$\nu~[\mathrm{THz}]$', labelpad=18, size=35)

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

  # Axis ticks
  ax.tick_params('x', length=12, width=1, which='major', labelsize='28', pad=10, direction="in")
  ax.tick_params('x', length=6,  width=1, which='minor', direction="in") 
  ax2.tick_params('x', length=12, width=1, which='major', labelsize='28', pad=10, direction="in")
  ax2.tick_params('x', length=6,  width=1, which='minor', direction="in") 
  ax.tick_params('y', length=12, width=1, which='major', labelsize='28', pad=10, direction="in", right="on")
  ax.tick_params('y', length=6,  width=1, which='minor', direction="in", right="on") 

  # Legend
  plt.plot([1, 1], [1, 1], lw=1, ls='-',  c=myDarkGray, label=r'NEP$_\mathsf{today}$, 1 day, $A_\mathrm{dish}/10$') 
  plt.plot([1, 1], [1, 1], lw=3, ls='-',  c=myDarkGray, label=r'NEP$_\mathsf{today}$, 10 days') 
  plt.plot([1, 1], [1, 1], lw=3, ls='-.', c=myDarkGray, label=r'NEP$_\mathsf{today}/100$, 10 days') 
  plt.plot([1, 1], [1, 1], lw=3, ls='--', c=myDarkGray, label=r'NEP$_\mathsf{today}/100$, 1000 days') 
  plt.legend(loc='lower right', prop={'size':15}, frameon=False, handlelength=2.8, borderpad=0.6)

  fig.text(0.36, 0.22, r'\textbf{BREAD} $A_\mathrm{dish} = ' + '{0}'.format(int(Adish)) + '~\mathrm{m}^2$', color=myDarkGray, size=text_size)
  
  if time == 8760:
    integrationT = ', $\Delta t_\mathrm{int} = 1~\mathrm{yr}$'
  elif time == 87600:
    integrationT = ', $\Delta t_\mathrm{int} = 10~\mathrm{yrs}$'
  else:
    integrationT = ', $\Delta t_\mathrm{int}' + ' = {0:.0f}'.format(time) + '~\mathrm{hrs}$'
  #fig.text(0.40, 0.175, r'$\mathrm{SNR}' + ' = {0:.0f}$'.format(snr) + ', $\epsilon_\mathrm{sig}' + ' = {0}$'.format(effic) + integrationT, color=myDarkGray, size=text_size*0.68)
  fig.text(0.36, 0.18, r'$\mathrm{SNR}' + ' = {0:.0f}$'.format(snr) + ', $\epsilon_\mathrm{sig}' + ' = {0}$'.format(effic), color=myDarkGray, size=text_size*0.9)

  # Plot margins
  plt.tight_layout(pad=0.3)
  plt.subplots_adjust( top=0.85,left=0.15 )

  # Save plot to pdf
  plt.savefig('fig_dark_photon_photosensors.pdf', format='pdf')

#__________________________________________
def calc_darkPhoton_coupling(nep, mirrorArea, minMass, maxMass, snr=5., effic=0.5, time=1., relicDensity = 0.45):
  '''
  Convert instrument parameters and detected signal power to dark photon coupling
   - nep          = overall noise equivalent power of photonsensor in units of Watts/sqrt(Hz)
   - mirrorArea   = area of dish antenna units is m^2
   - min/maxMass  = min and max mass to plot (input in meV, output in eV)
   - snr          = required signal to noise ratio 
   - effic        = overall signal power efficiency 
   - time         = integration time in hours
   - relicDensity = dark matter relic density in GeV/cm^3
  '''

  ratio    = ( snr / 5. )
  noise    = ( nep / 1.0e-21 )
  area     = ( 10.0 / mirrorArea  ) 
  dt       = math.sqrt( 1. / time )
  epsilon  = ( 0.5 / effic )
  rho      = ( 0.3 / relicDensity )

  kineticMixingSq = 3.7 * ratio * noise * area * dt * epsilon * rho
  kineticMixing   = math.sqrt( kineticMixingSq ) * 1e-14
 
  return [minMass*1e-3, maxMass*1e-3], [kineticMixing, kineticMixing]

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