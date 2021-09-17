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

  myMediumRed      = '#ef3b2c'
  myDarkRed        = '#a50f15'

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

  shuket = pd.read_csv('limits/common/DP_SHUKET.txt', sep=' ', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(shuket['x'], shuket['y'], 1e-1, edgecolor='none', facecolor=myDarkerBlue)  

  tokyo1 = pd.read_csv('limits/common/DP_Tokyo-Tomita.txt', sep=' ', lineterminator='\n', names=['x', 'y'])
  plt.plot(tokyo1['x'],   tokyo1['y'],  color=myDarkerBlue, lw=1, zorder=-2)
  plt.fill_between(tokyo1['x'], tokyo1['y'], 1e-1, edgecolor='none', facecolor=myDarkerBlue)  

  dmcosmo  = pd.read_csv('limits/common/DP_Combined.txt', sep=' ', lineterminator='\n', names=['x', 'y'])
  plt.plot(dmcosmo['x'],   dmcosmo['y'],  color=myLighterBlue, lw=1)
  plt.fill_between(dmcosmo['x'], dmcosmo['y'], 1e-1, edgecolor='none', facecolor=myLightestBlue)

  stellar  = pd.read_csv('limits/common/DP_Combined_Stellar.txt', sep=' ',lineterminator='\n', names=['x', 'y'])
  plt.plot(stellar['x'],   stellar['y'],  color=myLightBlue, lw=1, zorder=-1)
  plt.fill_between(stellar['x'], stellar['y'], 1e-1, edgecolor='none', facecolor=myLighterBlue)
  
  #lab      = pd.read_csv('limits/common/DP_Combined_Laboratory.txt', sep=' ', lineterminator='\n', names=['x', 'y'])
  #plt.plot(lab['x'], lab['y'], color=myLightBlue, lw=1)
  #plt.fill_between(lab['x'], lab['y'], 1e-1, edgecolor='none', facecolor=myLightestBlue)

  # Default values
  Adish = 10. # dish area in m^2
  snr   = 5.  # signal to noise
  effic = 0.5 # signal detection efficiency
  time  = 24000.  # integration time in hours

  #-----------------------------
  # Projections
  #-----------------------------
  x, y = calc_darkPhoton_coupling(1e-21, Adish, 0.05, 500, snr, effic, time)
  plt.plot(x, y, alpha=0.9,lw=4, ls='-', c=myMediumOrange, zorder=4)   
  x, y = calc_darkPhoton_coupling_dcr(1e-6, Adish, 100, 4000, snr, effic, time)
  plt.plot(x, y, alpha=0.9,lw=4, ls='--', c=myMediumOrange, zorder=4) 

  x, y = calc_darkPhoton_coupling(1e-20, Adish, 0.05, 500, snr, effic, time)
  plt.plot(x, y, alpha=0.9,lw=4, ls='-', c=myDarkRed, zorder=4) 
  x, y = calc_darkPhoton_coupling_dcr(1e-4, Adish, 100, 4000, snr, effic, time)
  plt.plot(x, y, alpha=0.9,lw=4, ls='--', c=myDarkRed, zorder=4) 

  x, y = calc_darkPhoton_coupling(1e-17, Adish, 0.05, 500, snr, effic, time)
  plt.plot(x, y, alpha=0.9,lw=4, ls='-', c=myDarkPurple, zorder=4) 
  x, y = calc_darkPhoton_coupling_dcr(1e-1, Adish, 100, 4000, snr, effic, time)
  plt.plot(x, y, alpha=0.9,lw=4, ls='--', c=myDarkPurple, zorder=4) 

  plt.plot([1, 1], [1, 1], lw=4, ls='-',  c=myDarkGray, label=r'$\mathrm{Bolometer~NEP}$') 
  plt.plot([1, 1], [1, 1], lw=4, ls='--', c=myDarkGray, label=r'$\mathrm{Photocounter~DCR}$') 
  plt.legend(loc='lower right', prop={'size':22}, frameon=False, handlelength=2.4, borderpad=0.8)

  # Existing constraint labels
  text_size = 23
  fig.text(0.18, 0.50, r'Haloscope',       color=myDarkerBlue, size=text_size)
  fig.text(0.35, 0.65, r'Dish',            color=myDarkerBlue, size=text_size)
  fig.text(0.19, 0.76, r'Cosmology',       color=myDarkBlue, size=text_size)
  fig.text(0.84, 0.68, r'Stellar',         color=myDarkBlue, size=text_size)
  #fig.text(0.23, 0.78, r"$\gamma \to A'$", color=myDarkBlue, size=text_size)

  fig.text(0.39, 0.525, r'$10^{-17}\,\mathrm{W}/\sqrt{\mathrm{Hz}}$ ',    color=myDarkPurple,   size=text_size)
  fig.text(0.39, 0.42,  r'$10^{-20}\,\mathrm{W}/\sqrt{\mathrm{Hz}}$ ',    color=myDarkRed,      size=text_size)
  fig.text(0.39, 0.325, r'$10^{-21}\,\mathrm{W}/\sqrt{\mathrm{Hz}}$ ',    color=myMediumOrange, size=text_size)

  fig.text(0.80, 0.43,  r'$0.1\,\mathrm{Hz}$ ',     color=myDarkPurple,   size=text_size)
  fig.text(0.66, 0.33,  r'$10^{-4}\,\mathrm{Hz}$ ', color=myDarkRed,      size=text_size)
  fig.text(0.66, 0.29,  r'$10^{-6}\,\mathrm{Hz}$ ', color=myMediumOrange, size=text_size)
  
  # Arrows for SHUKET and Tokyo dish antenna
  ax.annotate('', xy=(0.22, 0.74),  xycoords='axes fraction',
            xytext=(0.26, 0.74), textcoords='axes fraction',
            arrowprops=dict(facecolor=myDarkerBlue, width=1.2, shrink=0.05, lw=0),
            horizontalalignment='left', verticalalignment='top',)
  ax.annotate('', xy=(0.312, 0.815),  xycoords='axes fraction',
            xytext=(0.312, 0.76), textcoords='axes fraction',
            arrowprops=dict(facecolor=myDarkerBlue, width=1.2, shrink=0.05, lw=0),
            horizontalalignment='left', verticalalignment='top',)

  # Axis log scale
  ax.set_xscale('log')
  ax.set_yscale('log')
  # Axis limits
  ax.set_xlim(1e-6, 4.13567)
  ax.set_ylim(1e-18, 1e-8)
  # Axis labels
  x_txt = r"$m_{A'}~[\mathrm{meV}]$"
  y_txt = r'$\kappa$'
  # Rescale axis from eV to meV
  scale_x=1e3
  ticks_x = mplt.ticker.FuncFormatter(lambda x, pos: r'${0:g}$'.format(x*scale_x))
  ax.xaxis.set_major_formatter(ticks_x)
  # Axis label properties
  plt.xlabel(x_txt, labelpad=15, size=38) 
  plt.ylabel(y_txt, labelpad=5, size=38)

  # Draw the upper THz axis
  ax2 = ax.twiny()
  ax2.set_xlim(241.79e-6, 1000)
  ax2.set_xscale('log')
  ax2.set_xlabel(r'$\nu~[\mathrm{THz}]$', labelpad=18, size=38)
  ax2.get_xaxis().set_major_formatter(mplt.ticker.ScalarFormatter())
  ax2.get_xaxis().set_major_formatter(mplt.ticker.FormatStrFormatter(r'$%g$'))

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

  fig.text(0.16, 0.28, r'\textbf{BREAD}', color=myDarkGray, size=text_size*1.1)
  fig.text(0.16, 0.24, r'$A_\mathrm{dish} = ' + '{0}'.format(int(Adish)) + '~\mathrm{m}^2$', color=myDarkGray, size=text_size)
  
  if time == 8760:
    integrationT = ', $\Delta t_\mathrm{int} = 1~\mathrm{day}$'
  elif time == 87600:
    integrationT = ', $\Delta t_\mathrm{int} = 10~\mathrm{days}$'
  else:
    integrationT = ', $\Delta t_\mathrm{int}' + ' = {0:.0f}'.format(time) + '~\mathrm{days}$'
  #fig.text(0.40, 0.175, r'$\mathrm{SNR}' + ' = {0:.0f}$'.format(snr) + ', $\epsilon_\mathrm{sig}' + ' = {0}$'.format(effic) + integrationT, color=myDarkGray, size=text_size*0.68)
  fig.text(0.16, 0.20, r'$\mathrm{SNR}' + ' = {0:.0f}$'.format(snr) + ', $\epsilon_\mathrm{sig}' + ' = {0}$'.format(effic) + r', $\Delta t = 1000~\mathrm{days}$', color=myDarkGray, size=text_size*0.9)

  fig.text(0.16, 0.17, r'Photocounter: SNR $\rightarrow Z = S/\sqrt{N}$', color=myMediumGray, size=text_size*0.4)
 
  # Plot margins
  plt.tight_layout(pad=0.3)
  plt.subplots_adjust( top=0.85,left=0.13 )

  # Save plot to pdf
  plt.savefig('fig_dark_photon_general_varyNEP1000days.pdf', format='pdf')

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
  rho      = ( 0.45 / relicDensity )

  kineticMixingSq = 3.7 * ratio * noise * area * dt * epsilon * rho
  kineticMixing   = math.sqrt( kineticMixingSq ) * 1e-14
 
  return [minMass*1e-3, maxMass*1e-3], [kineticMixing, kineticMixing]

#__________________________________________
def calc_darkPhoton_coupling_dcr(dcr, mirrorArea, minMass, maxMass, Zsignif=5., effic=0.5, time=1., relicDensity = 0.45):
  '''
  Convert instrument parameters and detected signal power to dark photon coupling
   - dcr          = overall dark count rate of photonsensor in Hz
   - mirrorArea   = area of dish antenna units is m^2
   - min/maxMass  = min and max mass to plot (input in meV, output in eV)
   - snr          = required signal to noise ratio 
   - effic        = overall signal power efficiency 
   - time         = integration time in hours
   - relicDensity = dark matter relic density in GeV/cm^3
  '''

  ratio    = ( Zsignif / 5. )
  noise    = math.sqrt( dcr / 0.01 )
  area     = ( 10.0 / mirrorArea  ) 
  dt       = math.sqrt( 1. / time )
  epsilon  = ( 0.5 / effic )
  rho      = ( 0.45 / relicDensity )

  kineticMixingSqMin = 5.93 * ratio * noise * area * dt * epsilon * rho * minMass
  kineticMixingSqMax = 5.93 * ratio * noise * area * dt * epsilon * rho * maxMass
  kineticMixingMin   = math.sqrt( kineticMixingSqMin ) * 1e-15
  kineticMixingMax   = math.sqrt( kineticMixingSqMax ) * 1e-15
 
  return [minMass*1e-3, maxMass*1e-3], [kineticMixingMin, kineticMixingMax]

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