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

from common import *

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
  #haloscope= pd.read_csv('limits/common/DP_Combined_Haloscope.txt',sep=' ', lineterminator='\n', names=['x', 'y'])
  #plt.fill_between(haloscope['x'], haloscope['y'], 1e-1, edgecolor='none', facecolor=myLightBlue)

  # Axion haloscopes following https://github.com/cajohare/AxionLimits/blob/master/DarkPhoton.ipynb
  costh = math.sqrt(0.019)
  dat= pd.read_csv('limits/common/Ax_ADMX.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(dat['x'], ax2dp(dat['x'], dat['y'], 7.6, costh), 1e-1, edgecolor='none', facecolor=myLightBlue)
  dat= pd.read_csv('limits/common/Ax_ADMX2018.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(dat['x'], ax2dp(dat['x'], dat['y'], 6.8, costh), 1e-1, edgecolor='none', facecolor=myLightBlue)
  dat= pd.read_csv('limits/common/Ax_ADMX2019_1.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(dat['x'], ax2dp(dat['x'], dat['y'], 7.6, costh), 1e-1, edgecolor='none', facecolor=myLightBlue)
  dat= pd.read_csv('limits/common/Ax_ADMX2019_2.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(dat['x'], ax2dp(dat['x'], dat['y'], 7.6, costh), 1e-1, edgecolor='none', facecolor=myLightBlue)
  dat= pd.read_csv('limits/common/Ax_ADMX2021.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(dat['x'], ax2dp(dat['x'], dat['y'], 7.6, costh), 1e-1, edgecolor='none', facecolor=myLightBlue)

  dat= pd.read_csv('limits/common/Ax_HAYSTAC_highres.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(dat['x'], ax2dp(dat['x'], dat['y'], 9., costh), 1e-1, edgecolor='none', facecolor=myLightBlue)
  dat= pd.read_csv('limits/common/Ax_HAYSTAC_2020_highres.txt',sep=' ', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(dat['x'], ax2dp(dat['x'], dat['y'], 9., costh), 1e-1, edgecolor='none', facecolor=myLightBlue)

  dat= pd.read_csv('limits/common/Ax_CAPP-1.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(dat['x'], ax2dp(dat['x'], dat['y'], 7.3, costh), 1e-1, edgecolor='none', facecolor=myLightBlue)
  dat= pd.read_csv('limits/common/Ax_CAPP-2.txt',sep=' ', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(dat['x'], ax2dp(dat['x'], dat['y'], 7.8, costh), 1e-1, edgecolor='none', facecolor=myLightBlue)
  dat= pd.read_csv('limits/common/Ax_CAPP-3.txt',sep=' ', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(dat['x'], ax2dp(dat['x'], dat['y'], 7.9, costh), 1e-1, edgecolor='none', facecolor=myLightBlue)

  dat = pd.read_csv('limits/common/Ax_QUAX.txt',sep=' ', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(dat['x'], ax2dp(dat['x'], dat['y'], 8.1, costh), 1e-1, edgecolor='none', facecolor=myLightBlue)
  
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
  
  lab      = pd.read_csv('limits/common/DP_Combined_Laboratory.txt', sep=' ', lineterminator='\n', names=['x', 'y'])
  plt.plot(lab['x'], lab['y'], color=myLightBlue, lw=1)
  plt.fill_between(lab['x'], lab['y'], 1e-1, edgecolor='none', facecolor=myLightestBlue)

  # Default values
  do_axion = False
  reldens  = 0.45 # relic density GeV/cm^3
  Adish    = 10. # dish area in m^2
  snr      = 5.  # signal to noise
  effic    = 0.5 # signal detection efficiency
  time     = 240.  # integration time in hours

  #-----------------------------
  # 10x time of Upgrade 1
  #-----------------------------
  # Bolometer [1.65, 83] meV
  x, y = calc_coupling_nep(5e-16, Adish, 1., 0.24, 248, snr, effic, time*100, reldens, do_axion)
  plt.plot(x, y, alpha=0.4,lw=3, ls='--', c=myDarkGray, zorder=4) 
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myMediumGray, alpha=0.05)

  # TES/KID [0.2, 125] meV
  x, y = calc_coupling_nep(2e-21, Adish, 1., 0.2, 125, snr, effic, time*100, reldens, do_axion)
  plt.plot(x, y, alpha=0.4,lw=3, ls='--', c=myDarkPurple, zorder=4) 
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myDarkPurple, alpha=0.05)

  # SNSPD [207, 830] meV
  x, y = calc_coupling_dcr(1e-8, Adish, 1., 124, 830, snr, effic, time*100, reldens, do_axion)
  plt.plot(x, y, alpha=0.4,lw=3, ls='--', c=myDarkPink, zorder=4) 
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myDarkPink, alpha=0.05)

  # QCD [2, 125] meV
  x, y = calc_coupling_dcr(4e-4, Adish, 1., 2, 125, snr, effic, time*100, reldens, do_axion)
  plt.plot(x, y, alpha=0.4,lw=3, ls='--', c=myMediumOrange, zorder=4)
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myLightOrange, alpha=0.1, zorder=4)

  #-----------------------------
  # Upgrade 1
  #-----------------------------
  # Bolometer [1.65, 83] meV
  x, y = calc_coupling_nep(5e-14, Adish, 1., 0.24, 248, snr, effic, time*100, reldens, do_axion)
  plt.plot(x, y, alpha=0.6,lw=3, ls='-.', c=myDarkGray, zorder=4) 
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myMediumGray, alpha=0.05, zorder=4)

  # TES/KID [0.2, 125] meV
  x, y = calc_coupling_nep(2e-19, Adish, 1., 0.2, 125, snr, effic, time*100, reldens, do_axion)
  plt.plot(x, y, alpha=0.6, lw=3, ls='-.',c=myDarkPurple, zorder=4) 

  # SNSPD [207, 830] meV
  x, y = calc_coupling_dcr(1e-4, Adish, 1., 124, 830, snr, effic, time*100, reldens, do_axion)
  plt.plot(x, y, alpha=0.6,lw=3, ls='-.', c=myDarkPink, zorder=4) 
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myDarkPink, alpha=0.05, zorder=4)
  
  # QCD [2, 125] meV
  x, y = calc_coupling_dcr(4., Adish, 1., 2, 125, snr, effic, time*100, reldens, do_axion)
  plt.plot(x, y, alpha=0.6,lw=3, ls='-.', c=myMediumOrange, zorder=4)
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myLightOrange, alpha=0.1, zorder=4)

  #-----------------------------
  # Baseline
  #-----------------------------
  # Far-IR 1.6 K IR Labs bolometer
  x, y = calc_coupling_nep(5e-14, Adish, 1., 0.24, 248, snr, effic, time, reldens, do_axion)
  plt.plot(x, y, alpha=0.9,lw=3, c=myDarkGray, zorder=4) 
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myMediumGray, alpha=0.05, zorder=4)

  # TES/KID [0.2, 125] meV
  x, y = calc_coupling_nep(2e-19, Adish, 1., 0.2, 125, snr, effic, time, reldens, do_axion)
  plt.plot(x, y, alpha=0.9,lw=3, c=myDarkPurple, zorder=4) 
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myDarkPurple, alpha=0.05, zorder=4)

  # SNSPD [207, 830] meV
  x, y = calc_coupling_dcr(1e-4, Adish, 1., 124, 830, snr, effic, time, reldens, do_axion)
  plt.plot(x, y, alpha=0.9,lw=3, c=myDarkPink, zorder=4) 
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myLightPink, alpha=0.05, zorder=4)

  # QCD [2, 125] meV (3e-21W/sqrt(Hz) NEP = 3.7 Hz DCR)
  x, y = calc_coupling_dcr(4., Adish, 1., 2, 125, snr, effic, time, reldens, do_axion)
  plt.plot(x, y, alpha=0.9,lw=3, c=myMediumOrange, zorder=4) 
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myLightOrange, alpha=0.1, zorder=4)

  #-----------------------------
  # Pilot dish size and 1 day 
  #-----------------------------

  # Gentec [0.4, 120] meV
  x, y = calc_coupling_nep(1e-8, Adish*(0.7/10), 1, 0.4, 120, snr, effic, time/10.0, reldens, do_axion)
  plt.plot(x, y, lw=1, ls='-', c=myDarkGreen, zorder=4) 
  plt.fill_between(x, y, 1e-1, edgecolor='none', facecolor=myMediumGreen, alpha=0.2, zorder=4)

  # Far-IR 1.6 K IR Labs bolometer
  x, y = calc_coupling_nep(5e-14, Adish*(0.7/10), 1, 0.24, 248, snr, effic, time/10.0, reldens, do_axion)
  plt.plot(x, y, alpha=0.9,lw=1, c=myDarkGray, zorder=4) 

  # SNSPD [207, 830] meV
  x, y = calc_coupling_dcr(1e-4, Adish*(0.7/10), 1, 124, 830, snr, effic, time/10.0, reldens, do_axion)
  plt.plot(x, y, alpha=0.9,lw=1, c=myDarkPink, zorder=4) 

  # Existing constraint labels
  text_size = 23
  fig.text(0.18, 0.40, r'Haloscope',       color=myDarkerBlue, size=text_size)
  fig.text(0.35, 0.565,r'Dish',            color=myDarkerBlue, size=text_size)
  fig.text(0.19, 0.64, r'Cosmology',       color=myDarkBlue, size=text_size)
  fig.text(0.87, 0.63, r'Stellar',         color=myDarkBlue, size=text_size)
  fig.text(0.23, 0.78, r"$\gamma \to A'$", color=myDarkBlue, size=text_size)
  
  # Sensors labels
  fig.text(0.63,0.807,r'\textbf{Gentec}',     color=myDarkGreen, size=text_size)
  fig.text(0.63,0.782, r'(293 K commercial)', color=myDarkGreen, size=text_size*0.5)
  fig.text(0.63,0.650, r'\textbf{IR Labs}',   color=myDarkGray,  size=text_size)
  fig.text(0.63,0.629, r'(1.6 K commercial)', color=myDarkGray,  size=text_size*0.5)

  fig.text(0.775,0.273, r'\textbf{SNSPD}',     color=myDarkPink,    size=text_size, rotation=10)
  fig.text(0.42, 0.355, r'\textbf{KID/TES}',   color=myDarkPurple,  size=text_size)
  fig.text(0.57, 0.290, r'\textbf{QCDet}',     color=myMediumOrange,size=text_size, rotation=10)

  # Arrows for SHUKET and Tokyo dish antenna
  # Horizontal arrow for SHUKET
  ax.annotate('', xy=(0.22, 0.62),  xycoords='axes fraction',
            xytext=(0.26, 0.62), textcoords='axes fraction',
            arrowprops=dict(facecolor=myDarkerBlue, width=1.2, shrink=0.05, lw=0),
            horizontalalignment='left', verticalalignment='top',)
  # Vertical arrow for Tokyo
  ax.annotate('', xy=(0.312, 0.68),  xycoords='axes fraction',
            xytext=(0.312, 0.63), textcoords='axes fraction',
            arrowprops=dict(facecolor=myDarkerBlue, width=1.2, shrink=0.05, lw=0),
            horizontalalignment='left', verticalalignment='top',)

  # Axis log scale
  ax.set_xscale('log')
  ax.set_yscale('log')
  # Axis limits
  ax.set_xlim(1e-6, 4.13567)
  ax.set_ylim(1e-18, 1e-6)
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

  # Legend
  plt.plot([1, 1], [1, 1], lw=1, ls='-',  c=myDarkGray, label=r'1 day, $A_\mathrm{dish}$/14') 
  plt.plot([1, 1], [1, 1], lw=3, ls='-',  c=myDarkGray, label=r'10 days') 
  plt.plot([1, 1], [1, 1], lw=3, ls='-.', c=myDarkGray, label=r'1000 days') 
  plt.plot([1, 1], [1, 1], lw=3, ls='--', c=myDarkGray, label=r'1000 days, NEP/100') 
  plt.legend(loc='lower right', prop={'size':17}, frameon=False, handlelength=2.7, handletextpad=0.5, borderpad=0.6, ncol=2, columnspacing=0.6)

  fig.text(0.16, 0.28, r'\textbf{BREAD}', color=myDarkGray, size=text_size*1.1)
  fig.text(0.16, 0.24, r'$A_\mathrm{dish} = ' + '{0}'.format(int(Adish)) + '~\mathrm{m}^2$', color=myDarkGray, size=text_size)
  
  if time == 8760:
    integrationT = ', $\Delta t_\mathrm{int} = 1~\mathrm{yr}$'
  elif time == 87600:
    integrationT = ', $\Delta t_\mathrm{int} = 10~\mathrm{yrs}$'
  else:
    integrationT = ', $\Delta t_\mathrm{int}' + ' = {0:.0f}'.format(time) + '~\mathrm{hrs}$'
  #fig.text(0.40, 0.175, r'$\mathrm{SNR}' + ' = {0:.0f}$'.format(snr) + ', $\epsilon_\mathrm{sig}' + ' = {0}$'.format(effic) + integrationT, color=myDarkGray, size=text_size*0.68)
  fig.text(0.16, 0.20, r'$\mathrm{SNR}' + ' = {0:.0f}$'.format(snr) + ', $\epsilon_\mathrm{sig}' + ' = {0}$'.format(effic), color=myDarkGray, size=text_size*0.9)

  fig.text(0.16, 0.17, r'SNSPD/QCDet: NEP $\to$ DCR, SNR $\rightarrow Z = S/\sqrt{N}$', color=myMediumGray, size=text_size*0.4)
 
  # Plot margins
  plt.tight_layout(pad=0.3)
  plt.subplots_adjust( top=0.85,left=0.13 )

  # Save plot to pdf
  plt.savefig('fig_dark_photon_photosensors.pdf', format='pdf')

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