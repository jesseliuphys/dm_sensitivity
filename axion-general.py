#!/usr/bin/env python3
'''
Axion DM sensitivity plotter
for a general ideal broadband sensor
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

 
  # Figures
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

  myLightGreen     = '#e5f5e0'
  myMediumGreen    = '#41ab5d'
  myDarkGreen      = '#006d2c'
  myDarkerGreen    = '#00441b'

  myLightGray      = '#d9d9d9'
  myMediumGray     = '#969696'
  myDarkGray       = '#525252'

  # Existing constraints from https://arxiv.org/abs/2105.04565 https://github.com/cajohare/AxionLimits/ 
  # Haloscope
  admx= pd.read_csv('limits/common/Ax_ADMX.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(admx['x'], admx['y'], 1e-1, edgecolor='none', facecolor=myLightBlue)
  admx1= pd.read_csv('limits/common/Ax_ADMX2018.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(admx1['x'], admx1['y'], 1e-1, edgecolor='none', facecolor=myLightBlue)
  admx2= pd.read_csv('limits/common/Ax_ADMX2019_1.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(admx2['x'], admx2['y'], 1e-1, edgecolor='none', facecolor=myLightBlue)
  admx3= pd.read_csv('limits/common/Ax_ADMX2019_2.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(admx3['x'], admx3['y'], 1e-1, edgecolor='none', facecolor=myLightBlue)  
  admx4= pd.read_csv('limits/common/Ax_ADMX2021.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(admx4['x'], admx4['y'], 1e-1, edgecolor='none', facecolor=myLightBlue)

  capp1= pd.read_csv('limits/common/Ax_CAPP-1.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(capp1['x'], capp1['y'], 1e-1, edgecolor='none', facecolor=myLightBlue)
  capp2= pd.read_csv('limits/common/Ax_CAPP-2.txt',sep=' ', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(capp2['x'], capp2['y'], 1e-1, edgecolor='none', facecolor=myLightBlue)
  capp3= pd.read_csv('limits/common/Ax_CAPP-3.txt',sep=' ', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(capp3['x'], capp3['y'], 1e-1, edgecolor='none', facecolor=myLightBlue)

  haystac= pd.read_csv('limits/common/Ax_HAYSTAC_highres.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(haystac['x'], haystac['y'], 1e-1, edgecolor='none', facecolor=myLightBlue)
  haystac2= pd.read_csv('limits/common/Ax_HAYSTAC_2020_highres.txt',sep=' ', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(haystac2['x'], haystac2['y'], 1e-1, edgecolor='none', facecolor=myLightBlue)

  rbf = pd.read_csv('limits/common/Ax_RBF_UF_Haloscopes.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(rbf['x'], rbf['y'], 1e-1, edgecolor='none', facecolor=myLightBlue)

  organ = pd.read_csv('limits/common/Ax_ORGAN.txt',sep=' ', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(organ['x'], organ['y'], 1e-1, edgecolor='none', facecolor=myLightBlue)

  quax = pd.read_csv('limits/common/Ax_QUAX.txt',sep=' ', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(quax['x'], quax['y'], 1e-1, edgecolor='none', facecolor=myLightBlue)
  quax2 = pd.read_csv('limits/common/Ax_QUAX2.txt',sep=' ', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(quax2['x'], quax2['y'], 1e-1, edgecolor='none', facecolor=myLightBlue)

  rades = pd.read_csv('limits/common/Ax_RADES.txt',sep=' ', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(rades['x'], rades['y'], 1e-1, edgecolor='none', facecolor=myLightBlue)

  muse = pd.read_csv('limits/common/Ax_Telescopes_MUSE.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(muse['x'], muse['y'], 1e-1, edgecolor='none', facecolor=myLighterBlue, zorder=-1, alpha=0.8)
  vimos = pd.read_csv('limits/common/Ax_Telescopes_VIMOS.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(vimos['x'], vimos['y'], 1e-1, edgecolor='none', facecolor=myLighterBlue, zorder=-1, alpha=0.8)

  # Neutron star
  ns0 = pd.read_csv('limits/common/Ax_NeutronStar.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(ns0['x'], ns0['y'], 1e-1, edgecolor='none', facecolor=myLighterBlue)
  ns1 = pd.read_csv('limits/common/Ax_NeutronStars_Battye.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(ns1['x'], ns1['y']*100, 1e-1, edgecolor='none', facecolor=myLighterBlue) #x100 following https://github.com/cajohare/AxionLimits/blob/master/AxionPhoton_Closeups.ipynb
  ns2 = pd.read_csv('limits/common/Ax_NeutronStars_GreenBank.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(ns2['x'], ns2['y'], 1e-1, edgecolor='none', facecolor=myLighterBlue)
  ns3 = pd.read_csv('limits/common/Ax_NeutronStars_VLA.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.fill_between(ns3['x'], ns3['y'], 1e-1, edgecolor='none', facecolor=myLighterBlue)

  # Stellar
  hb= pd.read_csv('limits/common/Ax_HorizontalBranch.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.plot(hb['x'],   hb['y'],  color=myLightBlue, lw=1, zorder=0)
  plt.fill_between(hb['x'], hb['y'], 1e-1, edgecolor='none', facecolor=myLighterBlue, alpha=0.8)

  cast= pd.read_csv('limits/common/Ax_CAST.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.plot(cast['x'],   cast['y'],  color=myLightBlue, lw=1, zorder=2)
  plt.fill_between(cast['x'], cast['y'], 1e-1, edgecolor='none', facecolor=myLightestBlue)
  cast2= pd.read_csv('limits/common/Ax_CAST_highm.txt',sep='\t', lineterminator='\n', names=['x', 'y'])
  plt.plot(cast2['x'],   cast2['y'],  color=myLightBlue, lw=1, zorder=-1)
  plt.fill_between(cast2['x'], cast2['y'], 1e-1, edgecolor='none', facecolor=myLightestBlue)

  # QCD axion band
  x  = np.array([3.97923878e-8,9.977])
  y1 = np.array([1.0221156e-16,2.53104578e-8])
  y2 = np.array([2.4677434e-18,5.1074015e-10])
  plt.plot([2.74789415e-7,9.88553094657], [1.0471285e-16,3.84591782e-9], myDarkerGreen,linewidth=2,linestyle='-', zorder=-5) # KSVZ
  plt.plot([2.74789415e-7,9.88553094657], [0.391*1.0471285e-16,0.391*3.84591782e-9], myDarkGreen, alpha=0.7, linewidth=2,linestyle='-', zorder=-5) # DFSZ
  
  plt.fill_between(x, y1, y2, color=myMediumGreen,linewidth=0,alpha=0.2,zorder=-5)
  #plt.fill_between(x, y1, y2, myMediumGreen,alpha=0.6,linewidth=2,linestyle='-', zorder=-1)

  # Cogenesis cayy=1
  plt.plot([1.03535e-8, 0.999511], [5.79881e-14, 5.72150e-10], color=myMediumGreen, lw=2, ls='dotted', zorder=-1)

  # Default values
  reldens = 0.45 # relic density in GeV/cm^3
  Adish   = 10. # dish area in m^2
  Bfield  = 10. # Tesla
  snr     = 5.  # signal to noise
  effic   = 0.5 # signal detection efficiency
  time    = 24000.  # integration time in hours

  #-----------------------------
  # Projections
  #-----------------------------
  do_axion = True
  x, y = calc_coupling_nep(1e-22, Adish, Bfield, 0.05, 500, snr, effic, time, reldens, do_axion)
  plt.plot(x, y, alpha=0.9,lw=4, ls='-', c=myMediumOrange, zorder=4)   
  x, y = calc_coupling_dcr(1e-6, Adish, Bfield, 4, 4000, snr, effic, time, reldens, do_axion)
  plt.plot(x, y, alpha=0.9,lw=4, ls='--', c=myMediumOrange, zorder=4) 

  x, y = calc_coupling_nep(1e-20, Adish, Bfield, 0.05, 500, snr, effic, time, reldens, do_axion)
  plt.plot(x, y, alpha=0.9,lw=4, ls='-', c=myDarkRed, zorder=4) 
  x, y = calc_coupling_dcr(1e-4, Adish, Bfield, 4, 4000, snr, effic, time, reldens, do_axion)
  plt.plot(x, y, alpha=0.9,lw=4, ls='--', c=myDarkRed, zorder=4) 

  x, y = calc_coupling_nep(1e-18, Adish, Bfield, 0.05, 500, snr, effic, time, reldens, do_axion)
  plt.plot(x, y, alpha=0.9,lw=4, ls='-', c=myDarkPurple, zorder=4) 
  x, y = calc_coupling_dcr(1e-2, Adish, Bfield, 4, 4000, snr, effic, time, reldens, do_axion)
  plt.plot(x, y, alpha=0.9,lw=4, ls='--', c=myDarkPurple, zorder=4) 

  plt.plot([1, 1], [1, 1], lw=4, ls='-',  c=myDarkGray, label=r'$\mathrm{Bolometer~NEP}$') 
  plt.plot([1, 1], [1, 1], lw=4, ls='--', c=myDarkGray, label=r'$\mathrm{Photocounter~DCR}$') 
  plt.legend(loc='lower right', prop={'size':22}, frameon=False, handlelength=2.4, borderpad=0.8)
    
  # ------------------------------------------------------- 
  # Add plot text
  # ------------------------------------------------------- 
  text_size = 23
  fig.text(0.19, 0.45,  r'Haloscope',  color=myDarkerBlue, size=text_size)
  fig.text(0.22, 0.77,  r'CAST',       color=myDarkBlue, size=text_size)
  fig.text(0.85, 0.74,  r'Stellar',    color=myDarkBlue, size=text_size)
  fig.text(0.92, 0.53,  r'Telescope',  color=myDarkBlue, size=text_size, rotation=90)

  # QCD axions
  fig.text(0.24, 0.285,r'KSVZ', color=myDarkerGreen, size=text_size, rotation=34)
  fig.text(0.24, 0.20, r'DFSZ', color=myDarkGreen, size=text_size, rotation=34)
  fig.text(0.26, 0.17, r'QCD axion models', alpha=0.8, color=myMediumGreen, size=text_size, rotation=34)

  fig.text(0.31, 0.595,  r'Cogenesis', color=myMediumGreen, size=text_size*0.9, rotation=19)
  fig.text(0.34, 0.585, r'$c_{a\gamma\gamma} = 1$', color=myMediumGreen, size=text_size*0.5, rotation=19)

  fig.text(0.45, 0.67,  r'$10^{-18}\,\mathrm{W}/\sqrt{\mathrm{Hz}}$ ',    color=myDarkPurple,   size=text_size, rotation=34)
  fig.text(0.45, 0.57,  r'$10^{-20}\,\mathrm{W}/\sqrt{\mathrm{Hz}}$ ',    color=myDarkRed,      size=text_size, rotation=34)
  fig.text(0.45, 0.47,  r'$10^{-22}\,\mathrm{W}/\sqrt{\mathrm{Hz}}$ ',    color=myMediumOrange, size=text_size, rotation=34)

  fig.text(0.64, 0.64,  r'$10^{-2}\,\mathrm{Hz}$ ', color=myDarkPurple,   size=text_size, rotation=44)
  fig.text(0.51, 0.42,  r'$10^{-4}\,\mathrm{Hz}$ ', color=myDarkRed,      size=text_size, rotation=34)
  fig.text(0.51, 0.38,  r'$10^{-6}\,\mathrm{Hz}$ ', color=myMediumOrange, size=text_size, rotation=34)
  
  fig.text(0.38, 0.17, r'Photocounter: SNR $\rightarrow Z = S/\sqrt{N}$', color=myMediumGray, size=text_size*0.4)
 
  if time == 24:
    integrationT = '$\Delta t = 1~\mathrm{day}$'
  elif time == 24000:
    integrationT = '$\Delta t = 1000~\mathrm{days}$'
  else:
    integrationT = '$\Delta t' + ' = {0:.0f}'.format(time) + '~\mathrm{hrs}$'
  # Labels
  fig.text(0.62, 0.43, r'\textbf{BREAD}', color=myDarkGray, size=text_size*1.2)
  fig.text(0.62, 0.39, r'$A_\mathrm{dish} = ' + '{0}'.format(int(Adish)) + '~\mathrm{m}^2$'r', $B_\mathrm{ext} = ' + '{0}'.format(int(Bfield)) + '~\mathrm{T}$', color=myDarkGray, size=text_size)
  fig.text(0.62, 0.34, r'$\mathrm{SNR}' + ' = {0:.0f}$'.format(snr) + r', $\epsilon_\mathrm{sig}' + ' = {0}$'.format(effic), color=myDarkGray, size=text_size)
  fig.text(0.62, 0.29, integrationT, color=myDarkGray, size=text_size)
  
  # ------------------------------------------------------- 
  # Axis properties
  # ------------------------------------------------------- 
  # Axis log scale
  ax.set_xscale('log')
  ax.set_yscale('log')
  # Axis limits
  ax.set_xlim(1e-6, 4.13567)
  ax.set_ylim(1e-16, 1e-9)
  # Axis labels
  x_txt = r"$m_{a}~[\mathrm{meV}]$"
  y_txt = r'$|g_{a\gamma\gamma}|~[\mathrm{GeV}^{-1}]$'
  # Rescale axis from eV to meV
  scale_x=1e3
  ticks_x = mplt.ticker.FuncFormatter(lambda x, pos: r'${0:g}$'.format(x*scale_x))
  ax.xaxis.set_major_formatter(ticks_x)
  # Axis label properties
  plt.xlabel(x_txt, labelpad=15, size=38) 
  plt.ylabel(y_txt, labelpad=10, size=38)

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
   
  # Plot margins
  plt.tight_layout(pad=0.3)
  plt.subplots_adjust( top=0.85,left=0.16 )
  
  save_name = 'fig_axion_general_varyNEP1000days'
  
  print('Saving as {0}'.format(save_name))
  plt.savefig(save_name + '.pdf', format='pdf', dpi=50)
  #plt.savefig(save_name + '_photosensors.png', format='png', dpi=400)
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