#!/usr/bin/env python
'''
Axion DM sensitivity plotter
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

  # ----------------------------------------------------------
  # Input the contours to plot
  cosmological   = csv_to_lists( 'limits/dp_cosmological.csv' )
  cast           = csv_to_lists( 'limits/ax_cast.csv' )
  solar_lifetime = csv_to_lists( 'limits/dp_solar_lifetime.csv' )
  cavity1        = csv_to_lists( 'limits/ax_cavity1.csv' )
  cavity2        = csv_to_lists( 'limits/ax_cavity2.csv' )

  # Default values
  Adish  = 10. # dish area in m^2
  snr    = 5.  # signal to noise
  effic  = 0.5 # signal detection efficiency
  time   = 1.  # integration time in hours
  Bfield = 10. # Tesla


  #mk_plot( cast, cavity1, cavity2, Adish=10., Bfield=1.,  snr=5., effic=0.5, time=1. )
  mk_plot( cast, cavity1, cavity2, Adish=10.,  Bfield=10., snr=5., effic=0.5, time=24 )
  mk_plot( cast, cavity1, cavity2, Adish=10.,  Bfield=10., snr=5., effic=0.5, time=24000 )
  #mk_plot( cast, cavity1, cavity2, Adish=10., Bfield=10., snr=10.,effic=0.5, time=1. )
  #mk_plot( cast, cavity1, cavity2, Adish=10., Bfield=10., snr=5., effic=0.5, time=1. )

  '''
  
  l_time = [0.001, 0.01, 0.1, 1, 10, 100, 1e3, 8760]
  for time in l_time:
    mk_plot( cast, cavity1, cavity2, Adish, Bfield, snr, effic, time )
  time   = 10
  l_effic = [0.001, 0.01, 0.1, 0.2, 0.5, 0.9, 1.0]
  for effic in l_effic:
    mk_plot( cast, cavity1, cavity2, Adish, Bfield, snr, effic, time )
  '''

  # ----------------------------------------------------------

#__________________________________________
def mk_plot( cast, cavity1, cavity2, Adish=10., Bfield=10., snr=5., effic=0.5, time=1.):
  '''
  This plots the contours from the raw (x,y) values of each contour
  '''
  print('Plotting contours for' )
 
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

  myLightPink      = '#fcc5c0'
  myDarkPink       = '#ce1256'
  myMediumGreen    = '#41ab5d'
  myDarkGreen      = '#006d2c'

  myLightGray      = '#d9d9d9'
  myMediumGray     = '#969696'
  myDarkGray       = '#525252'

  # ------------------------------------------------------- 
  # Sensitivity lines and region
  # ------------------------------------------------------- 

  # Existing haloscopes
  plt.fill(cavity1['x'], cavity1['y'], myDarkBlue, linewidth=1, zorder=-1, edgecolor= myDarkerBlue) 
  plt.fill(cavity2['x'], cavity2['y'], myDarkBlue, linewidth=1, zorder=-1, edgecolor= myDarkerBlue)

  # QCD axion
  x  = np.array([-7.4002,0.9990])
  y1 = np.array([-15.9905,-7.5967])
  y2 = np.array([-17.6077,-9.2918])
  plt.plot([-6.561,0.9950], [-15.98,-8.415], myDarkGreen,linewidth=2,linestyle='-', zorder=-5) # KSVZ
  plt.fill_between(x, y1, y2, color=myMediumGreen,linewidth=0,alpha=0.2,zorder=-5)
  #plt.fill_between(x, y1, y2, myMediumGreen,alpha=0.6,linewidth=2,linestyle='-', zorder=-1)
  
  # HAYSTAC
  plt.plot([-4.627, -4.627], [-9, -13.6], myDarkBlue, linewidth=2, linestyle='-', zorder=-2) 
  #plt.fill(shuket['x'], shuket['y'], myLightPurple, linewidth=1, zorder=-1, edgecolor= myDarkPurple) 

  # Astro constraints
  plt.fill(cast['x'], cast['y'], myLightestBlue, linewidth=1, zorder=-1, edgecolor= myDarkBlue)
  #plt.fill(cosmological['x'], cosmological['y'], myLightBlue, linewidth=1, zorder=-1, edgecolor= myDarkerBlue)  
  #plt.fill(solar_lifetime['x'], solar_lifetime['y'], myLighterBlue, linewidth=1, zorder=-1, edgecolor= myDarkerBlue) 

  #-----------------------------
  # 100x run time of stage 1
  #-----------------------------
  # Bolometer [1.65, 83] meV
  x, y = calc_axion_coupling(1e-21, Adish, Bfield, 0.05, 1000, snr, effic, time)
  plt.plot(x, y, alpha=0.7,lw=1, ls='--', c=myDarkPurple, zorder=4) 
  plt.fill_between(x, y, [-5, -5], edgecolor='none', facecolor=myMediumPurple, alpha=0.1)

  # Bolometer [1.65, 83] meV
  x, y = calc_axion_coupling(1e-20, Adish, Bfield, 0.05, 1000, snr, effic, time)
  plt.plot(x, y, alpha=0.7,lw=3, ls='--', c=myDarkPurple, zorder=4) 
  #plt.fill_between(x, y, [-5, -5], edgecolor='none', facecolor=myDarkGray, alpha=0.1)
  
  #-----------------------------
  # Stage 1
  #-----------------------------
  # Bolometer [1.65, 83] meV
  x, y = calc_axion_coupling(1e-17, Adish, Bfield, 0.05, 1000, snr, effic, time)
  plt.plot(x, y, alpha=0.7,lw=1, ls='-', c=myDarkPurple, zorder=4) 
  #plt.fill_between(x, y, [-5, -5], edgecolor='none', facecolor=myDarkGray, alpha=0.1)

  #-----------------------------
  # Baseline
  #-----------------------------
  # Far-IR 1.6 K IR Labs bolometer
  x, y = calc_axion_coupling(1e-14, Adish, Bfield, 0.05, 1000, snr, effic, time)
  plt.plot(x, y, alpha=0.8,lw=3, ls='-', c=myDarkPurple, zorder=4) 
  #plt.fill_between(x, y, [-5, -5], edgecolor='none', facecolor=myDarkGray, alpha=0.15)

  #-----------------------------
  # General power eye-guides
  #-----------------------------

  plt.plot([1, 1], [1, 1], lw=3, ls='-',  c=myDarkPurple, label=r'NEP $= 10^{-14}$ W Hz$^{-1/2}$') 
  plt.plot([1, 1], [1, 1], lw=1, ls='-',  c=myDarkPurple, label=r'NEP $= 10^{-17}$ W Hz$^{-1/2}$') 
  plt.plot([1, 1], [1, 1], lw=3, ls='--', c=myDarkPurple, label=r'NEP $= 10^{-20}$ W Hz$^{-1/2}$') 
  plt.plot([1, 1], [1, 1], lw=1, ls='--', c=myDarkPurple, label=r'NEP $= 10^{-21}$ W Hz$^{-1/2}$') 
  plt.legend(loc='lower right', prop={'size':18}, frameon=False, handlelength=2.8, borderpad=0.8)

  # ------------------------------------------------------- 
  # Axis properties
  # ------------------------------------------------------- 

  plt.xlim(-6, 1)
  plt.ylim(-16, -6)

  # axes labels
  x_txt = r"$\log_{10}[m_{a} / \mathrm{eV}]$"
  y_txt = r'$\log_{10}[|g_{a\gamma\gamma}| / \mathrm{GeV}^{-1}]$'
    
  # Axis label properties
  plt.xlabel(x_txt, labelpad=20, size=35) 
  plt.ylabel(y_txt, labelpad=10, size=35)

  # Draw the upper THz axis
  ax2 = ax.twiny()
  ax2.set_xlim(-3.62, 3.38)
  ax2.set_xlabel(r'$\log_{10}[\nu/\mathrm{THz}]$', labelpad=10, size=35)
   
  # ------------------------------------------------------- 
  # Add plot text
  # ------------------------------------------------------- 
  text_size = 23

  #fig.text(0.22, 0.20, r'ADMX',    color=myDarkerBlue, size=text_size)
  fig.text(0.22, 0.285, r'Cavity', color=myDarkerBlue, size=text_size)
  #fig.text(0.34, 0.62, r'SHUKET', color=myMediumPurple, size=text_size, rotation=90)
  fig.text(0.325, 0.52, r'HAYSTAC',color=myDarkBlue, size=text_size, rotation=90)

  #fig.text(0.22, 0.76, r'Cosmology', color=myDarkerBlue, size=text_size)
  fig.text(0.25, 0.74, r'CAST', color=myDarkerBlue, size=text_size)
  fig.text(0.68, 0.56, r'KSVZ', color=myDarkGreen, size=text_size, rotation=28)
  fig.text(0.74, 0.612, r'QCD axion models', color=myMediumGreen, size=text_size, rotation=26)

  # Sensors
  #fig.text(0.62, 0.78, r'Bolometer',    color=myDarkGray,    size=text_size)
  #fig.text(0.62, 0.76, r'(Commercial)', color=myDarkGray,    size=text_size*0.5)
  #fig.text(0.80, 0.80, r'SNSPD',     color=myDarkPink,    size=text_size)
  #fig.text(0.54, 0.61, r'KID',       color=myMediumOrange,size=text_size)
  #fig.text(0.43, 0.56, r'TES',       color=myDarkPurple,  size=text_size)
  #fig.text(0.61, 0.63, r'QCD',       color=myDarkGray,    size=text_size, rotation=90)

  # Eye guides
  #fig.text(0.35, 0.65, r'$10^{-15}~\mathrm{W}$', color=myMediumGray, size=0.7*text_size, rotation=25)
  #fig.text(0.35, 0.47, r'$10^{-20}~\mathrm{W}$', color=myMediumGray, size=0.7*text_size, rotation=25)
  #fig.text(0.35, 0.30, r'$10^{-25}~\mathrm{W}$', color=myMediumGray, size=0.7*text_size, rotation=25)
  fig.text(0.70, 0.39, r'\textbf{BREAD}', color=myDarkGray, size=text_size)
  fig.text(0.70, 0.35, r'$A_\mathrm{dish} = ' + '{0}'.format(int(Adish)) + '~\mathrm{m}^2$'r', $B = ' + '{0}'.format(int(Bfield)) + '~\mathrm{T}$', color=myDarkGray, size=text_size*0.8)
  
  if time == 24:
    integrationT = '$\Delta t_\mathrm{int} = 1~\mathrm{days}$'
  elif time == 24000:
    integrationT = '$\Delta t_\mathrm{int} = 1000~\mathrm{days}$'
  else:
    integrationT = '$\Delta t_\mathrm{int}' + ' = {0:.0f}'.format(time) + '~\mathrm{hrs}$'
  fig.text(0.40, 0.22, r'$\mathrm{SNR}' + ' = {0:.0f}$'.format(snr) + r', $\epsilon_\mathrm{sig}' + ' = {0}$'.format(effic), color=myDarkGray, size=text_size*0.8)
  fig.text(0.40, 0.18, integrationT, color=myDarkGray, size=text_size*0.8)
  #fig.text(0.49, 0.18, r'$\epsilon_\mathrm{sig}' + ' = {0}$'.format(effic), color=myDarkGray, size=text_size*0.8)


  # Adjust axis ticks
  ax.minorticks_on()
  ax2.minorticks_on()

  ax.tick_params('x', length=12, width=1, which='major', labelsize='28', pad=10)
  ax.tick_params('x', length=6,  width=1, which='minor') 
  ax2.tick_params('x', length=12, width=1, which='major', labelsize='28', pad=10)
  ax2.tick_params('x', length=6,  width=1, which='minor') 
  ax.tick_params('y', length=12, width=1, which='major', labelsize='28', pad=10)
  ax.tick_params('y', length=6,  width=1, which='minor') 
   
  plt.tight_layout(pad=0.3)
  plt.subplots_adjust( top=0.85,left=0.15 )
  
  save_name = 'fig_axion_Adish{0}_Bfield{1}_snr{2}_effic{3}_time{4}'.format(Adish, Bfield, snr, effic, time)
  
  print('Saving as {0}'.format(save_name))
  plt.savefig(save_name + '_general.pdf', format='pdf', dpi=50)
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

#__________________________________________
def calc_axion_coupling(nep, mirrorArea, Bfield, minMass, maxMass, snr=5., effic=0.5, time=1., relicDensity = 0.45):
  '''
  Convert instrument parameters and detected signal power to axion coupling
   - power is units of Watts
   - mirrorArea units is m^2
   - Bfield units is Tesla
   - min/maxMass  = min and max mass to plot
   - snr          = required signal to noise ratio 
   - effic        = overall signal power efficiency 
   - time         = integration time in hours
   - relicDensity = dark matter relic density in GeV/cm^3
  Returns lowest [minCoupling, maxCoupling] coupling values probed 
  in units of 10^{-11}/GeV for [minMass, maxMass] 
  '''
  ratio    = ( snr / 5. )
  noise    = ( nep / 1.0e-21 )
  area     = ( 10. / mirrorArea  ) 
  dt       = math.sqrt( 1. / time )
  epsilon  = ( 0.5 / effic )
  rho      = ( 0.3 / relicDensity )
  magnet   = ( 10. / Bfield )**2 

  massOverCouplingSq = 1. / ( 5.376 * ratio * noise * area * dt * epsilon * rho * magnet )
  massOverCoupling = math.sqrt( massOverCouplingSq )

  # For the input min and max masses, calculate corresponding couplings
  minCoupling = ( minMass ) / massOverCoupling
  maxCoupling = ( maxMass ) / massOverCoupling

  # Convert meV to eV, coupling to 1/GeV
  log10minMass = np.log10(minMass * 10**(-3))
  log10maxMass = np.log10(maxMass * 10**(-3))
  log10minCoupling = np.log10(minCoupling * 10**(-12)) 
  log10maxCoupling = np.log10(maxCoupling * 10**(-12)) 

  return [log10minMass, log10maxMass], [log10minCoupling, log10maxCoupling]
  
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