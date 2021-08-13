#!/usr/bin/env python
'''
Dark photon DM sensitivity plotter
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
  cast           = csv_to_lists( 'limits/dp_cast.csv' )
  solar_lifetime = csv_to_lists( 'limits/dp_solar_lifetime.csv' )
  cavity1        = csv_to_lists( 'limits/dp_cavity1.csv' )
  cavity2        = csv_to_lists( 'limits/dp_cavity2.csv' )
  cavity3        = csv_to_lists( 'limits/dp_cavity3.csv' )
  shuket         = csv_to_lists( 'limits/dp_shuket.csv' )

  # Default values
  Adish = 10. # dish area in m^2
  snr   = 5.  # signal to noise
  effic = 0.5 # signal detection efficiency
  time  = 1.  # integration time in hours

  mk_plot( cosmological, cast, solar_lifetime, cavity1, cavity2, cavity3, shuket, Adish=10., snr=5., effic=0.5, time=240. )
  #mk_plot( cosmological, cast, solar_lifetime, cavity1, cavity2, cavity3, shuket, Adish=1.,  snr=5., effic=0.5, time=1. )
  #mk_plot( cosmological, cast, solar_lifetime, cavity1, cavity2, cavity3, shuket, Adish=10., snr=10.,effic=0.5, time=1. )
  #mk_plot( cosmological, cast, solar_lifetime, cavity1, cavity2, cavity3, shuket, Adish=10., snr=5., effic=0.5, time=0.1 )

  '''
  l_time = [0.001, 0.01, 0.1, 1, 10, 100, 1e3, 8760]
  for time in l_time:
    mk_plot( cosmological, cast, solar_lifetime, cavity1, cavity2, cavity3, shuket, Adish, snr, effic, time )
  time  = 10
  l_effic = [0.001, 0.01, 0.1, 0.2, 0.5, 0.9, 1.0]
  for effic in l_effic:
    mk_plot( cosmological, cast, solar_lifetime, cavity1, cavity2, cavity3, shuket, Adish, snr, effic, time )
  '''
  # ----------------------------------------------------------

#__________________________________________
def mk_plot( cosmological, cast, solar_lifetime, cavity1, cavity2, cavity3, shuket, Adish=10., snr=5., effic=0.5, time=1. ):
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
  # Sensitivity lines and regions
  # ------------------------------------------------------- 

  # Existing haloscopes
  plt.fill(cavity1['x'], cavity1['y'], myDarkBlue, linewidth=1, zorder=-1, edgecolor= myDarkerBlue) 
  plt.fill(cavity2['x'], cavity2['y'], myDarkBlue, linewidth=1, zorder=-1, edgecolor= myDarkerBlue)
  plt.fill(cavity3['x'], cavity3['y'], myDarkBlue, linewidth=1, zorder=-1, edgecolor= myDarkerBlue)

  plt.plot([-4.6045, -4.6045], [-10, -15], myDarkBlue, linewidth=1, linestyle='-', zorder=-2) 
  plt.fill(shuket['x'], shuket['y'], myMediumBlue, linewidth=1, zorder=-1, edgecolor= myDarkerBlue) 

  # Astro constraints
  plt.fill(cast['x'],           cast['y'],           myMediumBlue,  linewidth=1, zorder=-1, edgecolor= myDarkerBlue)
  plt.fill(cosmological['x'],   cosmological['y'],   myLightestBlue,linewidth=1, zorder=-1, edgecolor= myDarkerBlue)  
  plt.fill(solar_lifetime['x'], solar_lifetime['y'], myLighterBlue, linewidth=1, zorder=-1, edgecolor= myDarkerBlue) 

  #-----------------------------
  # 10x time of Upgrade 1
  #-----------------------------
  # Bolometer [1.65, 83] meV
  x, y = calc_darkPhoton_coupling(4e-17, Adish, 0.24, 248, snr, effic, time*100)
  plt.plot(x, y, alpha=0.7,lw=3, ls='--', c=myDarkGray, zorder=4) 
  plt.fill_between(x, y, [-5, -5], edgecolor='none', facecolor=myDarkGray, alpha=0.1)

  # KID [0.2, 5] meV
  x, y = calc_darkPhoton_coupling(3e-21, Adish, 0.2, 5, snr, effic, time*100)
  plt.plot(x, y, myLightOrange, alpha=0.7,lw=3, ls='--', c=myMediumOrange, zorder=4) 
  plt.fill_between(x, y, [-5, -5], edgecolor='none', facecolor=myLightOrange, alpha=0.2)

  # TES [0.2, 1.2] meV
  x, y = calc_darkPhoton_coupling(2e-21, Adish, 0.19, 1.2, snr, effic, time*100)
  plt.plot(x, y, alpha=0.8,lw=3, ls='--', c=myDarkPurple, zorder=4) 
  plt.fill_between(x, y, [-5, -5], edgecolor='none', facecolor=myDarkPurple, alpha=0.15)

  # SNSPD [207, 830] meV
  x, y = calc_darkPhoton_coupling(1e-21, Adish, 207, 830, snr, effic, time*100)
  plt.plot(x, y, alpha=0.3,lw=3, ls='--', c=myDarkPink, zorder=4) 
  plt.fill_between(x, y, [-5, -5], edgecolor='none', facecolor=myDarkPink, alpha=0.1)

  # QCD [6.2] meV
  x, y = calc_darkPhoton_coupling(1e-22, Adish, 6.18, 6.22, snr, effic, time*100)
  plt.fill_between(x, y, [-5, -5], color=myDarkGray, alpha=0.3)

  #-----------------------------
  # Upgrade 1
  #-----------------------------
  # Bolometer [1.65, 83] meV
  x, y = calc_darkPhoton_coupling(4e-17, Adish, 0.24, 248, snr, effic, time)
  plt.plot(x, y, alpha=0.7,lw=3, ls='-.', c=myDarkGray, zorder=4) 
  plt.fill_between(x, y, [-5, -5], edgecolor='none', facecolor=myDarkGray, alpha=0.1)

  # KID [0.2, 5] meV
  x, y = calc_darkPhoton_coupling(3e-21, Adish, 0.2, 5, snr, effic, time)
  plt.plot(x, y, myLightOrange, alpha=0.7,lw=3, ls='-.', c=myMediumOrange, zorder=4) 
  plt.fill_between(x, y, [-5, -5], edgecolor='none', facecolor=myLightOrange, alpha=0.2)
  
  # TES [0.2, 1.2] meV
  x, y = calc_darkPhoton_coupling(2e-21, Adish, 0.19, 1.2, snr, effic, time)
  plt.plot(x, y, alpha=0.8, lw=3, ls='-.',c=myDarkPurple, zorder=4) 

  # SNSPD [207, 830] meV
  x, y = calc_darkPhoton_coupling(1e-21, Adish, 207, 830, snr, effic, time)
  plt.plot(x, y, alpha=0.3,lw=3, ls='-.', c=myDarkPink, zorder=4) 
  plt.fill_between(x, y, [-5, -5], edgecolor='none', facecolor=myDarkPink, alpha=0.1)

  #-----------------------------
  # Baseline
  #-----------------------------
  '''
  # Standard 1.6 K IR Labs bolometer [0.24, 83] meV
  x, y = calc_darkPhoton_coupling(1e-13, Adish, 83.0, 248.0, snr, effic, time)
  plt.plot(x, y, alpha=0.8,lw=3, c=myDarkGray, zorder=4) 
  plt.fill_between(x, y, [-5, -5], edgecolor='none', facecolor=myDarkGray, alpha=0.15)

  # Standard 1.6 K IR Labs bolometer [0.24, 83] meV
  x, y = calc_darkPhoton_coupling(5e-14, Adish, 4.0, 83.0, snr, effic, time)
  plt.plot(x, y, alpha=0.8,lw=3, c=myDarkGray, zorder=4) 
  plt.fill_between(x, y, [-5, -5], edgecolor='none', facecolor=myDarkGray, alpha=0.15)
  '''
  # Far-IR 1.6 K IR Labs bolometer
  x, y = calc_darkPhoton_coupling(4e-15, Adish, 0.24, 248, snr, effic, time)
  plt.plot(x, y, alpha=0.8,lw=3, c=myDarkGray, zorder=4) 
  plt.fill_between(x, y, [-5, -5], edgecolor='none', facecolor=myDarkGray, alpha=0.15)

  # KID [0.2, 5] meV
  x, y = calc_darkPhoton_coupling(3e-19, Adish, 0.2, 5, snr, effic, time)
  plt.plot(x, y, alpha=0.8,lw=3, c=myMediumOrange, zorder=4) 
  plt.fill_between(x, y, [-5, -5], edgecolor='none', facecolor=myLightOrange, alpha=0.6)

  # TES [0.2, 1.2] meV
  x, y = calc_darkPhoton_coupling(2e-19, Adish, 0.19, 1.2, snr, effic, time)
  plt.plot(x, y, alpha=0.8,lw=3, c=myDarkPurple, zorder=4) 
  plt.fill_between(x, y, [-5, -5], edgecolor='none', facecolor=myDarkPurple, alpha=0.15)

  # SNSPD [207, 830] meV
  x, y = calc_darkPhoton_coupling(1e-18, Adish, 207, 830, snr, effic, time)
  plt.plot(x, y, alpha=0.8,lw=3, c=myDarkPink, zorder=4) 
  plt.fill_between(x, y, [-5, -5], edgecolor='none', facecolor=myLightPink, alpha=0.6)

  # QCD [6.2] meV
  x, y = calc_darkPhoton_coupling(1e-20, Adish, 6.18, 6.22, snr, effic, time)
  #plt.plot(x, y, myMediumGray, alpha=0.8,lw=3, c=myDarkGray, zorder=4, label=r'$\mathrm{QCD}~10^{-20}~\mathrm{W}~\mathrm{Hz}^{-1/2}$') 
  plt.fill_between(x, y, [-5, -5], color=myDarkGray, alpha=0.9)

  # Gentec [0.4, 120] meV
  x, y = calc_darkPhoton_coupling(1e-8, Adish, 0.4, 120, snr, effic, time)
  plt.plot(x, y, lw=3, ls='-', c=myDarkGreen) 
  plt.fill_between(x, y, [-5, -5], edgecolor='none', facecolor=myMediumGreen, alpha=0.3)

  #-----------------------------
  # Pilot dish size and 1 day 
  #-----------------------------

  # Far-IR 1.6 K IR Labs bolometer
  x, y = calc_darkPhoton_coupling(4e-15, Adish/10.0, 0.24, 248, snr, effic, time/10.0)
  plt.plot(x, y, alpha=0.8,lw=1, c=myDarkGray, zorder=4) 

  # SNSPD [207, 830] meV
  x, y = calc_darkPhoton_coupling(1e-18, Adish/10.0, 207, 830, snr, effic, time/10.0)
  plt.plot(x, y, alpha=0.8,lw=1, c=myDarkPink, zorder=4) 

  #-----------------------------
  # General power eye-guides
  #-----------------------------
  #x, y = calc_darkPhoton_coupling(1e-15, Adish, 1e-6, 1e4, snr=1., effic=1., time=1/3600.)
  #plt.plot(x, y, lw=1, ls=':', c=myMediumGray) 
  
  #x, y = calc_darkPhoton_coupling(1e-20, Adish, 1e-6, 1e4, snr=1., effic=1., time=1/3600.)
  #plt.plot(x, y, lw=1, ls=':', c=myMediumGray) 

  #x, y = calc_darkPhoton_coupling(1e-25, Adish, 1e-6, 1e4, snr=1., effic=1., time=1/3600.)
  #plt.plot(x, y, lw=1, ls=':', c=myMediumGray) 

  # Constant rate
  #x, y = calc_darkPhoton_coupling_rate(3600.0, Adish, 1e-6, 1e4)
  #plt.plot(x, y, lw=1, ls=':', c=myMediumGray) 

  #x, y = calc_darkPhoton_coupling_rate(1.0, Adish, 1e-6, 1e4)
  #plt.plot(x, y, lw=1, ls=':', c=myMediumGray) 

  #plt.plot([1, 1], [1, 1], lw=3, ls='-',  c=myDarkGray, label='Baseline') 
  #plt.plot([1, 1], [1, 1], lw=3, ls='-.', c=myDarkGray, label='Upgrade 1') 
  #plt.plot([1, 1], [1, 1], lw=3, ls='--', c=myDarkGray, label='Upgrade 2') 
  plt.plot([1, 1], [1, 1], lw=1, ls='-',  c=myDarkGray, label=r'NEP$_\mathsf{today}$, 1 day, $A_\mathrm{dish}/10$') 
  plt.plot([1, 1], [1, 1], lw=3, ls='-',  c=myDarkGray, label=r'NEP$_\mathsf{today}$, 10 days') 
  plt.plot([1, 1], [1, 1], lw=3, ls='-.', c=myDarkGray, label=r'NEP$_\mathsf{today}/100$, 10 days') 
  plt.plot([1, 1], [1, 1], lw=3, ls='--', c=myDarkGray, label=r'NEP$_\mathsf{today}/100$, 1000 days') 
  plt.legend(loc='lower right', prop={'size':16}, frameon=False, handlelength=2.8, borderpad=0.8)

  # ------------------------------------------------------- 
  # Axis properties
  # ------------------------------------------------------- 

  plt.xlim(-6, 1)
  plt.ylim(-18, -6.5)

  # axes labels
  x_txt = r"$\log_{10}[m_{A'} / \mathrm{eV}]$"
  y_txt = r'$\log_{10}\kappa$'
    
  # Axis label properties
  plt.xlabel(x_txt, labelpad=20, size=35) 
  plt.ylabel(y_txt, labelpad=20, size=35)

  # Draw the upper THz axis
  ax2 = ax.twiny()
  ax2.set_xlim(-3.62, 3.38)
  ax2.set_xlabel(r'$\log_{10}[\nu/\mathrm{THz}]$', labelpad=10, size=35)
   
  # ------------------------------------------------------- 
  # Add plot text
  # ------------------------------------------------------- 
  text_size = 23

  # Cavity haloscopes
  fig.text(0.18, 0.25,  r'Cavity',  color=myDarkerBlue, size=text_size)
  #fig.text(0.18, 0.22, r'ADMX,',  color=myDarkerBlue, size=text_size)
  #fig.text(0.18, 0.18, r'RBF, UF',color=myDarkerBlue, size=text_size)
  fig.text(0.34, 0.59,  r'SHUKET', color=myDarkBlue,   size=text_size, rotation=90)
  fig.text(0.325, 0.39, r'Qubit',  color=myDarkBlue,   size=text_size, rotation=90)

  # Astrophysics
  fig.text(0.22, 0.76, r'Cosmology', color=myDarkerBlue, size=text_size)
  fig.text(0.86, 0.68, r'CAST',      color=myDarkerBlue, size=text_size)
  fig.text(0.87, 0.76, r'Stellar',   color=myDarkerBlue, size=text_size)

  # Sensors
  fig.text(0.61, 0.79,r'Pyroelectric', color=myDarkGreen,   size=text_size)
  fig.text(0.61, 0.77,r'(Commercial)', color=myDarkGreen,   size=text_size*0.5)
  fig.text(0.63,0.64, r'IR Labs',    color=myDarkGray,    size=text_size)
  fig.text(0.63,0.62, r'(Commercial)',   color=myDarkGray,   size=text_size*0.5)

  fig.text(0.79, 0.51, r'SNSPD',        color=myDarkPink,    size=text_size)
  fig.text(0.54, 0.445,  r'KID',          color=myMediumOrange,size=text_size)
  fig.text(0.43, 0.40,  r'TES',          color=myDarkPurple,  size=text_size)
  fig.text(0.61, 0.39,  r'QCDet',        color=myDarkGray,    size=text_size, rotation=90)

  # Eye guides
  #fig.text(0.89, 0.60,  r'$10^{-15}~\mathrm{W}$', color=myMediumGray, size=0.7*text_size)
  #fig.text(0.89, 0.432, r'$10^{-20}~\mathrm{W}$', color=myMediumGray, size=0.7*text_size)
  #fig.text(0.89, 0.263, r'$10^{-25}~\mathrm{W}$', color=myMediumGray, size=0.7*text_size)

  #fig.text(0.87, 0.50, r'Photon/second', color=myMediumGray, size=0.5*text_size, rotation=15)
  #fig.text(0.87, 0.375, r'Photon/hour',    color=myMediumGray, size=0.5*text_size, rotation=15)

  fig.text(0.36, 0.22, r'\textbf{BREAD} $A_\mathrm{dish} = ' + '{0}'.format(int(Adish)) + '~\mathrm{m}^2$', color=myDarkGray, size=text_size)
  
  if time == 8760:
    integrationT = ', $\Delta t_\mathrm{int} = 1~\mathrm{yr}$'
  elif time == 87600:
    integrationT = ', $\Delta t_\mathrm{int} = 10~\mathrm{yrs}$'
  else:
    integrationT = ', $\Delta t_\mathrm{int}' + ' = {0:.0f}'.format(time) + '~\mathrm{hrs}$'
  #fig.text(0.40, 0.175, r'$\mathrm{SNR}' + ' = {0:.0f}$'.format(snr) + ', $\epsilon_\mathrm{sig}' + ' = {0}$'.format(effic) + integrationT, color=myDarkGray, size=text_size*0.68)
  fig.text(0.36, 0.18, r'$\mathrm{SNR}' + ' = {0:.0f}$'.format(snr) + ', $\epsilon_\mathrm{sig}' + ' = {0}$'.format(effic), color=myDarkGray, size=text_size*0.9)

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
  
  save_name = 'fig_dark_photon_Adish{0}_snr{1}_effic{2}_time{3}'.format(Adish, snr, effic, time)
  print('Saving as {0}.pdf'.format(save_name))
  plt.savefig(save_name + '_photosensors.pdf', format='pdf', dpi=50)
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
def calc_darkPhoton_coupling(nep, mirrorArea, minMass, maxMass, snr=5., effic=0.5, time=1., relicDensity = 0.45):
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

  ratio    = ( snr / 5. )
  noise    = ( nep / 1.0e-21 )
  area     = ( 10.0 / mirrorArea  ) 
  dt       = math.sqrt( 1. / time )
  epsilon  = ( 0.5 / effic )
  rho      = ( 0.3 / relicDensity )

  kineticMixingSq = 3.7 * ratio * noise * area * dt * epsilon * rho
  kineticMixing   = math.sqrt( kineticMixingSq )
  
  # Convert meV to eV, coupling to 1/GeV
  log10minMass = np.log10(minMass * 10**(-3))
  log10maxMass = np.log10(maxMass * 10**(-3))
  log10kineticMixing = np.log10(kineticMixing * 10**(-14)) 
  log10kineticMixing = np.log10(kineticMixing * 10**(-14)) 
 
  return [log10minMass, log10maxMass], [log10kineticMixing, log10kineticMixing]


#__________________________________________
def calc_darkPhoton_coupling_rate(photonsPerHour, mirrorArea, minMass, maxMass):
  '''
  Convert instrument parameters and detected signal power to axion coupling
   - photon rate in photons per hour
   - mirrorArea units is m^2
   - mass in meV
  '''
  relicDensity = 0.45 # GeV/cm^3
  #kineticMixingSq = ( power / float(4.5 * 10**(-21) ) ) * ( 0.3 / relicDensity ) * ( 10.0 / mirrorArea  ) 

  kineticMixingSqMin = (35.6/3600) * ( photonsPerHour ) * ( minMass * 10**(-3) ) * ( 0.3 / relicDensity ) * ( 10.0 / mirrorArea  ) 
  kineticMixingSqMax = (35.6/3600) * ( photonsPerHour ) * ( maxMass * 10**(-3) ) * ( 0.3 / relicDensity ) * ( 10.0 / mirrorArea  ) 
  
  kineticMixingMin = math.sqrt( kineticMixingSqMin )
  kineticMixingMax = math.sqrt( kineticMixingSqMax )
    
  # Convert meV to eV, coupling to 1/GeV
  log10minMass = np.log10(minMass * 10**(-3))
  log10maxMass = np.log10(maxMass * 10**(-3))
  log10kineticMixingMin = np.log10(kineticMixingMin * 10**(-13)) 
  log10kineticMixingMax = np.log10(kineticMixingMax * 10**(-13)) 
  
  return [log10minMass, log10maxMass], [log10kineticMixingMin, log10kineticMixingMax]
  
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