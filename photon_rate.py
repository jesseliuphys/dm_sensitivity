#!/usr/bin/env python
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
  x_txt = r'$\mathrm{Axion~mass~[eV}]$'
  y_txt = r'$\mathrm{Emitted~photons~[day}^{-1}]$'
    
  # Axis label properties
  plt.xlabel(x_txt, labelpad=20, size=35) 
  plt.ylabel(y_txt, labelpad=20, size=35)

  ax.set_xscale('log')
  ax.set_yscale('log')
  # Draw the upper THz axis
  ax2 = ax.twiny()
  ax2.set_xlim(2.42e-5, 242)
  ax2.set_xlabel(r'$\mathrm{Frequency}~[\mathrm{THz}]$', labelpad=20, size=35)
   
  ax2.set_xscale('log')

  # ------------------------------------------------------- 
  # Add plot text
  # ------------------------------------------------------- 


  #fig.text(0.89, 0.263, r'$10^{-25}~\mathrm{W}$', color=myMediumGray, size=0.7*text_size)

  #fig.text(0.87, 0.50, r'Photon/second', color=myMediumGray, size=0.5*text_size, rotation=15)
  #fig.text(0.87, 0.375, r'Photon/hour',    color=myMediumGray, size=0.5*text_size, rotation=15)

  #fig.text(0.77, 0.18, r'$A_\mathrm{antenna} = 10~\mathrm{m}^2$', color=myDarkGray, size=text_size)
  fig.text(0.21, 0.30, r'\textbf{BREAD}', color=myDarkGray, size=text_size)
  fig.text(0.21, 0.24, r'$A_\mathrm{dish} = 10~\mathrm{m}^2$', color=myDarkGray, size=text_size)
  fig.text(0.21, 0.18, r'$B_\mathrm{ext} = 10~\mathrm{T}$', color=myDarkGray, size=text_size)

  # Adjust axis ticks
  ax.minorticks_on()
  #ax2.minorticks_on()

  ax.tick_params('x', length=12, width=1, which='major', labelsize='28', pad=10)
  ax.tick_params('x', length=6,  width=1, which='minor') 
  ax2.tick_params('x', length=12, width=1, which='major', labelsize='28', pad=10)
  ax2.tick_params('x', length=6,  width=1, which='minor') 
  ax.tick_params('y', length=12, width=1, which='major', labelsize='28', pad=10)
  ax.tick_params('y', length=6,  width=1, which='minor') 
   
  plt.tight_layout(pad=0.3)
  plt.subplots_adjust( top=0.85,left=0.17 )
  
  save_name = 'photon_rate_qcdaxion'
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
def calc_axion_rate(mirrorArea, Bfield, rhoDM, minMass, maxMass, DFSZ=False):

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