#!/usr/bin/env python3
'''
Write out text file of limits for axion and dark photon
'''

import matplotlib as mplt
mplt.use('Agg') # So we can use without X forwarding

import numpy as np
import pandas as pd
import os, json, math, csv, argparse, datetime

from common import *

#__________________________________________
def main():

  # --------------------------------------
  # Axions
  # --------------------------------------
  # Default values
  do_axion = True
  reldens  = 0.45 # relic density GeV/cm^3
  Adish    = 10. # dish area in m^2
  Bfield   = 10. # Tesla
  snr      = 5.  # signal to noise
  effic    = 0.5 # signal detection efficiency
  time     = 240.  # integration time in hours

  l_x, l_y = [], []

  # TES [0.2, 1.2] meV
  x, y = calc_coupling_nep(2e-21, Adish, Bfield, 0.2, 2, snr, effic, time*100, reldens, do_axion)
  for xi in x: l_x.append(xi)
  for yi in y: l_y.append(yi)
  # QCD [2, 125] meV
  x, y = calc_coupling_dcr(4e-4, Adish, Bfield, 2, 124, snr, effic, time*100, reldens, do_axion)
  for xi in x: l_x.append(xi)
  for yi in y: l_y.append(yi)
  # SNSPD [207, 830] meV
  x, y = calc_coupling_dcr(1e-8, Adish, Bfield, 124, 830, snr, effic, time*100, reldens, do_axion)
  for xi in x: l_x.append(xi)
  for yi in y: l_y.append(yi)
  print('---------------------------------\nAxion\n--------------------------------')
  print('mass coupling')
  file_name = 'limits/common/Ax_BREAD.txt'
  with open(file_name, 'w') as fout:
    for xi, yi in zip(l_x, l_y):
      line = '{0:.4g} {1:.4g}'.format(xi, yi)
      fout.write(line + '\n')
      print(line)
  print('\nWrote summary curve to file: {0}'.format(file_name))

  # --------------------------------------
  # Dark photons
  # --------------------------------------
  l_x, l_y = [], []
  # Default values
  do_axion = False
  reldens  = 0.45 # relic density GeV/cm^3
  Adish    = 10. # dish area in m^2
  snr      = 5.  # signal to noise
  effic    = 0.5 # signal detection efficiency
  time     = 240.  # integration time in hours

  # TES/KID [0.2, 125] meV
  x, y = calc_coupling_nep(2e-21, Adish, 1., 0.2, 2, snr, effic, time*100, reldens, do_axion)
  for xi in x: l_x.append(xi)
  for yi in y: l_y.append(yi)
  # QCD [2, 125] meV
  x, y = calc_coupling_dcr(4e-4, Adish, 1., 2, 124, snr, effic, time*100, reldens, do_axion)
  for xi in x: l_x.append(xi)
  for yi in y: l_y.append(yi)
  # SNSPD [207, 830] meV
  x, y = calc_coupling_dcr(1e-8, Adish, 1., 124, 830, snr, effic, time*100, reldens, do_axion)
  for xi in x: l_x.append(xi)
  for yi in y: l_y.append(yi)
  print('\n---------------------------------\nDark photon\n---------------------------------')
  print('mass coupling')
  with open('limits/common/DP_BREAD.txt', 'w') as fout:
    for xi, yi in zip(l_x, l_y):
      line = '{0:.4g} {1:.4g}'.format(xi, yi)
      fout.write(line + '\n')
      print(line)
  
  print('\nWrote summary curve to file: {0}'.format(file_name))

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