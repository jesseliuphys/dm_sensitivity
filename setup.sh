#!/bin/bash
echo "Setting up LaTeX"
export PATH=/cvmfs/sft.cern.ch/lcg/external/texlive/2016/bin/x86_64-linux:$PATH

echo "Setting up pyanalysis 2.0 (matplotlib etc) from LCG 87"
source /cvmfs/sft.cern.ch/lcg/releases/LCG_87/pyanalysis/2.0/x86_64-slc6-gcc49-opt/pyanalysis-env.sh

