#!/bin/bash
echo "Setting up environment from cvmfs..."

echo "Setting up LaTeX (texlive 2016)"
export PATH=/cvmfs/sft.cern.ch/lcg/external/texlive/2016/bin/x86_64-linux:$PATH

echo "Setting up LCG 100 (x86_64-centos7-gcc8-opt)"
source /cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-gcc8-opt/setup.sh

