# DM sensitivity
Basic plotting scripts for projecting BREAD sensitivity in axion (spin 0) and dark photon (spin 1) scenarios.
Packages used are numpy (1.19.4), pandas (1.1.4), matplotlib (3.3.3), and a tex (TeX Live 2016) backend. 
These scripts create one plot per run. 

# Scripts for plots in paper

The following scripts are used to make the various plots in the paper [arXiv:2111.12103](https://arxiv.org/abs/2111.12103):
* Fig 2: the `dark_photon.py` and `axion.py` scripts make the main sensitivity plots for the paper. Each line for BREAD has to be manually coded in the scropt. The plot is created by running the script. 
* Fig 3: the single plot for QCD axion (KSVZ, DFSZ) emitted photon rate is made using `photon_rate.py`
* Fig 4: the sensitivity to the axion-photon coupling (normalized to DFSZ) vs experimental parameters are created by the
 `coupling_vs_ABsq.py`, `coupling_vs_efficiency.py`, `coupling_vs_nep.py` scripts.
* Fig 5: the `dark_photon-general.py` and `axion-general.py` make the generic bolometer/photocounter plots

* Table 2: The `calc_coupling_sensitivity.py` is a simple script to compute the coupling sensitivity (eq. 2 in paper) normalized to KSVZ and DFSZ for axions, and 1e-14 for dark photons displayed in this table. 

* The `print_summary_curve.py` script writes out the numerical values of summary limit curves to a space-separated text file. 

## References:
* BREAD concept paper: Liu et al https://arxiv.org/abs/2111.12103
* Dish antenna paper: Horns et al https://arxiv.org/abs/1212.2970
* Astrophysical constraints: Arias et al https://arxiv.org/abs/1201.5902
* QCD axion band: Di Luzio et al https://arxiv.org/abs/1610.07593
* Dark photon limits cookbook: Caputo et al https://arxiv.org/abs/2105.04565 
* Axion and dark photon limits: https://github.com/cajohare/AxionLimits/ 