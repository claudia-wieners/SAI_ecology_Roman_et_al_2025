# SAI_ecology_Roman_et_al_2025

This repository contains the codes needed to reproduce the figures in Roman de Miguel, Pflüger, Wieners, Gurevitch and Harrison, 2025: Beyond climate velocity: Indicators for ecological impacts of late, sudden SRM (accepted at OOCC). 

General remarks:
- The python notebooks (.ipynb) require the auxiliary python (.py) files ro run, these are assumed to be in the same folder as the main script. 
- Data for 2m-temperature (TREFHT) and precipitation (convective, PRECC, plus large-scale, PRECL), as well as some landmasks, are in the folder "data", assumed to be a subfolder of the folder where the notebooks are run. These data are not saved here due to file size, but can be found under: https://doi.org/10.5281/zenodo.17940938 
- The intermediate data for the matlab plotscripts (.m) are now in a folder called "intermediate data", but in the scripts it is assumed that these data are directly in the same folder as the matlab scripts.

Fig.1: Time series of global mean surface temperature
- Script: GMST_scenarios.ipynb
- Input data: data / <scenario> / temperature data files; here, <scenario> is Control, SAI 2020, or SAI 2080J. 

Fig.2: not applicable (infographics without data)

Fig.3: Climate velocities
We first computed the climate velocities, saved the values as intermediate data (.nc files), then ran a sceond script (matlab) to plot them. 
- Script for computing velocities: Climatevelocities_computation.ipynb
- Input data: data / <scenario> / temperature and precip data files; here, <scenario> is Control, SAI 2020, or SAI 2080. 
- Intermediate climate velocity data: v_T_xxxx.nc where xxxx is con (control, i.e. no SAI), sai2020 (SAI from 2020) and sai2080 (SAI from 2080).
- Auxiliary data for plotting: coastlines_0_360.mat (for plotting coastlines)
- Plotting script: climvel_plotter_projectedmap.m 

Fig. 4: Isotherm shifts
Script: Isolines_computation.ipynb
Input data: data / <scenario> / temperature and precip data files; here, <scenario> is Control, SAI 2020, or SAI 2080. 

Fig. 5a,b and S1: Maps of simplified Köppen climate zones 
We first computed the simplified Köppen climate classifications, saved them as intermediate data (.mat file), then ran a sceond script to plot them.
- Script for computing Köppen: Koppen_class_simple.ipynb. NOTE: need to restart the kernel and delete all cell output before making a new computation for another scenario or period. 
- Input data: data / <scenario> / landmask, temperature and precip data files; here, <scenario> is Control, SAI 2020, or SAI 2080.  
- Intermediate Köppen classification data: Koppendata_xxx.mat 
  (7 files; where xxx is the scenario and time. CONT=control (no SAI); S20 = SAI 2020; S80 = SAI 2080; 
  e stands for early (2050-2065), m for medium (2065-2080), l for late (2085-2100); REF is control 2020-2035.
- Auxiliary data for plotting: coastlines_0_360.mat (for plotting coastlines)
- Plotting scripts: Koppen_plotter_projectedmap.m  
  and Koppen_formatter.m (auxiliary function)

Fig. 5c-f and fig. S2 Maps of changes of the Köppen climate zones and matrices quantifying area changes
as fig. 5a, except that we use a different final script. 
- Plotting script: Koppen_plotter_changed_area_projectedmap.m
  and Koppen_formatter.m (auxiliary function)


Auxiliary files Python environment: 
- Roman_Ecology_specfile.txt (platform-specific)
- Roman_Ecology.yml 

