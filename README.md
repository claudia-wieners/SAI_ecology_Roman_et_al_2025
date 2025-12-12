# SAI_ecology_Roman_et_al_2025

This repository contains the codes needed to reproduce the figures in Roman de Miguel, Pflüger, Wieners, Gurevitch and Harrison, 2025: Beyond climate velocity: Indicators for ecological impacts of late, sudden SRM (accepted at OOCC). 

Fig.1: Time series of global mean surface temperature
Script: 
Input data: 

Fig.2: not applicable (infographics without data)

Fig.3: Climate velocities
We first computed the climate velocities, saved the values as intermediate data (.nc files), then ran a sceond script (matlab) to plot them. 
- Script for computing velocities:
- Input data: 
- Intermediate climate velocity data: v_T_xxxx.nc where xxxx is con (control, i.e. no SAI), sai2020 (SAI from 2020) and sai2080 (SAI from 2080).
- Auxiliary data for plotting: coastlines_0_360.mat (for plotting coastlines)
- Plotting script: climvel_plotter_projectedmap.m 

Fig. 4: Isotherm shifts
Script: 
Input data: 

Fig. 5a,b and S1: Maps of simplified Köppen climate zones 
We first computed the simplified Köppen climate classifications, saved them as intermediate data (.mat file), then ran a sceond script to plot them.
- Script for computing Köppen:
- Input data: 
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


 
Fig. S2: Changes in area covered by Köppen zones
Input: the intermediate data (Köppen classification) explained unter fig. 5
Script: 


Auxiliary files: 
yml file of the python environment
