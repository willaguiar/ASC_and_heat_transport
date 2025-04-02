# Antarctic Slope Current (ASC) control on cross-slope ehat transport

A collaborative project to investigate the relationship between ASC variability (spatial and temporal) and cross-slope heat transport variability. This will use ACCESS-OM2-01 IAF output.

This is a central repository where we can work on analysis scripts together, suggest new analysis directions, discuss figures etc.

## How it works
All aspects of the project are tracked through issues. Create an issue to represent each small task or analysis idea. Issues will develop to include discussion of analysis methods and figures associated with each task.

To start contributing to the code, make your own branch directly in this repository, edit away on your branch, and then submit pull requests between your branch and the master branch (or merge directly if noone else is working on the same code).

## Hackathon schedule (in Canberra/Sydney time zone)

Fortnightly Thursday mornings 9:30-11:30am.

## Write-a-ton (in Canberra/Sydney time zone)
To finish up the project we will run a in person Write-a-ton in Canberra. This will be happening April 3rd/2025, and April 4th/2025. The goal of the Write-a-ton is to have a final version of a manuscriupt, ready for submission by the last day.

## Data Plan
Use IAF cycle 3 for monthly,  and interannual analysis (we have daily temp, salt, uhrho_et, vhrho_nt). Seasonal analysis will be done using daily data.


We need:
* Contour data for 3 different isobaths: 650m, 1000m and 1500m. _Done. We chose to use the 1000m isobath only_
* Daily ASC speed on 1000m isobath in IAF experiment, in $\sigma_0$ coordinate, resampled monthly. _1000m now computed._
* Daily cross slope heat transport on 1000m isobath in IAF experiment, in $\sigma_0$ coordinate, subset to monthly data _Done_

All data used is binned into 10 degrees longitude bins, which smoothes out the high spatial variability of the CSHT

## Dimensions to explore

Using the 10 degree binned data in $\sigma_0$ we will test...

How do the correlations between ASC and cross-slope heat transport vary with:
* Different water masses (e.g. AASW, CDW and DSW).
* Different time frequencies. Interannual (use annual averages), seasonal (use monthly climatology), high frequency (use daily data and remove monthly climatology).
* Different regions. Try: individual gridpoints on the isobath, small single trough-scale regions, large sectors, circumpolar integrals/averages.


## Important updates
Up until now, we have decided to analyze the relationship between the CSHT and ASC speed in the Water masss space. We are defining three water masses sitting on the 1000m isobath based on their age, normalized for each longitude bin:
* _CDW_: Defined as waters in each logintude bin with normalized age above 0.95.
* _AASW_: Surface waters above the CDW
* _DSW_: Dense waters below CDW

#### Output locations

Daily, 10deg binned,  CSHT in $\sigma_0$ for the IAF simualtion (PS: _without Zonal Convergence_ you have to do monthly averages before using the data):
`/g/data/hh5/tmp/access-om/fbd581/ASC_project/Lon_binning/CSHT_potrho0_10degBins.nc` 

Daily, 10deg binned, CSHT+ZC in $\sigma_0$ for the IAF simualtion (PS: you have to do monthly averages before using the data):
`/g/data/hh5/tmp/access-om/fbd581/ASC_project/Lon_binning/CSHT_ZC_potrho0_10degBins.nc` 

Monthly Water masses masks for each 10 degree bins for the IAF simualtion:
`/g/data/v45/wf4500/ASC_project_files/**_age_mask.nc`

or
`/g/data/ik11/users/wf4500/ASC_project_files/**_age_mask.nc`

Daily _climatology_ of Water masses masks for each 10 degree bins for the _daily climatology_ of the IAF simualtion:
`/g/data/v45/wf4500/ASC_project_files/****_age_Climatology_mask.nc`

or
`/g/data/ik11/users/wf4500/ASC_project_files/****_age_Climatology_mask.nc`

Monthly mean ASC speed  in $\sigma_0$ for the IAF simualtion:
`/g/data/v45/wf4500/ASC_project_files/ASC_10deg_bin.nc`

or
`/g/data/ik11/users/wf4500/ASC_project_files/ASC_10deg_bin.nc`

Daily climatology of ASC speed in $\sigma_0$ for the IAF simualtion:
`/g/data/v45/wf4500/ASC_project_files/DailyU_climatology_rho_10deg_bins.nc`

or
`/g/data/ik11/users/wf4500/ASC_project_files/DailyU_climatology_rho_10deg_bins.nc`
