# ASC_and_heat_transport

A collaborative project to investigate the relationship between ASC variability (spatial and temporal) and cross-slope heat transport variability. This will initially use ACCESS-OM2-01 IAF output, potentially evolving to use output from the panan-01, panan-005, and panan-0025 simulations currently underway.

This is a central repository where we can work on analysis scripts together, suggest new analysis directions, discuss figures etc.

## How it works
All aspects of the project are tracked through issues. Create an issue to represent each small task or analysis idea. Issues will develop to include discussion of analysis methods and figures associated with each task.

To start contributing to the code, make your own branch directly in this repository, edit away on your branch, and then submit pull requests between your branch and the master branch (or merge directly if noone else is working on the same code).

## Hackathon schedule (in Canberra/Sydney time zone)

Fortnightly Thursday mornings 9:30-11:30am.

## Data Plan

Only use IAF cycle 3, we have daily temp, salt, uhrho_et, vhrho_nt.

We need:
* Contour data for 3 different isobaths: 650m, 1000m and 1500m. _Done._
* Daily ASC speed on all isobaths, with depth coordinate. Also subset this so we have monthly data. _1000m now computed._
* Daily cross slope heat transport on all isobaths, with depth coordinate. Also subset this so we have monthly data.

## Dimensions to explore

First focus on only East Antarctica, 0-160E.

How do the correlations between ASC and cross-slope heat transport vary with:
* Different isobaths.
* Different depth ranges (e.g. depth integrated vs below 100m vs below 300m). This has the advantage of removing the Ekman layer which is controlled by different dynamics rather than the ASC. Also heat transport in the surface layer is not so relevant for delivering heat to ice shelves.
* Different time frequencies. Interannual (use annual averages), seasonal (use monthly climatology), high frequency (use daily data and remove monthly climatology).
* Different regions. Try: individual gridpoints on the isobath, small single trough-scale regions, large sectors, circumpolar integrals/averages.
* Different measures of the ASC. How sensitive are correlations if we use different aspects of the ASC eg: depth average ASC speed, lower water column ASC speed, a measure of the ASC depth structure.

## Important updates

#### Output locations
Time-mean mask along 1km isobath defining the regimes : `/g/data/v45/txs156/ASC_classification/clustering_mask_time_mean.nc`

Monthly mask along 1km isobath defining the regimes: `/g/data/v45/txs156/ASC_classification/clustering_mask.nc`

Monthly vertically integrated CSHT, across the 1km isobath (IAF, cycle 3): `/g/data/v45/wf4500/ASC_project_files/Cross_slope_heat_transport/OM2_IAF/monthly/vertically_integrated/` 

Daily CSHT, across the 1km isobath , as a function of depth(IAF, cycle 3, all years): `/g/data/v45/wf4500/ASC_project_files/Cross_slope_heat_transport/OM2_IAF/daily_z/`

Daily *Deseasoned* CSHT, across the 1km isobath , as a function of depth(IAF, cycle 3, all years): `/g/data/v45/wf4500/ASC_project_files/Cross_slope_heat_transport/OM2_IAF/CSHT_daily_deseasoned.nc`

Daily ASC speed on XY stepwise grid, over 1km isobath, as a function of depth (IAF, cycle 3): `/g/data/v45/wf4500/ASC_project_files/ASC_speed/OM2_IAF/`

Daily *Deseasoned* ASC speed on XY stepwise grid, over 1km isobath, as a function of depth (IAF, cycle 3), binned into 3 longitude degrees bins: 
`/g/data/v45/wf4500/ASC_project_files/Cross_slope_heat_transport/OM2_IAF/Ualong_daily_deseasoned.nc`

Daily ASC speed on XY stepwise grid, over 1km isobath, as a function of depth (IAF, cycle 3), binned into 3 longitude degrees bins: 
`/g/data/v45/wf4500/ASC_project_files/Binned_ASC_speed/OM2_IAF/`

1500m contour (cleaned up version):
`/g/data/x77/wgh581/Post_Process/access-om2/Antarctic_slope_contour_1500m_no_loops.npz`



**We have decided to do the analysis using CSHT and Mass transport in $\sigma_0$ bins that are indicative of the [water mass structure in each shelf sector](https://github.com/willaguiar/ASC_and_heat_transport/issues/40#issuecomment-2151227185)** Here areteh outputs of CSHT and ASC speed in these new levels:

CSHT in Water mass defined sigma levels:
`/g/data/v45/wf4500/ASC_project_files/Cross_slope_heat_transport/OM2_IAF/WMbins_daily_rho/*`

ASC speed in Water mass defined sigma levels:
`/g/data/v45/wf4500/ASC_project_files/ASC_speed/daily_rho/WMbins_daily_rho/*`
