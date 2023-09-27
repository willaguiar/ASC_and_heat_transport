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
* Contour data for 3 different isobaths: 700m, 1000m and 1500m.
* Daily ASC speed on all isobaths, with depth coordinate. Also subset this so we have monthly data.
* Daily cross slope heat transport on all isobaths, with depth coordinate. Also subset this so we have monthly data.

## Dimensions to explore

How do the correlations between ASC and cross-slope heat transport vary with:
* Different isobaths.
* Different depth ranges (e.g. depth integrated vs below 100m vs below 300m). This has the advantage of removing the Ekman layer which is controlled by different dynamics rather than the ASC. Also heat transport in the surface layer is not so relevant for delivering heat to ice shelves.
* Different time frequencies. Interannual (use annual averages), seasonal (use monthly climatology), high frequency (use daily data and remove monthly climatology).
* Different regions. Try: individual gridpoints on the isobath, small single trough-scale regions, large sectors, circumpolar integrals/averages.
* Different measures of the ASC. How sensitive are correlations if we use different aspects of the ASC eg: depth average ASC speed, lower water column ASC speed, a measure of the ASC depth structure.

## Important updates

#### Output locations
Monthly vertically integrated CSHT, on the 1000m isobath (IAF, cycle 3): `/home/156/wf4500/x77_wf4500/ASC_project/cross_slope_transport/Monthly/01deg_jra55v140_iaf_cycle3`

Daily ASC speed on XY stepwise grid, as a function of depth (IAF, cycle 3): `/g/data/x77/wf4500/ASC_project/ASC_speed/daily_z/OM2_IAF_XYgrid/`

Daily ASC speed on XY stepwise grid, as a function of depth (IAF, cycle 3), binned into 3 longitude degrees bins: 
/g/data/nm03/wf4500/ASC_project/model_data/access-om2/binned_ASC_speed/daily/IAF/

