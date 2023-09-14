# ASC_and_heat_transport

A collaborative project to investigate the relationship between ASC variability (spatial and temporal) and cross-slope heat transport variability. This will initially use ACCESS-OM2-01 IAF output, potentially evolving to use output from the panan-01, panan-005, and panan-0025 simulations currently underway.

This is a central repository where we can work on analysis scripts together, suggest new analysis directions, discuss figures etc.

## How it works
All aspects of the project are tracked through issues. Create an issue to represent each small task or analysis idea. Issues will develop to include discussion of analysis methods and figures associated with each task.

To start contributing to the code, make your own branch directly in this repository, edit away on your branch, and then submit pull requests between your branch and the master branch (or merge directly if noone else is working on the same code).

## Hackathon schedule (in Canberra/Sydney time zone)

Fortnightly on TBD:

*June 8th, 9:30 am to 11:30 am

*June 22nd,  9:30 am to 11:30 am


## Task List

<del>Task 1 (Paul,Taimoor)- Find alternative ways to define the ASC (e.g., cluster analysis)
- We decided to use different isobaths to calculate ASC and heat transport on. Our attempts to find a different definition weren't promising and we didn't want to spend too much time on this part.

Task 2 (Wilton)- Local (point by point) correlations between the ASC speed and cross slope heat transport. [These correlations for the monthly and vertically averaged variables vary greatly along the contour](https://github.com/willaguiar/ASC_and_heat_transport/issues/1). The exception was a region in East Antarctica. 
- Next: Do the same correlations using [a] Daily data, [b] vars integrated solely below the mixed layer to remove the effect of the surface poleward heat transport by the Ekman transport during the summer.

Task 3-  Subtract the daily climatology of ASC from the daily timeseries to look at variability at frequencies higher than seasonal.

Task 4 (Fabio)- Interannual variability ( annual averages of speed, Heat transport in IAF simulation). Plot timeseries ( different regions) and correlations.

Task 5 (Wilton)-  Meridional heat convergence in the East Antarctic Peninsula (and Western Ross sea). [Check discussions on this task here](https://github.com/willaguiar/ASC_and_heat_transport/issues/10). Rerunning the cross-slope heat transport calculations using this array would use extra SUs, that we could save by analysing  the Antarctic Peninsula as a single sector and discounting the Meridional heat convergence in the sector. So we should stick to the contour with purely zonal convergence.


Task 6 (Wilton,Fabio)- Alternative way of calculating the cross-shelf heat transport. Specifically, transpose the zonal and meridional components to be along slope and cross slope. Then calculate the heat convergence along-slope.

Task 7 (Paul,Ellie,Fabio)- Literature review (Dynamics, and previous work on ASC vs cross slope heat transport).



## Important updates
#### Grid
We previously discussed wether we should calculate the cross-slope heat transports and ASC speeds on the U-grid points (Cell corners), or on the XY stepwise grid (Cell faces). Both would require some interpolation, however, interpolating ASC speed is more straightforward, so we decided to stick to the XY stepwise grid (Cell faces) for all calculations. The python codes to calculate ASC speed on the proper grid are `/Python/*_along_slope_velocities_XYgrid.py`

#### Output locations
Daily vertically integrated CSHT, on the 1000m isobath (RYF, yr 2170 to 2179): `/home/156/wf4500/x77_wf4500/ASC_project/cross_slope_transport/Daily/01deg_jra55v13_ryf9091`

Monthly vertically integrated CSHT, on the 1000m isobath (IAF, cycle 3): `/home/156/wf4500/x77_wf4500/ASC_project/cross_slope_transport/Monthly/01deg_jra55v140_iaf_cycle3`


