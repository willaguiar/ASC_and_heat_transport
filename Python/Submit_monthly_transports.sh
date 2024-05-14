#!/bin/bash
#PBS -N HTday_sigma
#PBS -P x77
#PBS -q hugemembw
#PBS -l walltime=4:30:00
#PBS -l mem=510GB
#PBS -l software=netcdf
#PBS -l ncpus=14
#PBS -l storage=gdata/v45+gdata/hh5+gdata/cj50+gdata/ik11+scratch/x77+scratch/v45+gdata/x77+gdata/nm03+scratch/nm03
#PBS -v month,year

#Calculates the monthly Along and cross isobath (1000m) speed monthly, for the specified year
## Note, run this with:
## qsub -v month=2,year=1972 Submit_monthly_transports.sh
## For running the along and cross slope transports, you will need more memory than 256 GB, Id say at least 1000GB

script_dir=/home/156/wf4500/x77_wf4500/ASC_project/GH/ASC_and_heat_transport/Python/
cd $script_dir

# load conda
module use /g/data/hh5/public/modules
module unload conda
#module load conda/analysis3-23.04
module load conda/analysis3-unstable
module list

# run
# python3 Monthly_cross_slope_heat_transport.py $month $year &>> ASC_job${month}.out


python3 Daily_cross_slope_heat_transport_IAF_from_offline_sigmaBinned.py $month $year &>> CSHTiaf_${year}_${month}.out
#python3 Binning_ASC_speed.py $month $year &>> ASCbinning_${year}_${month}.out
#python3 Daily_heat_trans_offline-Z.py $month $year &>> tempXY_${year}_${month}.out



