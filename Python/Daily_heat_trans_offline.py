# # Compute monthly heat transport along contour offline (when is not outputed by the model). DOes not discounts the freezing temperature

import cosima_cookbook as cc
from cosima_cookbook import distributed as ccd
import matplotlib.pyplot as plt
import numpy as np
import netCDF4 as nc
import xarray as xr
import glob,os
import cmocean.cm as cmocean
import os
import sys
import glob

import logging
logging.captureWarnings(True)
logging.getLogger('py.warnings').setLevel(logging.ERROR)
logging.getLogger('distributed.utils_perf').setLevel(logging.ERROR)

from dask.distributed import Client

import climtas.nci
import warnings # ignore these warnings
warnings.filterwarnings("ignore", category = FutureWarning)
warnings.filterwarnings("ignore", category = UserWarning)
warnings.filterwarnings("ignore", category = RuntimeWarning)

if __name__ == '__main__':

	climtas.nci.GadiClient()
	
	session = cc.database.create_session()
	
	#imon=int(sys.argv[1])
	#if imon<=9:
	#    month='0' + str(imon)
	#else:
	#    month=str(imon)
	
	
	monthdays = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
	month = str(int(sys.argv[1]))
	month = month.zfill(2)
	year = str(sys.argv[2])
	start_time=year+'-'+month  
	 #Start_time0 and end_time0 are for importing the daily transport, and it ahsthe number of days in the month    
	start_time0=year+'-'+month +'-01'     
	end_time0=year+'-'+month +'-' + str(monthdays[np.int(np.int(sys.argv[1])-1)])
	print(start_time0) 
	print(end_time0) 
	exp = '01deg_jra55v13_ryf9091'
	
	print("Start date =" + start_time) 
	year2=str(np.int(start_time[0:4])+1)
	month2=str(np.int(start_time[5:7])+1)
	month2=str(int(month2))
	month2 = month2.zfill(2)    
	print("month2 is =" + month2) 
	print("year2 is =" + str(year2))     
	
	imon = np.int(sys.argv[1])
	if imon <12:
	    end_time=year+'-'+month2
	else:
	    end_time=year2+'-01'
	    
	print("End date =" + end_time) 
	
	      
	time_period = str(int(start_time[:4]))+'-'+str(int(end_time[:4]))
	
	# reference density value:
	rho_0 = 1035.0
	# specific heat capacity of sea water:
	cp = 3992.1
	lat_range = slice(-90,-59)
	
	isobath_depth = 1000
	
	# pick a freezing point temp:
	temp_freezing = -3.82
	
	
	
	ndays= str(monthdays[np.int(np.int(sys.argv[1])-1)])
	start_time=year+'-'+month+'-01'
	end_time=year+'-'+month+'-'+ ndays
	
	start_time,end_time
	
	
	
	yt_ocean = cc.querying.getvar(exp,'yt_ocean',session,n=1)
	yt_ocean = yt_ocean.sel(yt_ocean=lat_range)
	yu_ocean = cc.querying.getvar(exp,'yu_ocean',session,n=1)
	yu_ocean = yu_ocean.sel(yu_ocean=lat_range)
	xt_ocean = cc.querying.getvar(exp,'xt_ocean',session,n=1)
	xu_ocean = cc.querying.getvar(exp,'xu_ocean',session,n=1)
	
	
	# interpolate temperature to vhrho grids:
	
	temp = cc.querying.getvar(exp,'temp',session,ncfile='%daily%',start_time=start_time,
	                           end_time=end_time).sel(time=slice(start_time,end_time))
	temp = temp.sel(yt_ocean=lat_range)
	
	temp_s = temp.copy()
	temp_s=temp_s.rename(yt_ocean='y_ocean',xt_ocean='x_ocean') #added by wilton
	temp_s.coords['y_ocean'] = np.array(yu_ocean)
	
	temp_n = temp.roll(yt_ocean=-1)
	temp_n=temp_n.rename(yt_ocean='y_ocean',xt_ocean='x_ocean') #added by wilton
	temp_n.coords['y_ocean'] = np.array(yu_ocean)
	# # temp_ygrid will be on the vhrho grid:
	temp_ygrid = (temp_s + temp_n)/2
	
	
	
	vhrho = cc.querying.getvar(exp,'vhrho_nt',session,ncfile='%daily%',start_time=start_time,end_time=end_time).sel(time=slice(start_time,end_time))
	
	# select latitude range and this month:
	vhrho = vhrho.sel(yt_ocean=lat_range)
	
	# Note that vhrho is defined as the transport across the northern edge of a tracer cell so its coordinates 
	#       should be (yu_ocean, xt_ocean).
	#  uhrho is defined as the transport across the eastern edge of a tracer cell so its coordinates should 
	#       be (yt_ocean, xu_ocean).
	#  However we will keep the actual name as simply y_ocean/x_ocean irrespective of the variable
	#       to make concatenation and sorting possible.
	vhrho = vhrho.rename({'yt_ocean':'y_ocean', 'xt_ocean':'x_ocean'})
	
	# vhrho and uhrho are saved with the wrong coords, change them:
	vhrho.coords['y_ocean'] = np.array(yu_ocean)
	
	
	
	# compute heat transport and depth sum:
	temp_yflux = (temp_ygrid*vhrho).sum('st_ocean')
	
	
	### Convert to transports 
	
	# First we also need to change coords on dxu, dyt, so we can multiply the transports:
	dxu = cc.querying.getvar(exp, 'dxu',session, n=1, ncfile = 'ocean_grid.nc')
	# select latitude range:
	dxu = dxu.sel(yu_ocean=lat_range)
	# remove geolon_t/geolat_t coordinates:
	dxu = dxu.reset_coords().dxu 
	dxu.coords['xu_ocean'] = xt_ocean.values
	dxu = dxu.rename({'yu_ocean':'y_ocean', 'xu_ocean':'x_ocean'}) 
	
	# convert to transport:
	temp_yflux = temp_yflux*dxu/rho_0
	
	
	# save temp_yflux for this year:
	ds = xr.Dataset({'temp_yflux': temp_yflux})
	#outpath = '/g/data/v45/akm157/model_data/access-om2/'+exp+'/daily_heat_transport/temp_yflux_from_daily_'+year+'.nc'
	outpath = '/g/data/x77/wf4500/ASC_project/model_data/access-om2/'+exp+'/daily_heat_transport/temp_yflux_from_daily_'+year+'_'+month+'.nc'
	
	
	ds.to_netcdf(outpath)
	
	
	
	
	# interpolate temperature to uhrho grids:
	
	temp = cc.querying.getvar(exp,'temp',session,ncfile='%daily%',start_time=start_time,end_time=end_time).sel(time=slice(start_time,end_time))
	temp = temp.sel(yt_ocean=lat_range)
	
	# interpolate to correct grid:
	temp = temp.rename({'yt_ocean':'y_ocean', 'xt_ocean':'x_ocean'}) 
	temp_w = temp.copy()
	temp_w.coords['x_ocean'] = xu_ocean.values
	temp_e = temp.roll(x_ocean=-1)
	temp_e.coords['x_ocean'] = xu_ocean.values
	# temp_xgrid will be on the uhrho grid:
	temp_xgrid = (temp_e + temp_w)/2
	
	
	uhrho = cc.querying.getvar(exp,'uhrho_et',session,ncfile='%daily%',start_time=start_time,end_time=end_time).sel(time=slice(start_time,end_time))
	
	# select latitude range and this month:
	uhrho = uhrho.sel(yt_ocean=lat_range)
	
	# Note that vhrho is defined as the transport across the northern edge of a tracer cell so its coordinates 
	#       should be (yu_ocean, xt_ocean).
	#  uhrho is defined as the transport across the eastern edge of a tracer cell so its coordinates should 
	#       be (yt_ocean, xu_ocean).
	#  However we will keep the actual name as simply y_ocean/x_ocean irrespective of the variable
	#       to make concatenation and sorting possible.
	uhrho = uhrho.rename({'yt_ocean':'y_ocean', 'xt_ocean':'x_ocean'})
	
	# vhrho and uhrho are saved with the wrong coords, change them:
	uhrho.coords['x_ocean'] = xu_ocean.values
	
	
	
	
	# compute heat transport and depth sum:
	temp_xflux = (temp_xgrid*uhrho).sum('st_ocean')
	
	
	### Convert to transports 
	
	# First we also need to change coords on dxu, dyt, so we can multiply the transports:
	dyt = cc.querying.getvar(exp, 'dyt',session, n=1, ncfile = 'ocean_grid.nc')
	# select latitude range:
	dyt = dyt.sel(yt_ocean=lat_range)
	# remove geolon_t/geolat_t coordinates:
	dyt = dyt.reset_coords().dyt
	dyt.coords['xt_ocean'] = xu_ocean.values
	dyt = dyt.rename({'yt_ocean':'y_ocean','xt_ocean':'x_ocean'})
	
	# convert to transports:
	temp_xflux = temp_xflux*dyt/rho_0
	
	
	
	# save temp_yflux for this year:
	ds = xr.Dataset({'temp_xflux': temp_xflux})
	outpath = '/g/data/x77/wf4500/ASC_project/model_data/access-om2/'+exp+'/daily_heat_transport/temp_xflux_from_daily_'+year+'_'+month+'.nc'
	ds.to_netcdf(outpath)
	print('temp_yflux and temp_xflux saved')
	
	
	
	quit()
