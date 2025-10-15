# # Compute monthly cross-slope heat transport on IAF simulations, using monthly output. Already discounted freezing temperature.

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
	end_time0=year+'-'+month +'-' + str(monthdays[int(int(sys.argv[1])-1)])
	print(start_time0) 
	print(end_time0) 
	exp = '01deg_jra55v140_iaf_cycle3'
	
	print("Start date =" + start_time) 
	year2=str(int(start_time[0:4])+1)
	month2=str(int(start_time[5:7])+1)
	month2=str(int(month2))
	month2 = month2.zfill(2)    
	print("month2 is =" + month2) 
	print("year2 is =" + str(year2))     
	
	imon = int(sys.argv[1])
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
	
	
	print("importing isobath mask")
	outfile = '/g/data/v45/akm157/model_data/access-om2/Antarctic_slope_contour_'+str(isobath_depth)+'m.npz'
	data = np.load(outfile)
	mask_y_transport = data['mask_y_transport']
	mask_x_transport = data['mask_x_transport']
	mask_y_transport_numbered = data['mask_y_transport_numbered']
	mask_x_transport_numbered = data['mask_x_transport_numbered']
	
	yt_ocean = cc.querying.getvar(exp,'yt_ocean',session,n=1)
	yt_ocean = yt_ocean.sel(yt_ocean=lat_range)
	yu_ocean = cc.querying.getvar(exp,'yu_ocean',session,n=1)
	yu_ocean = yu_ocean.sel(yu_ocean=lat_range)
	xt_ocean = cc.querying.getvar(exp,'xt_ocean',session,n=1)
	xu_ocean = cc.querying.getvar(exp,'xu_ocean',session,n=1)   
	
	mask_x_transport =xr.DataArray(data['mask_x_transport']).assign_coords({"dim_0": np.array(yt_ocean),"dim_1": np.array(xu_ocean)}).rename(dim_0="yt_ocean",dim_1="xu_ocean")
	mask_y_transport =xr.DataArray(data['mask_y_transport']).assign_coords({"dim_0": np.array(yt_ocean),"dim_1": np.array(xu_ocean)}).rename(dim_0="yu_ocean",dim_1="xt_ocean")
	mask_x_transport_numbered =xr.DataArray(data['mask_x_transport_numbered']).assign_coords({"dim_0": np.array(yt_ocean),"dim_1": np.array(xt_ocean)}).rename(dim_0="yt_ocean",dim_1="xt_ocean")
	mask_y_transport_numbered =xr.DataArray(data['mask_y_transport_numbered']).assign_coords({"dim_0": np.array(yt_ocean),"dim_1": np.array(xt_ocean)}).rename(dim_0="yt_ocean",dim_1="xt_ocean")
	
	num_points = int(np.maximum(np.max(mask_y_transport_numbered),np.max(mask_x_transport_numbered)))
	    
	lat_along_contour = np.zeros((num_points))
	lon_along_contour = np.zeros((num_points))
	
	# locations for zonal transport:
	x_indices_masked = mask_x_transport_numbered.stack().values
	x_indices = np.sort(x_indices_masked[x_indices_masked>0])
	for count in x_indices:
	    count = int(count)
	    jj = int(np.where(mask_x_transport_numbered==count)[0])
	    ii = int(np.where(mask_x_transport_numbered==count)[1])   
	    lon_along_contour[count-1] = xu_ocean[ii].values
	    lat_along_contour[count-1] = mask_x_transport_numbered.yt_ocean[jj].values
	    
	# locations for meridional transport:
	y_indices_masked = mask_y_transport_numbered.stack().values
	y_indices = np.sort(y_indices_masked[y_indices_masked>0])
	for count in y_indices:
	    count = int(count)
	    jj = np.where(mask_y_transport_numbered==count)[0]
	    ii = np.where(mask_y_transport_numbered==count)[1]
	    lon_along_contour[count-1] = mask_x_transport_numbered.xt_ocean[ii].values
	    lat_along_contour[count-1] = yu_ocean[jj].values
	
	# Note temp_yflux_adv is also positioned on north centre edge of t-cell.
	
	print("importing temp_yflux")
	temp_yflux = cc.querying.getvar(exp,'temp_yflux_adv',session,start_time=start_time0, end_time=end_time0)
	temp_xflux = cc.querying.getvar(exp,'temp_xflux_adv',session,start_time=start_time0, end_time=end_time0)
	
	# select latitude range:
	temp_yflux = temp_yflux.sel(yu_ocean=lat_range).sel(time=slice(start_time0,end_time0))
	temp_xflux = temp_xflux.sel(yt_ocean=lat_range).sel(time=slice(start_time0,end_time0))
	
	temp_yflux = temp_yflux.sum('st_ocean')
	temp_xflux = temp_xflux.sum('st_ocean')
	
	temp_yflux = temp_yflux.load()
	temp_xflux = temp_xflux.load()
	
	print("Saving uhrho_vhrho")
	# save a long term average of vhrho_nt and uhrho_et:
	outpath = '/g/data/x77/wf4500/ASC_project/uhrho_vhrho_'+start_time+'.nc'
	# check if already exists:
	if os.path.exists(outpath):
	    average_transports = xr.open_dataset(outpath)
	    # extract arrays from dataset:
	    uhrho_et = average_transports.uhrho_et
	    vhrho_nt = average_transports.vhrho_nt
	else:
	    vhrho_nt = cc.querying.getvar(exp,'vhrho_nt',session,start_time=start_time, end_time=end_time)
	    uhrho_et = cc.querying.getvar(exp,'uhrho_et',session,start_time=start_time, end_time=end_time)
	#    vhrho_nt = cc.querying.getvar(exp,'vhrho_nt',session,frequency='1 daily')
	#    uhrho_et = cc.querying.getvar(exp,'uhrho_et',session,frequency='1 daily')
	
	    vhrho_nt = vhrho_nt.sel(yt_ocean=lat_range).sel(time=slice(start_time0,end_time0))
	    uhrho_et = uhrho_et.sel(yt_ocean=lat_range).sel(time=slice(start_time0,end_time0))
	
	    vhrho_nttime=vhrho_nt.time
	    vhrho_nt = vhrho_nt.mean('time')
	    uhrho_et = uhrho_et.mean('time')
	    
	    
	    outpath = '/g/data/x77/wf4500/ASC_project/model_data/access-om2/'+exp+'/Antarctic_cross_slope/Monthly/uhrho_vhrho_'+start_time+'.nc'
	    ds = xr.Dataset({'vhrho_nt': vhrho_nt,'uhrho_et':uhrho_et})
	    ds.to_netcdf(outpath)
	    ds.close()
	print('time being calculated is ...')      
	print(vhrho_nttime)
	import os
	outpath = '/g/data/x77/wf4500/ASC_project/model_data/access-om2/'+exp+'/Antarctic_cross_slope/Monthly/uhrho_vhrho_'+start_time+'.nc'
	ds = xr.open_dataset(outpath)
	vhrho_nt = ds['vhrho_nt']
	uhrho_et = ds['uhrho_et']
	
	# subtract freezing point heat transport:
	yt_ocean = cc.querying.getvar('01deg_jra55v13_ryf9091','yt_ocean',session,n=1)
	dxu = cc.querying.getvar('01deg_jra55v13_ryf9091','dxu',session,n=1)
	dyt = cc.querying.getvar('01deg_jra55v13_ryf9091','dyt',session,n=1)
	# give dxu and dyt correct coordinates:
	dxu.coords['nj'] = yt_ocean.values
	dxu.coords['ni'] = xt_ocean['xt_ocean'].values
	dxu = dxu.rename(({'ni':'xt_ocean', 'nj':'yt_ocean'}))
	dyt.coords['nj'] = yt_ocean.values
	dyt.coords['ni'] = xt_ocean['xt_ocean'].values
	dyt = dyt.rename(({'ni':'xt_ocean', 'nj':'yt_ocean'}))
	# select latitude range:
	dxu = dxu.sel(yt_ocean=lat_range)
	dyt = dyt.sel(yt_ocean=lat_range)
	
	# Note that in newer mom5 versions this could also be done with ty_trans_int_z, 
	# but there is a problem with this diagnostic in older runs, and even
	# using ty_trans, there is a slight difference. Not sure why?
	
	# Note vhrho_nt is v*dz*1035 and is positioned on north centre edge of t-cell.
	# sum in depth:
	vhrho_nt = vhrho_nt.sum('st_ocean')
	uhrho_et = uhrho_et.sum('st_ocean')
	# convert to transport:
	vhrho_nt = vhrho_nt*dxu/rho_0
	uhrho_et = uhrho_et*dyt/rho_0
	
	# overwrite coords, so we can add the freezing point (with uhrho_et and vhrho_nt) without problems:
	yu_ocean = cc.querying.getvar(exp,'yu_ocean',session,n=1)
	yu_ocean = yu_ocean.sel(yu_ocean=lat_range)
	vhrho_nt.coords['yt_ocean'] = yu_ocean.values
	vhrho_nt = vhrho_nt.rename(({'yt_ocean':'yu_ocean'}))
	uhrho_et.coords['xt_ocean'] = xu_ocean.values
	uhrho_et = uhrho_et.rename(({'xt_ocean':'xu_ocean'}))
	
	freezing_point_heat_trans_zonal = cp*rho_0*uhrho_et*temp_freezing
	freezing_point_heat_trans_meridional = cp*rho_0*vhrho_nt*temp_freezing
	
	print("Removing freezing point heat transport")
	# compare both ways:
	temp_yflux = temp_yflux - freezing_point_heat_trans_meridional
	temp_xflux = temp_xflux - freezing_point_heat_trans_zonal
	    
	temp_yflux_with_mask = np.multiply(np.array(temp_yflux),np.array(mask_y_transport))
	temp_xflux_with_mask = np.multiply(np.array(temp_xflux),np.array(mask_x_transport))
	
	temp_yflux_with_mask = temp_yflux_with_mask[0]
	temp_xflux_with_mask =temp_xflux_with_mask[0]
	
	# multiply by mask to get correct direction into or out of isobath contour:
	heat_trans_across_contour = np.zeros((num_points))
	
	# locations for zonal transport, already calculated indices above:
	for count in x_indices:
	    count = int(count)
	    jj = int(np.where(mask_x_transport_numbered==count)[0])
	    ii = int(np.where(mask_x_transport_numbered==count)[1])
	    heat_trans_across_contour[count-1] += temp_xflux_with_mask[jj,ii]
	    
	# locations for meridional transport, already calculated indices above:
	for count in y_indices:
	    count = int(count)
	    jj = int(np.where(mask_y_transport_numbered==count)[0])
	    ii = int(np.where(mask_y_transport_numbered==count)[1])
	    heat_trans_across_contour[count-1] += temp_yflux_with_mask[jj,ii]
	    
	# convert to longitude coordinate and average into 3 degree longitude bins:
	# in degrees:
	bin_width = 3
	bin_spacing = 0.25
	lon_west = -280
	lon_east = 80
	
	# new coordinate and midpoints of longitude bins:
	full_lon_coord = np.arange(lon_west,lon_east+bin_spacing,bin_spacing)
	lon_bin_midpoints = np.arange(lon_west+bin_width/2,lon_east-bin_width/2,bin_spacing)
	n_bin_edges = len(full_lon_coord)
	
	# sum into longitude bins:
	# need to be very careful of loops, we can't just mask over longitude values, but instead pick indices 
	# on the isobath contour and sum continously along contour between defined indices.
	# (i.e. lon_along_contour is not monotonic)
	# find points on contour to define edges of longitude bins:
	bin_edge_indices = np.zeros(n_bin_edges)
	for lon_bin in range(n_bin_edges-1):
	    # find first isobath point that has the right longitude:
	    first_point = np.where(lon_along_contour>=full_lon_coord[lon_bin])[0][0]
	    # then find all other isobath points with the same longitude as that first point:
	    same_lon_points = np.where(lon_along_contour==lon_along_contour[first_point])[0]
	    # we want the most southerly of these points on the same longitude line:
	    bin_edge_indices[lon_bin] = same_lon_points[np.argmin(lat_along_contour[same_lon_points])]
	    
	# define east/west edges:
	bin_edge_indices = bin_edge_indices.astype(int)
	bin_edge_indices_west = bin_edge_indices[:-int(bin_width/bin_spacing)-1]
	bin_edge_indices_east = bin_edge_indices[int(bin_width/bin_spacing):-1]
	n_bins = len(bin_edge_indices_west)
	
	# sum heat transport from isobath coord into new longitude coord:
	cross_slope_heat_trans = np.zeros(n_bins)
	for lon_bin in range(n_bins):
	    heat_trans_this_bin = heat_trans_across_contour[bin_edge_indices_west[lon_bin]:bin_edge_indices_east[lon_bin]]
	    cross_slope_heat_trans[lon_bin] = np.sum(heat_trans_this_bin)
	    
	# find average latitude of each bin, so we can plot back on the isobath:
	lat_bin_midpoints = np.zeros(n_bins)
	for lon_bin in range(n_bins):
	    # find nearest isobath point:
	    lon_index = np.where(lon_along_contour>=lon_bin_midpoints[lon_bin])[0][0]
	    lat_bin_midpoints[lon_bin] = lat_along_contour[lon_index]
	
	
	temp_xflux=temp_xflux[0]
	temp_yflux=temp_yflux[0]
	
	print('Calculating Zonal heat convergence')
	# Need to make sure the zonal boundaries here match exactly with the zonal boundaries used for the 
	# longitude averaging above, by using same bin_edge_indices.
	# Just check if isobath point is on x or y grid. if x, sum up to this point from south.
	# if on y grid, sum zonal transport on nearest u grid point to west.
	
	zonal_heat_trans_west = np.zeros(n_bins)
	for lon_bin in range(n_bins):
	    # west limit:
	    # reset these to False:
	    on_x_grid = False
	    on_y_grid = False
	    # mask_x_transport_numbered etc indexing starts from 1 not 0, so add 1:
	    isobath_west_index = int(bin_edge_indices_west[lon_bin]+1)
	    # check if the point is on the x or y transport grid:
	    if len(np.where(mask_x_transport_numbered==isobath_west_index)[0])>0:
	        on_x_grid = True
	        jj = int(np.where(mask_x_transport_numbered==isobath_west_index)[0])
	        ii = int(np.where(mask_x_transport_numbered==isobath_west_index)[1])
	    elif len(np.where(mask_y_transport_numbered==isobath_west_index)[0])>0:
	        on_y_grid = True
	        jj = int(np.where(mask_y_transport_numbered==isobath_west_index)[0])
	        ii = int(np.where(mask_y_transport_numbered==isobath_west_index)[1])
	    if on_x_grid == True:
	        zonal_heat_trans_west[lon_bin] = np.sum(temp_xflux[:jj,ii])
	    # in this case we want transport half a grid point to the west:
	    elif on_y_grid == True:
	        # careful if ii=0, then we need heat trans from lon=80, because at limit of zonal grid
	        if ii==0:
	            zonal_heat_trans_west[lon_bin] = np.sum(temp_xflux[:jj+1,-1])
	        else:
	            zonal_heat_trans_west[lon_bin] = np.sum(temp_xflux[:jj+1,ii-1])
	
	zonal_heat_trans_east = np.zeros(n_bins)
	for lon_bin in range(n_bins):
	    # east limit:
	    # reset these to False:
	    on_x_grid = False
	    on_y_grid = False
	    # mask_x_transport_numbered etc indexing starts from 1 not 0, so add 1:
	    isobath_east_index = int(bin_edge_indices_east[lon_bin]+1)
	    # check if the point is on the x or y transport grid:
	    if len(np.where(mask_x_transport_numbered==isobath_east_index)[0])>0:
	        on_x_grid = True
	        jj = int(np.where(mask_x_transport_numbered==isobath_east_index)[0])
	        ii = int(np.where(mask_x_transport_numbered==isobath_east_index)[1])
	    elif len(np.where(mask_y_transport_numbered==isobath_east_index)[0])>0:
	        on_y_grid = True
	        jj = int(np.where(mask_y_transport_numbered==isobath_east_index)[0])
	        ii = int(np.where(mask_y_transport_numbered==isobath_east_index)[1])
	    if on_x_grid == True:
	        zonal_heat_trans_east[lon_bin] = np.sum(temp_xflux[:jj,ii])
	    # in this case we want transport half a grid point to the west:
	    elif on_y_grid == True:
	        # wrap around to east side of grid:
	        if ii==0:
	            zonal_heat_trans_east[lon_bin] = np.sum(temp_xflux[:jj+1,-1])
	        else:
	            zonal_heat_trans_east[lon_bin] = np.sum(temp_xflux[:jj+1,ii-1])
	
	zonal_convergence = zonal_heat_trans_east - zonal_heat_trans_west
	
	# convert to data arrays, so we can save as netcdf:
	zonal_convergence = xr.DataArray(zonal_convergence, coords = [('lon_bin_midpoints', lon_bin_midpoints)])
	cross_slope_heat_trans = xr.DataArray(cross_slope_heat_trans, coords = [('lon_bin_midpoints', lon_bin_midpoints)])
	heat_trans_across_contour = xr.DataArray(heat_trans_across_contour,coords = [('lon_along_contour',lon_along_contour)])
	
	
	#After running this time, take off the 2 in the end of the .nc file
	#outpath = '/g/data/x77/wf4500/ASC_project/Cross_slope_transport/Montlhy/Ant_cross_slope_heat_terms_online_'+str(isobath_depth)+'m_'+start_str+'.nc'
	outpath = '/g/data/x77/wf4500/ASC_project/cross_slope_transport/Monthly/' + exp + '/Ant_cross_slope_heat_terms_online_'+str(isobath_depth)+'m_'+ start_time+'.nc'
	
	print('Saving heat transports...')
	ds = xr.Dataset({'zonal_convergence': zonal_convergence,'cross_slope_heat_trans':cross_slope_heat_trans,
	                'lon_bin_midpoints':lon_bin_midpoints,'lat_bin_midpoints':lat_bin_midpoints,
	                 'bin_width':bin_width,'bin_spacing':bin_spacing,'heat_trans_across_contour':heat_trans_across_contour,
	                'lon_along_contour':lon_along_contour})
	ds.to_netcdf(outpath)
	
	print('Finished successful')
	quit()    
