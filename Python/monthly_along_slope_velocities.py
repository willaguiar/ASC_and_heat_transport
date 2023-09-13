# # Compute monthly ASC along the U-grid contour and saves it



import cosima_cookbook as cc
import matplotlib.pyplot as plt
import netCDF4 as nc
import xarray as xr
import numpy as np
import sys,os
from dask.distributed import Client
#from cosima_cookbook import distributed as ccd
# Optional modules
#import xarray.ufuncs as xu
import xgcm

import climtas.nci

import warnings # ignore these warnings
warnings.filterwarnings("ignore", category = FutureWarning)
warnings.filterwarnings("ignore", category = UserWarning)
warnings.filterwarnings("ignore", category = RuntimeWarning)

if __name__ == '__main__':

	climtas.nci.GadiClient()
	# Load database
	#DO NOT FORGET TO SET THE EXPERIMENT NAME (EXPT) AND TO ALTER LINES [44,51] AND [103,113]
	session = cc.database.create_session()
	# Define experiment in database
	expt = '01deg_jra55v13_ryf9091'
	
	#expt0 = '01deg_jra55v140_iaf' #Line below only necessary for iaf. if ryf repeat name
	
	lat_slice = slice(-80, -59)
	
	# Import bathymetry
	hu = cc.querying.getvar(expt, 'hu', session, n=1)
	
	# Import grid cell length
	dxu = cc.querying.getvar(expt, 'dxu', session, n=1)
	dyu = cc.querying.getvar(expt, 'dyu', session, n=1)
	
	##THis part below doesnt need to be run for IAF experiments, BUT SHOULD BE UNCOMMENTED FOR RYF EXPERIMENTS  
	## Change coordinate name
	dxu.coords['ni'] = hu['xu_ocean'].values
	dxu.coords['nj'] = hu['yu_ocean'].values
	dyu.coords['ni'] = hu['xu_ocean'].values
	dyu.coords['nj'] = hu['yu_ocean'].values
	dxu = dxu.rename(({'ni':'xu_ocean', 'nj':'yu_ocean'}))
	dyu = dyu.rename(({'ni':'xu_ocean', 'nj':'yu_ocean'}))
	
	# Select latitude range
	hu = hu.sel(yu_ocean=lat_slice)
	hu = hu.load()
	dxu = dxu.sel(yu_ocean=lat_slice)
	dxu = dxu.load()
	dyu = dyu.sel(yu_ocean=lat_slice)
	dyu = dyu.load()
	
	# Load model grid information
	path_to_folder = '/g/data/ik11/outputs/access-om2-01/01deg_jra55v13_ryf9091/output000/ocean/'
	grid = xr.open_mfdataset(path_to_folder+'ocean_grid.nc', combine='by_coords')

	ds = xr.merge([hu, grid])
	ds.coords['xt_ocean'].attrs.update(axis='X')
	ds.coords['xu_ocean'].attrs.update(axis='X', c_grid_axis_shift=0.5)
	ds.coords['yt_ocean'].attrs.update(axis='Y')
	ds.coords['yu_ocean'].attrs.update(axis='Y', c_grid_axis_shift=0.5)
	
	grid = xgcm.Grid(ds, periodic=['X'])
	
	# Take gradient and move to u grid
	# Simple gradient over one grid cell. 
	# In latitudinal direction, we need to specify what happens at the boundary.
	dhu_dx = grid.interp( grid.diff(ds.hu, 'X') / grid.interp(ds.dxu, 'X'), 'X')#, 'Y', boundary='extend')
	dhu_dy = grid.interp( grid.diff(ds.hu, 'Y', boundary='extend') / grid.interp(ds.dyt, 'X'), 'Y', boundary='extend')# 'X')
	
	# Same as: (other option)
	#dhu_dx = grid.interp(grid.diff(ds.hu, 'X'), 'X') / ds.dxu
	#dhu_dy = grid.interp(grid.diff(ds.hu, 'Y', boundary='extend'), 'Y', boundary='extend') / ds.dyu
	
	# Select latitude slice
	dhu_dx = dhu_dx.sel(yu_ocean=lat_slice)
	dhu_dy = dhu_dy.sel(yu_ocean=lat_slice)
	
	# Calculate the magnitude of the topographic slope
	slope = (dhu_dx**2 + dhu_dy**2)**0.5 #before it was : xu.sqrt(dhu_dx**2 + dhu_dy**2)
	
	# Load velocity data
	#yearinit=2177
	yearinit = sys.argv[1]
	monthseries=['01','02','03','04','05','06','07','08','09','10','11','12']
	mi=0
	isobath_depth = 1000
	test = np.int(yearinit)+1
#	print(test)  
# 	print('isobath depth is ' + str(isobath_depth) + ' m')
	for n in np.arange(0,12): 
		start_time=str(yearinit) + '-' + monthseries[n]
		if n <= 10:
		    end_time=str(yearinit) + '-' + monthseries[n+1]
		    
		if n >= 11:   
		    end_time=str(np.int(yearinit)+1) + '-01'
		
		
		print(end_time)
#		if n >= 10:
#		    end_time=str(int(np.int((yearinit)+1)) + '-01'
#		    print(str(n))        
#		
#		
#		     end_time=str(yearinit) + '-' + monthseries[n+1]       
#		else:
#		    end_time=str(int(np.int((yearinit)+1)) + '-01'
#		
#		print('End time is ' + end_time)
#		
		u = cc.querying.getvar(expt, 'u', session, start_time=start_time, end_time=end_time, ncfile="ocean.nc")
		u = u.sel(time=slice(start_time,end_time)).sel(yu_ocean=lat_slice).isel(time=0)
		v = cc.querying.getvar(expt, 'v', session, start_time=start_time, end_time=end_time, ncfile="ocean.nc")
		v = v.sel(time=slice(start_time,end_time)).sel(yu_ocean=lat_slice).isel(time=0)
		
# 	    #     #if IAF run lines this way
# 	    # u = cc.querying.getvar(expt, 'u', session, start_time=start_time, end_time=end_time, frequency='1 monthly')
# 	    # u = u.sel(time=slice(start_time,end_time)).sel(yu_ocean=lat_slice).isel(time=0)
# 	    # v = cc.querying.getvar(expt, 'v', session, start_time=start_time, end_time=end_time, frequency='1 monthly')
# 	    # v = v.sel(time=slice(start_time,end_time)).sel(yu_ocean=lat_slice).isel(time=0)
		
		# Along-slope velocity
		u_along = u*dhu_dy/slope - v*dhu_dx/slope
		# Cross-slope velocity
		v_cross = u*dhu_dx/slope + v*dhu_dy/slope
		 
		
		# import edges of st_ocean and add lat/lon dimensions:
		st_edges_ocean = cc.querying.getvar(expt, 'st_edges_ocean', session, start_time=start_time, end_time=end_time, n=1)
		st_edges_array = st_edges_ocean.expand_dims({'yu_ocean':u.yu_ocean,'xu_ocean':u.xu_ocean}, axis=[1,2])
		
		# adjust edges at bottom for partial thickness:
		st_edges_with_partial = st_edges_array.where(st_edges_array<hu, other=hu)
		thickness = st_edges_with_partial.diff(dim='st_edges_ocean')
		
		# change coordinate of thickness to st_ocean (needed for multipling with other variables):
		st_ocean = cc.querying.getvar(expt, 'st_ocean', session, n=1)
		thickness['st_edges_ocean'] = st_ocean.values
		thickness = thickness.rename(({'st_edges_ocean':'st_ocean'}))
		thickness = thickness
		
		# Barotropic velocity
		u_btr = (u_along * thickness).sum('st_ocean') / thickness.sum('st_ocean')
		v_btr = (v_cross * thickness).sum('st_ocean') / thickness.sum('st_ocean')
		
		# Load the data (this is necessary for selecting the points along the contour, otherwise the code will run forever)
		u_along = u_along.load()
		u_btr   = u_btr.load()
		
		v_cross = v_cross.load()
		v_btr   = v_btr.load()
		
		# Load the isobath information
		file_path = '/g/data/x77/wgh581/Post_Process/Antarctic_slope_contour_u_t_' + str(isobath_depth) + 'm_plt_approach.npz'
		data = np.load(file_path)
		
		# We are interested in the u-grid:
		num_points = data['num_points_u']
		distance_along_contour = data['distance_along_contour_u']
		lon_along_contour = data['lon_along_contour_u']
		lat_along_contour = data['lat_along_contour_u']
		
		# Total velocity along contour
		u_along_contour = np.zeros((1, len(u_along.st_ocean), num_points))
		for ii in range(num_points):
		        u_along_contour[:,:,ii] = u_along.sel(xu_ocean=lon_along_contour[ii], method='nearest').sel(yu_ocean=lat_along_contour[ii], method='nearest')
		
		# Barotropic velocity along contour
		u_btr_contour = np.zeros((1,num_points))
		for ii in range(num_points):
		    u_btr_contour[:,ii] = u_btr.sel(xu_ocean=lon_along_contour[ii], method='nearest').sel(yu_ocean=lat_along_contour[ii], method='nearest')
		 
		# Total velocity along contour
		v_cross_contour = np.zeros((1, len(v_cross.st_ocean), num_points))
		for ii in range(num_points):
		    v_cross_contour[:,:,ii] = v_cross.sel(xu_ocean=lon_along_contour[ii], method='nearest').sel(yu_ocean=lat_along_contour[ii], method='nearest')
		 
		# Barotropic velocity along contour
		v_btr_contour = np.zeros((1,num_points))
		for ii in range(num_points):
		    v_btr_contour[:,ii] = v_btr.sel(xu_ocean=lon_along_contour[ii], method='nearest').sel(yu_ocean=lat_along_contour[ii], method='nearest')
		
		thickness_contour = np.zeros((len(thickness.st_ocean), num_points))
		for ii in range(num_points):
		    thickness_contour[:,ii] = thickness.sel(xu_ocean=lon_along_contour[ii], method='nearest').sel(yu_ocean=lat_along_contour[ii], method='nearest')
		
		# PREVIOUS SAVEDIR TO LOCATE THE FILES save_dir  = '/g/data/x77/wf4500/ASC_project/ASC_speed/monthly/OM2_RYF/'
		save_dir  = 'add here dir you want to save the files'        
		file_name = 'Antarctic_slope_contour_1km_velocities_'
		
		data_u=xr.DataArray((u_along_contour[0]),name="u_along_contour",dims=["st_ocean","distance_along_contour"])
		data_v=xr.DataArray((v_cross_contour[0]),name="v_cross_contour",dims=["st_ocean","distance_along_contour"])
		data_ubt=xr.DataArray((u_btr_contour[0]),name="u_btr_contour",dims=["distance_along_contour"])
		data_vbt=xr.DataArray((v_btr_contour[0]),name="v_btr_contour",dims=["distance_along_contour"])
		data_lon=xr.DataArray((lon_along_contour),name="lon_along_contour",dims=["distance_along_contour"])
		data_lat=xr.DataArray((lat_along_contour),name="lat_along_contour",dims=["distance_along_contour"])
		data_thick=xr.DataArray((thickness_contour),name="thickness_contour",dims=["st_ocean","distance_along_contour"])
		
		data_to_saveuv = xr.merge([data_u,data_v,data_ubt,data_vbt,data_lon,data_lat,data_thick])
		data_to_saveuv['distance_along_contour']=np.array(distance_along_contour)
		data_to_saveuv['st_ocean']=np.array(st_ocean)
		data_to_saveuv
		data_to_saveuv.to_netcdf(save_dir + file_name + start_time +"_uv.nc")
