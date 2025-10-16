# # Compute daily ASC along contour, on x,y grid and save in netcdfs for each month. U,V converted from uhrho and vhrho



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
import geopy.distance

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
	expt = '01deg_jra55v140_iaf_cycle3'
	
	#expt0 = '01deg_jra55v140_iaf' #Line below only necessary for iaf. if ryf repeat name
	
	lat_slice = slice(-80, -57.6)
	cutout_latind=510
	
	# Import bathymetry
	hu = cc.querying.getvar(expt, 'hu', session, n=1)
	
	# Import grid cell length
	dxu = cc.querying.getvar(expt, 'dxu', session, n=1)
	dyu = cc.querying.getvar(expt, 'dyu', session, n=1)
	
	##THis part below doesnt need to be run for IAF experiments, BUT SHOULD BE UNCOMMENTED FOR RYF EXPERIMENTS  
	## Change coordinate name
	dxu.coords['xu_ocean'] = hu['xu_ocean'].values                            
	dxu.coords['yu_ocean'] = hu['yu_ocean'].values                            
	dyu.coords['xu_ocean'] = hu['xu_ocean'].values                            
	dyu.coords['yu_ocean'] = hu['yu_ocean'].values                            

	
	# Select latitude range
	hu = hu.sel(yu_ocean=lat_slice)
	hu = hu.load()
	dxu = dxu.sel(yu_ocean=lat_slice)
	dxu = dxu.load()
	dyu = dyu.sel(yu_ocean=lat_slice)
	dyu = dyu.load()
	
	# Load model grid information
	path_to_folder = '/g/data/ik11/outputs/access-om2-01/01deg_jra55v13_ryf9091/output000/ocean/'
	grid = xr.open_mfdataset('/g/data/x77/wf4500/ASC_project/GH/GH_AP2025/ASC_and_heat_transport/Python/OM2_grid/ocean_grid.nc', combine='by_coords')

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
	

	
	# Select latitude slice
	dhu_dx = dhu_dx.sel(yu_ocean=lat_slice)
	dhu_dy = dhu_dy.sel(yu_ocean=lat_slice)
	
	# Calculate the magnitude of the topographic slope
	slope = (dhu_dx**2 + dhu_dy**2)**0.5 #before it was : xu.sqrt(dhu_dx**2 + dhu_dy**2)
	
	    
	    
	 
	###################################################################  
	#interpolating the slope to x,y grid    
	Slopedata=xr.open_dataset('/g/data/x77/wf4500/ASC_project/Post_process/slope1000m.nc') 
	x_s=slope.xu_ocean
	y_s=slope.yu_ocean
	
	# Load model grid information
	path_to_folder = '/g/data/ik11/outputs/access-om2-01/01deg_jra55v13_ryf9091/output000/ocean/'
	gridx = xr.open_mfdataset('/g/data/x77/wf4500/ASC_project/GH/GH_AP2025/ASC_and_heat_transport/Python/OM2_grid/ocean_grid.nc', combine='by_coords')
	
	dsx = xr.merge([Slopedata.slope, gridx])
	dsx.coords['xt_ocean'].attrs.update(axis='X')
	dsx.coords['xu_ocean'].attrs.update(axis='X', c_grid_axis_shift=0.5)
	dsx.coords['yt_ocean'].attrs.update(axis='Y')
	dsx.coords['yu_ocean'].attrs.update(axis='Y', c_grid_axis_shift=0.5)
	
	gridx = xgcm.Grid(dsx, periodic=['X'])
	
	
	#y
	gridy = xr.open_mfdataset('/g/data/x77/wf4500/ASC_project/GH/GH_AP2025/ASC_and_heat_transport/Python/OM2_grid/ocean_grid.nc', combine='by_coords')
	dsy = xr.merge([Slopedata.slope, gridy])
	dsy.coords['xt_ocean'].attrs.update(axis='X')
	dsy.coords['xu_ocean'].attrs.update(axis='X', c_grid_axis_shift=0.5)
	dsy.coords['yt_ocean'].attrs.update(axis='Y')
	dsy.coords['yu_ocean'].attrs.update(axis='Y', c_grid_axis_shift=0.5)
	
	gridy = xgcm.Grid(dsy, periodic=['X'])
	
	# Calculating the slope on the x and y grid
	slope_xg=gridx.interp(dsx.slope, 'Y').isel(yt_ocean=slice(0,cutout_latind))
	slope_yg=gridx.interp(dsx.slope, 'X').isel(yu_ocean=slice(0,cutout_latind))
	    
	
	###################################################################    
	#Interpolating the dhu into the x and y grids
	x_s=slope.xu_ocean
	y_s=slope.yu_ocean
	# Load model grid information
	path_to_folder = '/g/data/ik11/outputs/access-om2-01/01deg_jra55v13_ryf9091/output000/ocean/'
	gridx = xr.open_mfdataset('/g/data/x77/wf4500/ASC_project/GH/GH_AP2025/ASC_and_heat_transport/Python/OM2_grid/ocean_grid.nc', combine='by_coords')
	
	dsx = xr.merge([Slopedata.dhu_dy,Slopedata.dhu_dx, gridx])
	dsx.coords['xt_ocean'].attrs.update(axis='X')
	dsx.coords['xu_ocean'].attrs.update(axis='X', c_grid_axis_shift=0.5)
	dsx.coords['yt_ocean'].attrs.update(axis='Y')
	dsx.coords['yu_ocean'].attrs.update(axis='Y', c_grid_axis_shift=0.5)
	gridx = xgcm.Grid(dsx, periodic=['X'])
	
	
	#y
	gridy = xr.open_mfdataset('/g/data/x77/wf4500/ASC_project/GH/GH_AP2025/ASC_and_heat_transport/Python/OM2_grid/ocean_grid.nc', combine='by_coords')
	dsy = xr.merge([Slopedata.dhu_dy,Slopedata.dhu_dx, gridy])
	dsy.coords['xt_ocean'].attrs.update(axis='X')
	dsy.coords['xu_ocean'].attrs.update(axis='X', c_grid_axis_shift=0.5)
	dsy.coords['yt_ocean'].attrs.update(axis='Y')
	dsy.coords['yu_ocean'].attrs.update(axis='Y', c_grid_axis_shift=0.5)
	gridy = xgcm.Grid(dsy, periodic=['X'])
	
	# Calculating the slope on the x and y grid
	dhu_dy_xg=gridx.interp(dsx.dhu_dy, 'Y').isel(yt_ocean=slice(0,cutout_latind))
	dhu_dy_yg=gridx.interp(dsx.dhu_dy, 'X').isel(yu_ocean=slice(0,cutout_latind))
	dhu_dx_xg=gridx.interp(dsx.dhu_dx, 'Y').isel(yt_ocean=slice(0,cutout_latind))
	dhu_dx_yg=gridx.interp(dsx.dhu_dx, 'X').isel(yu_ocean=slice(0,cutout_latind))
	
	
	
	  
	# Load velocity data
	#yearinit=2177
	yearinit = sys.argv[2]
	imonth = int(sys.argv[1])
	month = str(int(sys.argv[1])) 
	month = month.zfill(2)
	monthseries=['01','02','03','04','05','06','07','08','09','10','11','12']
	dayseries=['31','28','31','30','31','30','31','31','30','31','30','31'] 
	mi=0
	isobath_depth = 1000
	test = int(yearinit)+1
	#Accounting for Leap years
	long_feby=np.arange(1960,2028,4)
	if yearinit in long_feby: dayseries=['31','29','31','30','31','30','31','31','30','31','30','31'] ; print('this is a leap year')
#	print(test)  
# 	print('isobath depth is ' + str(isobath_depth) + ' m')
	for n in np.arange(0,1): 
        
		start_time=str(yearinit) + '-' + month + '-01'
		end_time=str(yearinit) + '-' + month + '-' + dayseries[imonth-1]

		
		print(start_time)		
		print(end_time)
#		
		uhrho = cc.querying.getvar(expt, 'uhrho_et', session,start_time=start_time, end_time=end_time,frequency='1 daily')
		uhrho = uhrho.sel(time=slice(start_time,end_time)).sel(yt_ocean=lat_slice)
		vhrho = cc.querying.getvar(expt, 'vhrho_nt', session,start_time=start_time, end_time=end_time,frequency='1 daily')
		vhrho = vhrho.sel(time=slice(start_time,end_time)).sel(yt_ocean=lat_slice)

		rho=1035
		ht = cc.querying.getvar(expt, 'dzt', session,start_time=start_time, end_time=end_time,frequency='1 monthly')
		ht=ht.sel(time=slice(start_time,end_time)).sel(yt_ocean=lat_slice).isel(time=0,drop=True)	

		u=(uhrho/ht)/rho
		v=(vhrho/ht)/rho
		u.name='u'
		v.name='v'
		
		# #reducing u,v time samples for testing - remove this line after fully testing the data
		# u=u.isel(time=slice(0,2))
		# v=v.isel(time=slice(0,2))
		
		
		####################Interpolating U to x,y grid!!
		gridx = xr.open_mfdataset('/g/data/x77/wf4500/ASC_project/GH/GH_AP2025/ASC_and_heat_transport/Python/OM2_grid/ocean_grid.nc', combine='by_coords')
		gridx=gridx.isel(time=0,drop=True)                           
		dsx = xr.merge([u, gridx])
		dsx.coords['xt_ocean'].attrs.update(axis='X')
		dsx.coords['xu_ocean'].attrs.update(axis='X', c_grid_axis_shift=0.5)
		dsx.coords['yt_ocean'].attrs.update(axis='Y')
		dsx.coords['yu_ocean'].attrs.update(axis='Y', c_grid_axis_shift=0.5)
		dsx=dsx.sel(time=slice(start_time,end_time))
		
		gridx = xgcm.Grid(dsx, periodic=['X'])
		
		
		#y
		gridy = xr.open_mfdataset('/g/data/x77/wf4500/ASC_project/GH/GH_AP2025/ASC_and_heat_transport/Python/OM2_grid/ocean_grid.nc', combine='by_coords')
		gridy=gridy.isel(time=0,drop=True)                          
		dsy = xr.merge([u, gridy])
		dsy.coords['xt_ocean'].attrs.update(axis='X')
		dsy.coords['xu_ocean'].attrs.update(axis='X', c_grid_axis_shift=0.5)
		dsy.coords['yt_ocean'].attrs.update(axis='Y')
		dsy.coords['yu_ocean'].attrs.update(axis='Y', c_grid_axis_shift=0.5)
		dsy=dsy.sel(time=slice(start_time,end_time)) 
		
		gridy = xgcm.Grid(dsy, periodic=['X'])
		
		# Calculating the slope on the x and y grid
		u_xg=gridx.interp(dsx.u, 'X').isel(yt_ocean=slice(0,cutout_latind))
		u_yg=gridx.interp(dsy.u, 'Y').isel(yu_ocean=slice(0,cutout_latind))
		        
		        
		
		        
		        
		        
		####################Interpolating V to x,y grid!!
		gridx = xr.open_mfdataset('/g/data/x77/wf4500/ASC_project/GH/GH_AP2025/ASC_and_heat_transport/Python/OM2_grid/ocean_grid.nc', combine='by_coords')
		gridx=gridx.isel(time=0,drop=True)                       
		dsx = xr.merge([v, gridx])
		dsx.coords['xt_ocean'].attrs.update(axis='X')
		dsx.coords['xu_ocean'].attrs.update(axis='X', c_grid_axis_shift=0.5)
		dsx.coords['yt_ocean'].attrs.update(axis='Y')
		dsx.coords['yu_ocean'].attrs.update(axis='Y', c_grid_axis_shift=0.5)
		dsx=dsx.sel(time=slice(start_time,end_time))
		gridx = xgcm.Grid(dsx, periodic=['X'])
		
		
		#y
		gridy = xr.open_mfdataset('/g/data/x77/wf4500/ASC_project/GH/GH_AP2025/ASC_and_heat_transport/Python/OM2_grid/ocean_grid.nc', combine='by_coords')
		gridy=gridy.isel(time=0,drop=True)                     
		dsy = xr.merge([v, gridy])
		dsy.coords['xt_ocean'].attrs.update(axis='X')
		dsy.coords['xu_ocean'].attrs.update(axis='X', c_grid_axis_shift=0.5)
		dsy.coords['yt_ocean'].attrs.update(axis='Y')
		dsy.coords['yu_ocean'].attrs.update(axis='Y', c_grid_axis_shift=0.5)
		dsy=dsy.sel(time=slice(start_time,end_time))
		gridy = xgcm.Grid(dsy, periodic=['X'])
		
		# Calculating the slope on the x and y grid
		v_xg=gridx.interp(dsx.v, 'X').isel(yt_ocean=slice(0,cutout_latind))
		v_yg=gridx.interp(dsy.v, 'Y').isel(yu_ocean=slice(0,cutout_latind))
		        
		        
		
		############################ Along-slope velocity
		Uxg_along = u_xg*dhu_dy_xg/slope_xg - v_xg*dhu_dx_xg/slope_xg
		Uyg_along = u_yg*dhu_dy_yg/slope_yg - v_yg*dhu_dx_yg/slope_yg
		############################ Cross-slope velocity       
		Vxg_cross = u_xg*dhu_dx_xg/slope_xg + v_xg*dhu_dy_xg/slope_xg        
		Vyg_cross = u_yg*dhu_dx_yg/slope_yg + v_yg*dhu_dy_yg/slope_yg
		
		 
		
		# import edges of st_ocean and add lat/lon dimensions:
		st_edges_ocean = cc.querying.getvar(expt, 'st_edges_ocean', session, start_time=start_time, end_time=end_time, n=1)
		st_edges_array = st_edges_ocean.expand_dims({'yu_ocean':hu.yu_ocean,'xu_ocean':hu.xu_ocean}, axis=[1,2])
		
		# adjust edges at bottom for partial thickness:
		st_edges_with_partial = st_edges_array.where(st_edges_array<hu, other=hu)
		thickness = st_edges_with_partial.diff(dim='st_edges_ocean')
		
		# change coordinate of thickness to st_ocean (needed for multipling with other variables):
		st_ocean = cc.querying.getvar(expt, 'st_ocean', session, n=1)
		thickness['st_edges_ocean'] = st_ocean.values
		thickness = thickness.rename(({'st_edges_ocean':'st_ocean'}))
		thickness = thickness
		
		         
		
		                
		                
		########### Import Adeles 1km contour on X,Y grid                
		#outfile = '/g/data/ik11/grids/Antarctic_slope_contour_1000m.npz'
		# outfile = '/g/data/ik11/grids/Antarctic_slope_contour_1000m.npz'
		# data = np.load(outfile)
		data = xr.open_dataset('/g/data/x77/wf4500/ASC_project/GH/GH_AP2025/ASC_and_heat_transport/Python/OM2_grid/Antarctic_slope_contour_650m_MOM6_01deg.nc')
		data = data.rename({'xh':'xt_ocean','yh':'yt_ocean','xq':'xu_ocean','yq':'yu_ocean'})
		mask_y_transport = data['mask_y_transport']; 
		mask_y_transport=mask_y_transport[:cutout_latind,:]
		mask_x_transport = data['mask_x_transport']; 
		mask_x_transport=mask_x_transport[:cutout_latind,:]
		mask_y_transport_numbered = data['mask_y_transport_numbered']
		mask_y_transport_numbered=mask_y_transport_numbered[:cutout_latind,:]
		mask_x_transport_numbered = data['mask_x_transport_numbered']
		mask_x_transport_numbered=mask_x_transport_numbered[:cutout_latind,:]
		
		
		#cutting the matrices
		ylength= np.shape(mask_x_transport)[0]
		
		
		yt_ocean0 = cc.querying.getvar(expt,'yt_ocean',session,n=1)
		yt_ocean0 = yt_ocean0.isel(yt_ocean=slice(0,ylength))
		yu_ocean0 = cc.querying.getvar(expt,'yu_ocean',session,n=1)
		yu_ocean0 = yu_ocean0.isel(yu_ocean=slice(0,ylength))
		xt_ocean0 = cc.querying.getvar(expt,'xt_ocean',session,n=1)
		xu_ocean0 = cc.querying.getvar(expt,'xu_ocean',session,n=1)
		
		
		#changing mask X transport
		# mask_x_transport =xr.DataArray(data['mask_x_transport']).assign_coords({"dim_0": np.array(yt_ocean0),"dim_1": np.array(xu_ocean0)}).rename(dim_0="y_ocean",dim_1="x_ocean")
		mask_x_transport = mask_x_transport.rename(yt_ocean="y_ocean",xu_ocean="x_ocean")
		# mask_x_transport['x_ocean'] = xu_ocean0.values ; mask_x_transport['y_ocean'] = yt_ocean0.values
		mask_x_transport=mask_x_transport[:cutout_latind,:]       ##############################
		
		
		#changing mask y transport
		# mask_y_transport =xr.DataArray(data['mask_y_transport']).assign_coords({"dim_0": np.array(yt_ocean0),"dim_1": np.array(xu_ocean0)}).rename(dim_0="y_ocean",dim_1="x_ocean")
		mask_y_transport = mask_y_transport.rename(yu_ocean="y_ocean",xt_ocean="x_ocean")
		# mask_y_transport['x_ocean'] = xt_ocean0.values ; mask_y_transport['y_ocean'] = yu_ocean0.values
		mask_y_transport=mask_y_transport[:cutout_latind,:]
		
		#changing mask X transport numbered
		# mask_x_transport_numbered =xr.DataArray(data['mask_x_transport_numbered']).assign_coords({"dim_0": np.array(yt_ocean0),"dim_1": np.array(xu_ocean0)}).rename(dim_0="y_ocean",dim_1="x_ocean")
		mask_x_transport_numbered = mask_x_transport_numbered.rename(yt_ocean="y_ocean",xu_ocean="x_ocean")
		# mask_x_transport_numbered['x_ocean'] = xu_ocean0.values ; mask_x_transport_numbered['y_ocean'] = yt_ocean0.values
		mask_x_transport_numbered=mask_x_transport_numbered[:cutout_latind,:]
		
		#changing mask Y transport numbered
		# mask_y_transport_numbered =xr.DataArray(data['mask_y_transport_numbered']).assign_coords({"dim_0": np.array(yu_ocean0),"dim_1": np.array(xt_ocean0)}).rename(dim_0="y_ocean",dim_1="x_ocean")
		mask_y_transport_numbered = mask_y_transport_numbered.rename(yu_ocean="y_ocean",xt_ocean="x_ocean")
		# mask_y_transport_numbered ['x_ocean'] = xt_ocean0.values ; mask_y_transport_numbered ['y_ocean'] = yu_ocean0.values
		mask_y_transport_numbered=mask_y_transport_numbered[:cutout_latind,:]
		
		num_points = int(np.maximum(np.max(mask_y_transport_numbered),np.max(mask_x_transport_numbered)))
		
		
		
		
		# Create the contour order data-array. Note that in this procedure the x-grid counts have x-grid
		#   dimensions and the y-grid counts have y-grid dimensions, but these are implicit, the dimension 
		#   *names* are kept general across the counts, the generic y_ocean, x_ocean, so that concatening works
		#   but we dont double up with numerous counts for one lat/lon point.
		
		# stack contour data into 1d:
		mask_x_numbered_1d = mask_x_transport_numbered.stack(contour_index = ['y_ocean', 'x_ocean'])
		mask_x_numbered_1d = mask_x_numbered_1d.where(mask_x_numbered_1d > 0, drop = True)
		
		mask_y_numbered_1d = mask_y_transport_numbered.stack(contour_index = ['y_ocean', 'x_ocean'])
		mask_y_numbered_1d = mask_y_numbered_1d.where(mask_y_numbered_1d > 0, drop = True)
		
		contour_ordering = xr.concat((mask_x_numbered_1d, mask_y_numbered_1d), dim = 'contour_index')
		contour_ordering = contour_ordering.sortby(contour_ordering)
		contour_index_array = np.arange(1, len(contour_ordering)+1)
		
		
		############### Getting the lon,lat along contour in the X,Y contour
		lat_along_contour = np.zeros((num_points))
		lon_along_contour = np.zeros((num_points))
		# locations for zonal transport:
		x_indices_masked = mask_x_transport_numbered.stack().values
		x_indices = np.sort(x_indices_masked[x_indices_masked>0])
		for count in x_indices:
		    count = int(count)
		    jj = int(np.where(mask_x_transport_numbered==count)[0])
		    ii = int(np.where(mask_x_transport_numbered==count)[1])   
		    lon_along_contour[count-1] = data.xu_ocean[ii].values
		    lat_along_contour[count-1] = mask_x_transport_numbered.y_ocean[jj].values
		    
		# locations for meridional transport:
		y_indices_masked = mask_y_transport_numbered.stack().values
		y_indices = np.sort(y_indices_masked[y_indices_masked>0])
		for count in y_indices:
		    count = int(count)
		    jj = np.where(mask_y_transport_numbered==count)[0]
		    ii = np.where(mask_y_transport_numbered==count)[1]
		    lon_along_contour[count-1] = mask_x_transport_numbered.x_ocean[ii].values
		    lat_along_contour[count-1] = data.yu_ocean[jj].values
		
		
		print('loading transposed Us')		
		# Importing the Along-slope speed on x and y grid
		Uyg_along0=Uyg_along.load()
		Uxg_along0=Uxg_along.load()
		Vyg_cross0=Vyg_cross.load()
		Vxg_cross0=Vxg_cross.load()

		Uxg_along0=Uxg_along0.rename(yt_ocean='y_ocean',xu_ocean='x_ocean')
		Uyg_along0=Uyg_along0.rename(yu_ocean='y_ocean',xt_ocean='x_ocean')
		Vxg_cross0 = Vxg_cross0.rename(yt_ocean='y_ocean',xu_ocean='x_ocean')
		Vyg_cross0 = Vyg_cross0.rename(yu_ocean='y_ocean',xt_ocean='x_ocean')
		 
		

		# stack transports into 1d and drop any points not on contour:
		x_along_1d = Uxg_along0.stack(contour_index = ['y_ocean', 'x_ocean'])
		x_along_1d = x_along_1d.where(mask_x_numbered_1d>0, drop = True)
		y_along_1d = Uyg_along0.stack(contour_index = ['y_ocean', 'x_ocean'])
		y_along_1d = y_along_1d.where(mask_y_numbered_1d>0, drop = True)

		# combine all points on contour:
		Ualong = xr.concat((x_along_1d, y_along_1d), dim = 'contour_index')
		Ualong = Ualong.sortby(contour_ordering)
		Ualong.coords['contour_index'] = contour_index_array
		Ualong = Ualong.load()	


		# stack transports into 1d and drop any points not on contour:
		x_cross_1d = Vxg_cross0.stack(contour_index = ['y_ocean', 'x_ocean'])
		x_cross_1d = x_cross_1d.where(mask_x_numbered_1d>0, drop = True)
		y_cross_1d = Vyg_cross0.stack(contour_index = ['y_ocean', 'x_ocean'])
		y_cross_1d = y_cross_1d.where(mask_y_numbered_1d>0, drop = True)
		
		# combine all points on contour:
		Ucross = xr.concat((x_cross_1d, y_cross_1d), dim = 'contour_index')
		Ucross = Ucross.sortby(contour_ordering)
		Ucross.coords['contour_index'] = contour_index_array
		Ucross = Ucross.load()

		 
		#Interpolating Thickness to X,Y grid
		gridx = xr.open_mfdataset('/g/data/x77/wf4500/ASC_project/GH/GH_AP2025/ASC_and_heat_transport/Python/OM2_grid/ocean_grid.nc', combine='by_coords')
		gridx=gridx.isel(time=0,drop=True) 
		dsx = xr.merge([thickness, gridx])
		dsx.coords['xt_ocean'].attrs.update(axis='X')
		dsx.coords['xu_ocean'].attrs.update(axis='X', c_grid_axis_shift=0.5)
		dsx.coords['yt_ocean'].attrs.update(axis='Y')
		dsx.coords['yu_ocean'].attrs.update(axis='Y', c_grid_axis_shift=0.5)
		gridx = xgcm.Grid(dsx, periodic=['X'])
		
		
		#y
		gridy = xr.open_mfdataset('/g/data/x77/wf4500/ASC_project/GH/GH_AP2025/ASC_and_heat_transport/Python/OM2_grid/ocean_grid.nc', combine='by_coords')
		gridy=gridy.isel(time=0,drop=True)         
		dsy = xr.merge([thickness, gridy])
		dsy.coords['xt_ocean'].attrs.update(axis='X')
		dsy.coords['xu_ocean'].attrs.update(axis='X', c_grid_axis_shift=0.5)
		dsy.coords['yt_ocean'].attrs.update(axis='Y')
		dsy.coords['yu_ocean'].attrs.update(axis='Y', c_grid_axis_shift=0.5)
		gridy = xgcm.Grid(dsy, periodic=['X'])
		
		# Calculating the slope on the x and y grid
		thickness_xg=gridx.interp(dsx.st_edges_ocean, 'Y').isel(yt_ocean=slice(0,cutout_latind))
		thickness_yg=gridx.interp(dsy.st_edges_ocean, 'X').isel(yu_ocean=slice(0,cutout_latind))
		        
		       
		############################ Loading thicknesses     
		thickness_yg0=thickness_yg.load()
		thickness_xg0=thickness_xg.load()       
		        
		# Calculating the slope on the x and y grid
		thickness_xg=thickness_xg.rename(yt_ocean='y_ocean',xu_ocean='x_ocean')
		thickness_yg=thickness_yg.rename(yu_ocean='y_ocean',xt_ocean='x_ocean')
		
		# stack transports into 1d and drop any points not on contour:
		x_thick_1d = thickness_xg.stack(contour_index = ['y_ocean', 'x_ocean'])
		x_thick_1d = x_thick_1d.where(mask_x_numbered_1d>0, drop = True)
		y_thick_1d = thickness_yg.stack(contour_index = ['y_ocean', 'x_ocean'])
		y_thick_1d = y_thick_1d.where(mask_y_numbered_1d>0, drop = True)
		
		# combine all points on contour:
		Uthick = xr.concat((x_thick_1d, y_thick_1d), dim = 'contour_index')
		Uthick = Uthick.sortby(contour_ordering)
		Uthick.coords['contour_index'] = contour_index_array
		Uthick = Uthick.load()
		
		
		distance_along_contour=np.ones(num_points)
		for n in range(num_points-1):
		    coords_1 = (lat_along_contour[n], lon_along_contour[n])
		    coords_2 = (lat_along_contour[n+1], lon_along_contour[n+1])
		    distance_along_contour[n]=geopy.distance.geodesic(coords_1, coords_2).km
		
		    
		#for the last casE CORNER
		coords_1 = (lat_along_contour[-1], lon_along_contour[-1])
		coords_2 = (lat_along_contour[0], lon_along_contour[0])
		distance_along_contour[num_points-1]=geopy.distance.geodesic(coords_1, coords_2).km
		
		distance_along_contour=np.cumsum(distance_along_contour)
		        
         
		        
		# PREVIOUS SAVEDIR TO LOCATE THE FILES save_dir  = '/g/data/x77/wf4500/ASC_project/ASC_speed/monthly/OM2_RYF/'
		save_dir  = '/g/data/x77/wf4500/ASC_project/ASC_speed/daily_z/OM2_IAF_XYgrid_650m/'
		file_name = 'Antarctic_slope_contour_650m_velocities_'
		print('saving file on... ' + save_dir + file_name + start_time +"_uv.nc")	
		
		data_u=xr.DataArray((Ualong),name="u_along_contour",dims=["time","st_ocean","contour_index"])
		data_v=xr.DataArray((Ucross),name="u_cross_contour",dims=["time","st_ocean","contour_index"])
		data_lon=xr.DataArray((lon_along_contour),name="lon_along_contour",dims=["contour_index"])
		data_lat=xr.DataArray((lat_along_contour),name="lat_along_contour",dims=["contour_index"])
		data_dist=xr.DataArray((distance_along_contour),name="distance_along_contour",dims=["contour_index"])


		data_thick=xr.DataArray((Uthick),name="thickness_contour",dims=["st_ocean","contour_index"])
		data_time=xr.DataArray(Ualong.time,name="time",dims=["time"])
		
		data_to_saveuv = xr.merge([data_u,data_v,data_lon,data_lat,data_thick,data_dist])
		data_to_saveuv['contour_index']=Ualong.contour_index
		data_to_saveuv['st_ocean']=np.array(st_ocean)
		data_to_saveuv['time']=Ualong.time
		data_to_saveuv
		data_to_saveuv.to_netcdf(save_dir + file_name + start_time +"_uv.nc")



