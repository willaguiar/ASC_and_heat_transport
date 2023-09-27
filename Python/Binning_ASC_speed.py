# # Compute daily ASC along contour, on x,y grid and save in netcdfs for each month. U,V converted from uhrho and vhrho



import cosima_cookbook as cc
import matplotlib.pyplot as plt
import netCDF4 as nc
import xarray as xr
import numpy as np
import sys,os
from dask.distributed import Client
import glob
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
	
	#importing the speed along the contour
	import glob
	src = '/g/data/x77/wf4500/ASC_project/ASC_speed/daily_z/OM2_IAF_XYgrid/Antarctic_slope_contour_1km_velocities_*'
	gl =glob.glob(src + "****_uv.nc")
	gl=sorted(gl); gl
	print("%i netcdf data files"%len(gl))
	##importing data from ocean.nc
	Along_slope_v=xr.open_mfdataset(gl, concat_dim='time', combine='nested').mean('time')
	Along_slope_v_time=xr.open_mfdataset(gl, concat_dim='time', combine='nested')
	
	
	#yearinit=2177
	yearinit = sys.argv[2]
	imonth = int(sys.argv[1])
	month = str(int(sys.argv[1])) 
	month = month.zfill(2)
	monthseries=['01','02','03','04','05','06','07','08','09','10','11','12']
	dayseries=['31','28','31','30','31','30','31','31','30','31','30','31'] 
	mi=0
	#Accounting for Leap years
	long_feby=np.arange(1960,2028,4)
	if yearinit in long_feby: dayseries=['31','29','31','30','31','30','31','31','30','31','30','31'] ; print('this is a leap year')
	for n in np.arange(0,1): 
        
		start_time=str(yearinit) + '-' + month + '-01'
		end_time=str(yearinit) + '-' + month + '-' + dayseries[imonth-1]
		start_time00=str(yearinit) + '_' + month
		
		Along_slope_v_time = Along_slope_v_time.sel(time=slice(start_time,end_time))
		ndays=np.size(Along_slope_v_time.time)        

		lon_along_contour=np.array(Along_slope_v_time.lon_along_contour.isel(time=0))
		lat_along_contour=np.array(Along_slope_v_time.lat_along_contour.isel(time=0))
		Along_slope_v_time=Along_slope_v_time.load()

		# convert to longitude coordinate and average into 3 degree longitude bins:
		bin_width = 3
		bin_spacing = 0.25
		lon_west = -280
		lon_east = 80
        
		surf_U=Along_slope_v_time.u_along_contour.isel(st_ocean=0,time=0)
        
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
        
		# find average latitude of each bin, so we can plot back on the isobath:
		lat_bin_midpoints = np.zeros(n_bins)
		for lon_bin in range(n_bins):
			# find nearest isobath point:
			lon_index = np.where(lon_along_contour>=lon_bin_midpoints[lon_bin])[0][0]
			lat_bin_midpoints[lon_bin] = lat_along_contour[lon_index]
        


		from joblib import Parallel, delayed
		def binning_Ualong(Along_slope_v_time,lon_along_contour,lat_along_contour, nn):
        	#nn will be ii*jj, just to make the paralelizarion more efficient    
        	#defining ii and tt
			iilist=np.tile(np.arange(75),31)
			ttlist=np.sort(np.tile(np.arange(31),75))  
        	# convert to longitude coordinate and average into 3 degree longitude bins:
        	# in degrees:
			bin_width = 3
			bin_spacing = 0.25
			lon_west = -280
			lon_east = 80
        	
			ii=iilist[nn];tt=ttlist[nn]
			surf_U=Along_slope_v_time.u_along_contour.isel(st_ocean=ii,time=tt)
        	
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
			surf_U_binned = np.zeros(n_bins)
			for lon_bin in range(n_bins):
				surf_U_this_bin = surf_U[bin_edge_indices_west[lon_bin]:bin_edge_indices_east[lon_bin]]
				surf_U_binned[lon_bin] = np.mean(surf_U_this_bin)
        
        	        
			return surf_U_binned
        
        
		surf_U_binned = Parallel(n_jobs=-1)(delayed(binning_Ualong)(Along_slope_v_time,lon_along_contour,lat_along_contour, nn) for nn in range(75*ndays))

		surf_U_binned_3d=np.reshape(surf_U_binned,[ndays,75,1428]);



		# PREVIOUS SAVEDIR TO LOCATE THE FILES save_dir  = '/g/data/x77/wf4500/ASC_project/ASC_speed/monthly/OM2_RYF/'
		#save_dir  = '/g/data/x77/wf4500/ASC_project/ASC_speed/daily_z/OM2_IAF_XYgrid/'  
		save_dir='/g/data/nm03/wf4500/ASC_project/model_data/access-om2/binned_ASC_speed/daily/IAF/'
		file_name = 'Binned_Antarctic_slope_contour_1km_velocities'
		print('saving file on... ' + save_dir + file_name + start_time00 +".nc")	
		
		data_along=xr.DataArray((surf_U_binned_3d),name="u_along_contour_binned",dims=["time","st_ocean","lon_bin_midpoints"])
		data_lat=xr.DataArray((lat_bin_midpoints),name="lat_bin_midpoints",dims=["lon_bin_midpoints"])



		
		data_to_saveuv = xr.merge([data_along,data_lat])
		data_to_saveuv['st_ocean']=Along_slope_v_time.st_ocean
		data_to_saveuv['time']=Along_slope_v_time.time
		data_to_saveuv['lon_bin_midpoints']=lon_bin_midpoints
		data_to_saveuv
		data_to_saveuv.to_netcdf(save_dir + file_name + start_time00 +".nc")



