
import xarray as xr
import pyinterp
import numpy as np
import netCDF4
import logging 



def write_stat(nc, group_name, binning):
    
    grp = nc.createGroup(group_name)
    grp.createDimension('lon', len(binning.x))
    grp.createDimension('lat', len(binning.y))
    
    longitude = grp.createVariable('lon', 'f4', 'lon', zlib=True)
    longitude[:] = binning.x
    latitude = grp.createVariable('lat', 'f4', 'lat', zlib=True)
    latitude[:] = binning.y
    
    stats = ['min', 'max', 'sum', 'sum_of_weights', 'variance', 'mean', 'count', 'kurtosis', 'skewness']
    for variable in stats:
        
        var = grp.createVariable(variable, binning.variable(variable).dtype, ('lat','lon'), zlib=True)
        var[:, :] = binning.variable(variable).T 


def compute_stats_map(ds, 
                      output_filename,
                      lon_bins=np.arange(0, 360, 1), 
                      lat_bins=np.arange(-90, 91, 1)
                      ):
                      
    ncfile = netCDF4.Dataset(output_filename,'w')
    
    binning = pyinterp.Binning2D(pyinterp.Axis(lon_bins, is_circle=True), pyinterp.Axis(lat_bins))

    # binning alongtrack
    binning.push(ds.longitude, ds.latitude, ds.sla_unfiltered, simple=True)
    write_stat(ncfile, 'alongtrack', binning)
    binning.clear()

    # binning map interp
    binning.push(ds.longitude, ds.latitude, ds.msla_interpolated, simple=True)
    write_stat(ncfile, 'maps', binning)
    binning.clear()

    # binning diff sla-msla
    binning.push(ds.longitude, ds.latitude, ds.sla_unfiltered - ds.msla_interpolated, simple=True)
    write_stat(ncfile, 'diff', binning)
    binning.clear()

    # add rmse
    diff2 = (ds.sla_unfiltered - ds.msla_interpolated)**2
    binning.push(ds.longitude, ds.latitude, diff2, simple=True)
    var = ncfile.groups['diff'].createVariable('rmse', binning.variable('mean').dtype, ('lat','lon'), zlib=True)
    var[:, :] = np.sqrt(binning.variable('mean')).T 
    
    # add rms alongtrack
    diff2 = (ds.sla_unfiltered)**2
    binning.push(ds.longitude, ds.latitude, diff2, simple=True)
    var = ncfile.groups['alongtrack'].createVariable('rms', binning.variable('mean').dtype, ('lat','lon'), zlib=True)
    var[:, :] = np.sqrt(binning.variable('mean')).T 
    
    # add rms map interp
    diff2 = (ds.msla_interpolated)**2
    binning.push(ds.longitude, ds.latitude, diff2, simple=True)
    var = ncfile.groups['maps'].createVariable('rms', binning.variable('mean').dtype, ('lat','lon'), zlib=True)
    var[:, :] = np.sqrt(binning.variable('mean')).T 
    
    ncfile.close()
    
    logging.info(f'  Results saved in: {output_filename}')
    
    
def compute_stats_timeseries(ds_interp, 
                             output_filename, 
                             freq='1D'):
    
    lat_min = -60
    lat_max = 60
    
    diff = ds_interp['sla_unfiltered'] - ds_interp['msla_interpolated']
    
    # resample 
    da_resample = diff.where((ds_interp.latitude >= lat_min) & (ds_interp.latitude<=lat_max), drop=True).resample(time=freq)
    
    # compute stats
    vmean = da_resample.mean()
    vminimum = da_resample.min()
    vmaximum = da_resample.max()
    vcount = da_resample.count()
    vvariance = da_resample.var()
    #vmedian = da_resample.median()
    vrms = np.sqrt(np.square(diff.where((ds_interp.latitude >= lat_min) & (ds_interp.latitude<=lat_max), drop=True)).resample(time=freq).mean())
    
    rmse = np.copy(vrms)
    
    # save stat to dataset
    ds = xr.Dataset(
        {
            "mean": (("time"), vmean.values),
            "min": (("time"), vminimum.values),
            "max": (("time"), vmaximum.values),
            "count": (("time"), vcount.values),
            "variance": (("time"), vvariance.values),
            #"median": (("time"), vmedian.values),
            "rms": (("time"), vrms.values),            
        },
        {"time": vmean['time']},
    )
    
    ds.to_netcdf(output_filename, group='diff')
    
    
    # resample 
    da_resample = ds_interp['sla_unfiltered'].where((ds_interp.latitude >= lat_min) & (ds_interp.latitude <= lat_max), drop=True).resample(time=freq)
    
    # compute stats
    vmean = da_resample.mean()
    vminimum = da_resample.min()
    vmaximum = da_resample.max()
    vcount = da_resample.count()
    vvariance = da_resample.var()
    #vmedian = da_resample.median()
    vrms = np.sqrt(np.square(ds_interp['sla_unfiltered'].where((ds_interp.latitude >= lat_min) & (ds_interp.latitude <= lat_max), drop=True)).resample(time=freq).mean())
    
    rms_alongtrack = np.copy(vrms)
    
    # save stat to dataset
    ds = xr.Dataset(
        {
            "mean": (("time"), vmean.values),
            "min": (("time"), vminimum.values),
            "max": (("time"), vmaximum.values),
            "count": (("time"), vcount.values),
            "variance": (("time"), vvariance.values),
            #"median": (("time"), vmedian.values),
            "rms": (("time"), vrms.values),            
        },
        {"time": vmean['time']},
    )
    
    ds.to_netcdf(output_filename, group='alongtrack', mode='a')

    
    # resample 
    da_resample = ds_interp['msla_interpolated'].where((ds_interp.latitude >= lat_min) & (ds_interp.latitude <= lat_max), drop=True).resample(time=freq)
    
    # compute stats
    vmean = da_resample.mean()
    vminimum = da_resample.min()
    vmaximum = da_resample.max()
    vcount = da_resample.count()
    vvariance = da_resample.var()
    #vmedian = da_resample.median()
    vrms = np.sqrt(np.square(ds_interp['msla_interpolated'].where((ds_interp.latitude >= lat_min) & (ds_interp.latitude <= lat_max), drop=True)).resample(time=freq).mean())
    
    # save stat to dataset
    ds = xr.Dataset(
        {
            "mean": (("time"), vmean.values),
            "min": (("time"), vminimum.values),
            "max": (("time"), vmaximum.values),
            "count": (("time"), vcount.values),
            "variance": (("time"), vvariance.values),
            #"median": (("time"), vmedian.values),
            "rms": (("time"), vrms.values),            
        },
        {"time": vmean['time']},
    )
    
    ds.to_netcdf(output_filename, group='maps', mode='a')
    
    logging.info(' ')
    logging.info(f'  Results saved in: {output_filename}')
    
    rmse_score = 1. - rmse/rms_alongtrack
    # mask score if nb obs < nb_min_obs
    nb_min_obs = 10
    rmse_score = np.ma.masked_where(vcount.values < nb_min_obs, rmse_score)
    
    mean_rmse = np.ma.mean(np.ma.masked_invalid(rmse_score))
    std_rmse = np.ma.std(np.ma.masked_invalid(rmse_score))
    
    logging.info(' ')
    logging.info(f'  MEAN RMSE Score = {mean_rmse}')
    logging.info(' ')
    logging.info(f'  STD RMSE Score = {std_rmse}')
    
    return mean_rmse, std_rmse






def compute_current_stats_map(ds, 
                      output_filename,
                      lon_bins=np.arange(0, 360, 3), 
                      lat_bins=np.arange(-90, 93, 3)
                      ):
                      
    ncfile = netCDF4.Dataset(output_filename,'w')
    
    binning = pyinterp.Binning2D(pyinterp.Axis(lon_bins, is_circle=True), pyinterp.Axis(lat_bins))

    # binning drifter U
    binning.push(ds.longitude, ds.latitude, ds.EWCT, simple=True)
    write_stat(ncfile, 'u_drifter', binning)
    binning.clear()

    # binning U map interp
    binning.push(ds.longitude, ds.latitude, ds.ugos_interpolated, simple=True)
    write_stat(ncfile, 'u_maps', binning)
    binning.clear()

    # binning diff u_drifter-u_map
    binning.push(ds.longitude, ds.latitude, ds.EWCT - ds.ugos_interpolated, simple=True)
    write_stat(ncfile, 'u_diff', binning)
    binning.clear()

    # add rmse
    diff2 = (ds.EWCT - ds.ugos_interpolated)**2
    binning.push(ds.longitude, ds.latitude, diff2, simple=True)
    var = ncfile.groups['u_diff'].createVariable('rmse', binning.variable('mean').dtype, ('lat','lon'), zlib=True)
    var[:, :] = np.sqrt(binning.variable('mean')).T 
    
    # add rms Udrifter
    diff2 = (ds.EWCT)**2
    binning.push(ds.longitude, ds.latitude, diff2, simple=True)
    var = ncfile.groups['u_drifter'].createVariable('rms', binning.variable('mean').dtype, ('lat','lon'), zlib=True)
    var[:, :] = np.sqrt(binning.variable('mean')).T 
    
    # add rms Umap interp
    diff2 = (ds.ugos_interpolated)**2
    binning.push(ds.longitude, ds.latitude, diff2, simple=True)
    var = ncfile.groups['u_maps'].createVariable('rms', binning.variable('mean').dtype, ('lat','lon'), zlib=True)
    var[:, :] = np.sqrt(binning.variable('mean')).T 
    
    
    
    
    # binning drifter V
    binning.push(ds.longitude, ds.latitude, ds.NSCT, simple=True)
    write_stat(ncfile, 'v_drifter', binning)
    binning.clear()

    # binning V map interp
    binning.push(ds.longitude, ds.latitude, ds.vgos_interpolated, simple=True)
    write_stat(ncfile, 'v_maps', binning)
    binning.clear()

    # binning diff v_drifter-u_map
    binning.push(ds.longitude, ds.latitude, ds.NSCT - ds.vgos_interpolated, simple=True)
    write_stat(ncfile, 'v_diff', binning)
    binning.clear()

    # add rmse
    diff2 = (ds.NSCT - ds.vgos_interpolated)**2
    binning.push(ds.longitude, ds.latitude, diff2, simple=True)
    var = ncfile.groups['v_diff'].createVariable('rmse', binning.variable('mean').dtype, ('lat','lon'), zlib=True)
    var[:, :] = np.sqrt(binning.variable('mean')).T 
    
    # add rms Vdrifter
    diff2 = (ds.NSCT)**2
    binning.push(ds.longitude, ds.latitude, diff2, simple=True)
    var = ncfile.groups['v_drifter'].createVariable('rms', binning.variable('mean').dtype, ('lat','lon'), zlib=True)
    var[:, :] = np.sqrt(binning.variable('mean')).T 
    
    # add rms Vmap interp
    diff2 = (ds.vgos_interpolated)**2
    binning.push(ds.longitude, ds.latitude, diff2, simple=True)
    var = ncfile.groups['v_maps'].createVariable('rms', binning.variable('mean').dtype, ('lat','lon'), zlib=True)
    var[:, :] = np.sqrt(binning.variable('mean')).T 
    
    
    ncfile.close()
    
    logging.info(f'  Results saved in: {output_filename}')
