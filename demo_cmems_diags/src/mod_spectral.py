import scipy.signal
from netCDF4 import Dataset
import numpy as np
import pyinterp
import xarray as xr
import hvplot.xarray
import logging

def compute_median_dx(dataset):
        
    return 0.001*np.median(pyinterp.geodetic.coordinate_distances(dataset['longitude'][:-1].values,
                                                              dataset['latitude'][:-1].values,
                                                              dataset['longitude'][1:].values,
                                                              dataset['latitude'][1:].values
                                                             ))


def compute_median_lon_lat(vlon, vlat, sub_segment_point, npt):
    
    # Near Greenwhich case
    if ((vlon[sub_segment_point + npt - 1] < 50.)
        and (vlon[sub_segment_point] > 320.)) \
            or ((vlon[sub_segment_point + npt - 1] > 320.)
                and (vlon[sub_segment_point] < 50.)):
        
        # make lon negative for lon > 180.
        tmp_lon = np.where(vlon[sub_segment_point:sub_segment_point + npt] > 180,
                           vlon[sub_segment_point:sub_segment_point + npt] - 360,
                           vlon[sub_segment_point:sub_segment_point + npt])
        # compute median longitude
        median_lon_sub_segment = np.median(tmp_lon)%360.
    
    # Far from greenwich
    else:
        median_lon_sub_segment = np.median(vlon[sub_segment_point:sub_segment_point + npt])

    median_lat_sub_segment = np.median(vlat[sub_segment_point:sub_segment_point + npt])
    
    return median_lon_sub_segment, median_lat_sub_segment


def compute_segment(ds_interp, lenght_scale):
    
    ds_interp = ds_interp.sortby('time').dropna('time')
    
    lon_along_track = ds_interp['longitude'].values
    lat_along_track = ds_interp['latitude'].values
    sla = ds_interp['sla_unfiltered'].values
    msla = ds_interp['msla_interpolated'].values
    
    delta_x = compute_median_dx(ds_interp) # in km
    npt = int(lenght_scale / delta_x)
    
    segment_overlapping = 0.25
    
    indi = np.where((np.diff(ds_interp['time']) > np.timedelta64(2, 's')))[0]
    track_segment_lenght = np.insert(np.diff(indi), [0], indi[0])
    selected_track_segment = np.where(track_segment_lenght >= npt)[0]

    list_lat_segment = []
    list_lon_segment = []
    list_sla_segment = []
    list_msla_segment = []

    if selected_track_segment.size > 0:

        for track in selected_track_segment:

            if track-1 >= 0:
                index_start_selected_track = indi[track-1]
                index_end_selected_track = indi[track]
            else:
                index_start_selected_track = 0
                index_end_selected_track = indi[track]

            start_point = index_start_selected_track
            end_point = index_end_selected_track

            for sub_segment_point in range(start_point, end_point - npt, int(npt*segment_overlapping)):
            
                mean_lon_sub_segment, mean_lat_sub_segment = compute_median_lon_lat(lon_along_track, 
                                                                                    lat_along_track,
                                                                                    sub_segment_point,
                                                                                    npt)

                sla_segment = sla[sub_segment_point:sub_segment_point + npt]
                msla_segment = msla[sub_segment_point:sub_segment_point + npt]
            
                list_sla_segment.append(sla_segment)
                list_lon_segment.append(mean_lon_sub_segment)
                list_lat_segment.append(mean_lat_sub_segment)
                list_msla_segment.append(msla_segment)
                
    return np.asarray(list_lon_segment), np.asarray(list_lat_segment), np.asarray(list_sla_segment), np.asarray(list_msla_segment), delta_x, npt


def spectral_computation(lon_segment, lat_segment, ref_segments, study_segments, delta_x, npt):


    delta_lat = 5.
    delta_lon = 5.
    nb_min_segment = 10
    vlon = np.arange(0., 360., 1.)
    vlat = np.arange(-80., 91., 1)

    list_mean_psd_ref = []
    list_nb_segment = []
    list_mean_frequency = []
    list_mean_psd_study = []
    list_mean_psd_diff_study_ref = []
    list_mean_coherence = []
    list_mean_cross_spectrum = []

    fs = 1.0 / delta_x

    # Loop over output lon/lat boxes and selection of the segment within the box plus/minus delta_lon/lat
    for ilat in vlat:

        lat_min = ilat - 0.5*delta_lat
        lat_max = ilat + 0.5*delta_lat

        selected_lat_index = np.where(np.logical_and(lat_segment >= lat_min, lat_segment <= lat_max))[0]
        ref_segments_tmp = ref_segments[selected_lat_index]
        if study_segments is not None:
            study_segments_tmp = study_segments[selected_lat_index]
        else:
            study_segments_tmp = None

        for ilon in vlon:

            lon_min = ilon - 0.5*delta_lon
            lon_max = ilon + 0.5*delta_lon

            if (lon_min < 0.) and (lon_max > 0.):
                selected_segment = np.where(np.logical_or(lon_segment[selected_lat_index] % 360. >= lon_min + 360.,
                                                          lon_segment[selected_lat_index] % 360. <= lon_max))[0]
            elif (lon_min > 0.) and (lon_max > 360.):
                selected_segment = np.where(np.logical_or(lon_segment[selected_lat_index] % 360. >= lon_min,
                                                          lon_segment[selected_lat_index] % 360. <= lon_max - 360.))[0]
            else:
                selected_segment = np.where(np.logical_and(lon_segment[selected_lat_index] % 360. >= lon_min,
                                                           lon_segment[selected_lat_index] % 360. <= lon_max))[0]

            if len(selected_segment) > nb_min_segment:
                selected_ref_segments = np.ma.masked_where(ref_segments_tmp[selected_segment].flatten() > 1.E10,
                                                           ref_segments_tmp[selected_segment].flatten())

                # Power spectrum density reference field
                wavenumber, psd_ref = scipy.signal.welch(selected_ref_segments,
                                                         fs=fs,
                                                         nperseg=npt,
                                                         scaling='density',
                                                         noverlap=0)

                list_mean_frequency.append(wavenumber)
                list_mean_psd_ref.append(psd_ref)
                list_nb_segment.append(selected_segment.size)


                selected_study_segments = np.ma.masked_where(
                        study_segments_tmp[selected_segment].flatten() > 1.E10,
                        study_segments_tmp[selected_segment].flatten())

                # Compute diff study minus ref
                diff_study_ref = selected_study_segments - selected_ref_segments

                # Power spectrum density of the error between to field
                wavenumber, psd_diff_study_ref = scipy.signal.welch(diff_study_ref,
                                                                        fs=fs,
                                                                        nperseg=npt,
                                                                        scaling='density',
                                                                        noverlap=0)

                # Power spectrum density study field
                wavenumber, psd_study = scipy.signal.welch(selected_study_segments,
                                                               fs=fs,
                                                               nperseg=npt,
                                                               scaling='density',
                                                               noverlap=0)

                # Magnitude square coherence between the ref and study field
                wavenumber, coherence = scipy.signal.coherence(selected_study_segments,
                                                                   selected_ref_segments,
                                                                   fs=fs,
                                                                   nperseg=npt,
                                                                   noverlap=0)

                # Cross spectrum
                wavenumber, cross_spectrum = scipy.signal.csd(selected_study_segments,
                                                                  selected_ref_segments,
                                                                  fs=fs,
                                                                  nperseg=npt,
                                                                  noverlap=0)

                list_mean_psd_study.append(psd_study)
                list_mean_psd_diff_study_ref.append(psd_diff_study_ref)
                list_mean_coherence.append(coherence)
                list_mean_cross_spectrum.append(cross_spectrum)

            else:

                list_mean_frequency.append(np.zeros((int(npt / 2) + 1)))
                list_mean_psd_ref.append(np.zeros((int(npt / 2) + 1)))
                list_nb_segment.append(0.)

                list_mean_psd_study.append(np.zeros((int(npt / 2) + 1)))
                list_mean_psd_diff_study_ref.append(np.zeros((int(npt / 2) + 1)))
                list_mean_coherence.append(np.zeros((int(npt / 2) + 1)))
                list_mean_cross_spectrum.append(np.zeros((int(npt / 2) + 1)))

    wavenumber = np.ma.mean(np.ma.masked_invalid(np.ma.masked_where(np.asarray(list_mean_frequency) == 0,
                                                              np.asarray(list_mean_frequency))), axis=0).filled(0.)
    
    nb_segment = np.asarray(list_nb_segment).reshape((vlat.size, vlon.size))
    psd_ref = np.transpose(np.asarray(list_mean_psd_ref)).reshape((wavenumber.size, vlat.size, vlon.size))
    psd_study = np.transpose(np.asarray(list_mean_psd_study)).reshape((wavenumber.size, vlat.size, vlon.size))
    psd_diff = np.transpose(np.asarray(list_mean_psd_diff_study_ref)).reshape((wavenumber.size, vlat.size, vlon.size))
    coherence = np.transpose(np.asarray(list_mean_coherence)).reshape((wavenumber.size, vlat.size, vlon.size))
    cross_spectrum = np.transpose(np.asarray(list_mean_cross_spectrum)).reshape((wavenumber.size, vlat.size, vlon.size))
    
    return wavenumber, vlat, vlon, nb_segment, psd_ref, psd_study, psd_diff, coherence, cross_spectrum



def compute_crossing(array, wavenumber, threshold=0.5):
    """
    :param array:
    :param wavenumber:
    :param threshold:
    :return:
    """


    flag_multiple_crossing = False
    zero_crossings = np.where(np.diff(np.sign(array - threshold)))[0]
    
    if len(zero_crossings) > 1:
        #print('Multiple crossing', len(zero_crossings))
        flag_multiple_crossing = True
        # MB add for large scale bais
        zero_crossings[0] = zero_crossings[-1]
         
    if len(zero_crossings) > 0:
        if zero_crossings[0] + 1 < array.size:

            array1 = array[zero_crossings[0]] - threshold
            array2 = array[zero_crossings[0] + 1] - threshold
            dist1 = np.log(wavenumber[zero_crossings[0]])
            dist2 = np.log(wavenumber[zero_crossings[0] + 1])
            log_wavenumber_crossing = dist1 - array1 * (dist1 - dist2) / (array1 - array2)
            resolution_scale = 1. / np.exp(log_wavenumber_crossing)

        else:
            resolution_scale = 0.

    else:
        resolution_scale = 0.

    return resolution_scale, flag_multiple_crossing


def compute_resolution(lon, lat, wavenumber, psd_diff, psd_ref):
    
    ratio = psd_diff/psd_ref

    resolution = np.empty((lat.size, lon.size))
    for jj in range(lat.size):
        for ii in range(lon.size):
            if not np.ma.is_masked(ratio[:, jj, ii]):
                resolution[jj, ii], flag = compute_crossing(ratio[:, jj, ii], wavenumber)
    
    return resolution


def write_psd_output(output_netcdf_file, wavenumber, vlat, vlon, nb_segment, psd_ref, psd_study, psd_diff_ref_study, coherence, cross_spectrum):

    nc_out = Dataset(output_netcdf_file, 'w', format='NETCDF4')
    
    nc_out.createDimension('wavenumber', wavenumber.size)
    nc_out.createDimension('lat', vlat.size)
    nc_out.createDimension('lon', vlon.size)

    wavenumber_out = nc_out.createVariable('wavenumber', 'f8', 'wavenumber', zlib=True)
    wavenumber_out.units = "1/km"
    wavenumber_out.axis = 'T'
    wavenumber_out[:] = wavenumber

    nb_segment_out = nc_out.createVariable('nb_segment', 'f8', ('lat', 'lon'), zlib=True)
    nb_segment_out.long_name = "number of segment used in spectral computation"
    nb_segment_out[:, :] = np.ma.masked_where(nb_segment == 0., nb_segment).filled(np.nan)

    lat_out = nc_out.createVariable('lat', 'f4', 'lat', zlib=True)
    lat_out[:] = vlat
    lat_out.units = 'degree_north'
    lon_out = nc_out.createVariable('lon', 'f4', 'lon', zlib=True)
    lon_out[:] = vlon
    lon_out.units = 'degree_east'

    psd_ref_out = nc_out.createVariable('psd_ref', 'f8', ('wavenumber', 'lat', 'lon'), zlib=True)
    psd_ref_out.units = 'm2/km'
    psd_ref_out.coordinates = "wavenumber lat lon"
    psd_ref_out.long_name = "power spectrum density reference field"
    psd_ref_out[:, :, :] = np.ma.masked_invalid(np.ma.masked_where(psd_ref == 0, psd_ref)).filled(np.nan)
    
    psd_study_out = nc_out.createVariable('psd_study', 'f8', ('wavenumber', 'lat', 'lon'), zlib=True)
    psd_study_out.units = 'm2/km'
    psd_study_out.coordinates = "wavenumber lat lon"
    psd_study_out.long_name = "power spectrum density study field"
    psd_study_out[:, :, :] = np.ma.masked_invalid(np.ma.masked_where(psd_study == 0, psd_study)).filled(np.nan)
    
    psd_diff = nc_out.createVariable('psd_diff', 'f8', ('wavenumber', 'lat', 'lon'), zlib=True)
    psd_diff.units = 'm2/km'
    psd_diff.coordinates = "wavenumber lat lon"
    psd_diff.long_name = "power spectrum density of difference study minus reference field"
    psd_diff[:, :, :] = np.ma.masked_invalid(np.ma.masked_where(psd_diff_ref_study == 0, psd_diff_ref_study))
    
    resolution = compute_resolution(vlon, vlat, wavenumber, psd_diff_ref_study, psd_ref)
    resolution_out = nc_out.createVariable('effective_resolution', 'f8', ('lat', 'lon'), zlib=True)
    resolution_out.units = 'km'
    resolution_out.coordinates = "lat lon"
    resolution_out.long_name = "Effective resolution computed as wavelenght where psd_err/psd_ref=0.5"
    resolution_out[:, :] = np.ma.masked_invalid(np.ma.masked_where(resolution == 0, resolution)).filled(np.nan)

    coherence_out = nc_out.createVariable('coherence', 'f8', ('wavenumber', 'lat', 'lon'), zlib=True)
    coherence_out[:, :, :] = np.ma.masked_invalid(np.ma.masked_where(coherence == 0, coherence)).filled(np.nan)
    coherence_out.coordinates = "wavenumber lat lon"
    coherence_out.long_name = "magnitude squared coherence between reference and study fields"

    cross_spectrum_real_out = nc_out.createVariable('cross_spectrum_real', 'f8', ('wavenumber', 'lat', 'lon'),
                                                        zlib=True)
    cross_spectrum_real_out[:, :, :] = np.ma.masked_invalid(np.ma.masked_where(np.real(cross_spectrum) == 0., np.real(cross_spectrum))).filled(np.nan)
    cross_spectrum_real_out.coordinates = "wavenumber lat lon"
    cross_spectrum_real_out.long_name = "real part of cross_spectrum between reference and study fields"
    cross_spectrum_imag_out = nc_out.createVariable('cross_spectrum_imag', 'f8', ('wavenumber', 'lat', 'lon'),
                                                        zlib=True)
    cross_spectrum_imag_out[:, :, :] = np.ma.masked_invalid(np.ma.masked_where(np.imag(cross_spectrum) == 0., np.imag(cross_spectrum))).filled(np.nan)
    cross_spectrum_imag_out.coordinates = "wavenumber lat lon"
    cross_spectrum_imag_out.long_name = "imaginary part of cross_spectrum between reference and study fields"

    nc_out.close()
    
    
def compute_psd_scores(ds_interp, output_filename, lenght_scale=1500. ):
    
    logging.info('Segment computation...')
    lon_segment, lat_segment, sla_segment, msla_segment, delta_x, npt = compute_segment(ds_interp, lenght_scale)
    
    logging.info('Spectral analysis...')
    wavenumber, vlat, vlon, nb_segment, psd_ref, psd_study, psd_diff, coherence, cross_spectrum = spectral_computation(lon_segment,
                                                                                                                       lat_segment,
                                                                                                                       sla_segment,
                                                                                                                       msla_segment, 
                                                                                                                       delta_x,
                                                                                                                       npt
                                                                                                                      )
    logging.info('Saving ouput...')
    write_psd_output(output_filename, wavenumber, vlat, vlon, nb_segment, psd_ref, psd_study, psd_diff, coherence, cross_spectrum)
    
    

def plot_psd_scores(filename):
    
    ds_psd = xr.open_dataset(filename)
    ds_psd['wavelenght'] = 1./ds_psd['wavenumber']
    ds_psd = ds_psd.assign_coords(wavelenght=ds_psd['wavelenght'])
    
    fig1 = ds_psd.psd_ref.hvplot.line(x='wavelenght', logx=True, logy=True, label='PSD ref', color='k', flip_xaxis=True, line_width=3)*\
    ds_psd.psd_study.hvplot.line(x='wavelenght', logx=True, logy=True, label='PSD study', color='r', flip_xaxis=True)*\
    ds_psd.psd_diff.hvplot.line(x='wavelenght', logx=True, logy=True, label='PSD err', color='grey', flip_xaxis=True)
    
    fig2 = ds_psd.coherence.hvplot.line(x='wavelenght', logx=True, logy=False, label='Coherence', color='k', flip_xaxis=True, line_width=2)*\
    (ds_psd.psd_diff/ds_psd.psd_ref).hvplot.line(x='wavelenght', logx=True, logy=False, label='PSD err/ PSD ref', color='r', flip_xaxis=True, line_width=2)
    
    return (ds_psd.effective_resolution.hvplot.quadmesh(x='lon', y='lat', clim=(100, 500), cmap='jet') +\
ds_psd.nb_segment.hvplot.quadmesh(x='lon', y='lat', cmap='jet') +\
fig1 + fig2).cols(2) 
            