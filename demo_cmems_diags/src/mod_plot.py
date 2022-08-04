import hvplot.xarray
from matplotlib import pyplot as plt
import numpy as np
import xarray as xr
import cartopy.crs as ccrs

def plot_maps(dataset, title=''):
    
    return dataset.hvplot.quadmesh(x='longitude', y='latitude', title=title, clim=(-0.8, 0.8), cmap='RdYlBu_r', rasterize=True)


def display_identification(ace, ce, lon_min=-180, lon_max=180, lat_min=-80, lat_max=90):
    
    # https://py-eddy-tracker.readthedocs.io/en/latest/python_module/02_eddy_identification/pet_display_id.html#sphx-glr-python-module-02-eddy-identification-pet-display-id-py
    
    fig = plt.figure(figsize=(15, 8))
    ax = fig.add_axes([0.03, 0.03, 0.90, 0.94])
    ax.set_aspect("equal")
    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)
    kwargs = dict(extern_only=True, color="k", lw=1)
    ace.display(ax, **kwargs), ce.display(ax, **kwargs)
    ace.filled(ax, "amplitude", cmap="magma_r", vmin=0, vmax=0.5)
    m = ce.filled(ax, "amplitude", cmap="magma_r", vmin=0, vmax=0.5)
    colorbar = plt.colorbar(m, cax=ax.figure.add_axes([0.95, 0.03, 0.02, 0.94]))
    colorbar.set_label("Amplitude (m)")
    
    plt.show()
    
    
def compare_identification_v0(ace1, ce1, title1, ace2, ce2, title2, lon_min=-180, lon_max=180, lat_min=-80, lat_max=90, ds_sst=None):
    
    fig = plt.figure(figsize=(15, 8))
    ax = fig.add_axes((0.05, 0.05, 0.9, 0.9))
    ax.set_aspect("equal")
    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)
    if ds_sst is not None:
        ax.pcolormesh(ds_sst['lat'], ds_sst['lat'], ds_sst['analysed_sst'])
    ace1.display(ax, label=f"Anticyclonic contour ({title1})", color="k", lw=2)
    ace2.display(ax, label=f"Anticyclonic contour ({title2})", color="lime", lw=2)
    ce1.display(ax, label=f"Cyclonic contour ({title1})", color="aqua", lw=2)
    ce2.display(ax, label=f"Cyclonic contour ({title2})", color="magenta", lw=2)
    _ = ax.legend(loc="upper right")
    

def compare_identification(ace1, ce1, title1, ace2, ce2, title2, lon_min=-180, lon_max=180, lat_min=-80, lat_max=90, ds_sst=None, ds_chloro=None):
    
    fig = plt.figure(figsize=(15, 16))
    #ax = fig.add_axes((0.05, 0.05, 0.9, 0.9))
    ax = plt.subplot(211)
    ax.set_aspect("equal")
    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)
    if ds_sst is not None:
        ds_sel = ds_sst.sel(lon=slice(lon_min, lon_max), lat=slice(lat_min, lat_max))
        vmin = np.nanpercentile(ds_sel['analysed_sst'][0, :, :], 5.)
        vmax = np.nanpercentile(ds_sel['analysed_sst'][0, :, :], 95.)
        c = ax.contourf(ds_sel['lon'], ds_sel['lat'], ds_sel['analysed_sst'][0, :, :], cmap='Spectral_r', levels=np.linspace(vmin, vmax, 100), extend='both', antialiased=True)
        plt.colorbar(c, label='SST Ostia [deg. K]')
    ace1.display(ax, label=f"Anticyclonic contour ({title1})", color="k", lw=1.5)
    ace2.display(ax, label=f"Anticyclonic contour ({title2})", color="lime", lw=1.5)
    ce1.display(ax, label=f"Cyclonic contour ({title1})", color="aqua", lw=1.5)
    ce2.display(ax, label=f"Cyclonic contour ({title2})", color="magenta", lw=1.5)
    _ = ax.legend(loc="upper right")
    
    
    ax = plt.subplot(212)
    ax.set_aspect("equal")
    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)
    if ds_chloro is not None:
        ds_sel = ds_chloro.where((ds_chloro.lon>=lon_min) & (ds_chloro.lon<=lon_max),drop=True)
        ds_sel = ds_sel.where((ds_sel.lat>=lat_min) & (ds_sel.lat<=lat_max),drop=True)
        vmin = np.nanpercentile(ds_sel['CHL'], 5.)
        vmax = np.nanpercentile(ds_sel['CHL'], 95.)
        c = ax.contourf(ds_sel['lon'], ds_sel['lat'], ds_sel['CHL'][:, :], cmap='viridis', levels=np.linspace(vmin, vmax, 100), extend='both', antialiased=True)
        plt.colorbar(c, label=' Chloro [milligram m-3]')
    ace1.display(ax, label=f"Anticyclonic contour ({title1})", color="k", lw=1.5)
    ace2.display(ax, label=f"Anticyclonic contour ({title2})", color="lime", lw=1.5)
    ce1.display(ax, label=f"Cyclonic contour ({title1})", color="aqua", lw=1.5)
    ce2.display(ax, label=f"Cyclonic contour ({title2})", color="magenta", lw=1.5)
    _ = ax.legend(loc="upper right")
    
    
        
def compare_identification_stats(ac1, c1, title1, ac2, c2, title2):
    fig = plt.figure(figsize=(12, 7))
    
    kwargs_ac1 = dict(label=f"Anticyclonic {title1}", color="k", histtype="step", density=True)
    kwargs_c1 = dict(label=f"Cyclonic {title1}", color="aqua", histtype="step", density=True)
    
    kwargs_ac2 = dict(label=f"Anticyclonic {title2}", color="lime", histtype="step", density=True)
    kwargs_c2 = dict(label=f"Cyclonic {title2}", color="magenta", histtype="step", density=True)

    for x0, name, title, xmax, factor, bins in zip(
        (0.4, 0.72, 0.08),
        ("speed_radius", "speed_average", "amplitude"),
        ("Speed radius (km)", "Speed average (cm/s)", "Amplitude (cm)"),
        (100, 50, 20),
        (0.001, 100, 100),
        (np.arange(0, 2000, 1), np.arange(0, 1000, 0.5), np.arange(0.0005, 1000, 0.2)),
    ):
        ax_hist = fig.add_axes((x0, 0.24, 0.27, 0.35))
        nb_ac1, _, _ = ax_hist.hist(ac1[name] * factor, bins=bins, **kwargs_ac1)
        nb_c1, _, _ = ax_hist.hist(c1[name] * factor, bins=bins, **kwargs_c1)
        nb_ac2, _, _ = ax_hist.hist(ac2[name] * factor, bins=bins, **kwargs_ac2)
        nb_c2, _, _ = ax_hist.hist(c2[name] * factor, bins=bins, **kwargs_c2)
        
        ax_hist.set_xticklabels([])
        ax_hist.set_xlim(0, xmax)
        ax_hist.grid()

        ax_cum = fig.add_axes((x0, 0.62, 0.27, 0.35))
        ax_cum.hist(ac1[name] * factor, bins=bins, cumulative=-1, **kwargs_ac1)
        ax_cum.hist(c1[name] * factor, bins=bins, cumulative=-1, **kwargs_c1)
        ax_cum.hist(ac2[name] * factor, bins=bins, cumulative=-1, **kwargs_ac2)
        ax_cum.hist(c2[name] * factor, bins=bins, cumulative=-1, **kwargs_c2)
        ax_cum.set_xticklabels([])
        ax_cum.set_title(title)
        ax_cum.set_xlim(0, xmax)
        ax_cum.set_ylim(0, 1)
        ax_cum.grid()

        ax_ratio = fig.add_axes((x0, 0.06, 0.27, 0.15))
        ax_ratio.set_xlim(0, xmax)
        ax_ratio.set_ylim(0, 2)
        ax_ratio.plot((bins[1:] + bins[:-1]) / 2, nb_c1 / nb_ac1)
        ax_ratio.plot((bins[1:] + bins[:-1]) / 2, nb_c2 / nb_ac2)
        ax_ratio.axhline(1, color="k")
        ax_ratio.grid()
        ax_ratio.set_xlabel(title)

    ax_cum.set_ylabel("Cumulative\npercent distribution")
    ax_hist.set_ylabel("Percent of observations")
    ax_ratio.set_ylabel("Ratio percent\nCyc/Acyc")
    ax_cum.legend()
    
    plt.show()
    
    
    
### PLOT SCORE
def plot_map_scores(filename):
    
    ds_score = xr.open_dataset(filename, group='diff')
    rmse_map = ds_score['rmse']
    error_variance_map = ds_score['variance']
    std_alongtrack = np.sqrt(error_variance_map)
    
    ds_score = xr.open_dataset(filename, group='alongtrack')
    rms_alongtrack = ds_score['rms']
    variance_alongtrack = ds_score['variance']
    mean_alongtrack = ds_score['mean']
    
    ds_score = xr.open_dataset(filename, group='maps')
    mean_map = ds_score['mean']
    
    
    
    rmse_score = 1. - rmse_map/rms_alongtrack
    rmse_score.name = 'rmse_score'
    
    variance_score = 1. - error_variance_map/variance_alongtrack
    variance_score.name = 'variance_score'
    
    bias = mean_map - mean_alongtrack
    bias.name = 'bias'
    vmin_bias = np.nanpercentile(bias, 5)
    vmax_bias = np.nanpercentile(bias, 95)
    
    normalized_bias = (mean_map - mean_alongtrack)/std_alongtrack
    normalized_bias.name = 'normalized_bias'
    vmin_normalized_bias = np.nanpercentile(normalized_bias, 5)
    vmax_normalized_bias = np.nanpercentile(normalized_bias, 95)
    
    return (rmse_score.hvplot.quadmesh(x='lon', y='lat', clim=(0, 1), cmap='RdYlGn', title='RMSE score') + \
            variance_score.hvplot.quadmesh(x='lon', y='lat', clim=(0, 1), cmap='RdYlGn', title='Variance score') + \
            bias.hvplot.quadmesh(x='lon', y='lat', clim=(vmin_bias, vmax_bias), cmap='coolwarm', title='Bias') +\
            normalized_bias.hvplot.quadmesh(x='lon', y='lat', clim=(vmin_normalized_bias, vmax_normalized_bias), cmap='coolwarm', title='Normalized Bias')
           ).cols(2)


def plot_timeseries_scores(filename):
    
    ds_score = xr.open_dataset(filename, group='diff')
    rmse_map = ds_score['rms']
    error_variance_map = ds_score['variance']
    std_alongtrack = np.sqrt(error_variance_map)
    
    ds_score = xr.open_dataset(filename, group='alongtrack')
    rms_alongtrack = ds_score['rms']
    variance_alongtrack = ds_score['variance']
    mean_alongtrack = ds_score['mean']
    
    ds_score = xr.open_dataset(filename, group='maps')
    mean_map = ds_score['mean']
    
    
    
    rmse_score = 1. - rmse_map/rms_alongtrack
    rmse_score.name = 'rmse_score'
    
    variance_score = 1. - error_variance_map/variance_alongtrack
    variance_score.name = 'variance_score'
    
    bias = mean_map - mean_alongtrack
    bias.name = 'bias'
    
    normalized_bias = (mean_map - mean_alongtrack)/std_alongtrack
    normalized_bias.name = 'normalized_bias'
    
    return (rmse_score.hvplot.line(x='time', title='RMSE score', ylim=(0, 1), grid=True, color='k') + \
            variance_score.hvplot.line(x='time', title='Variance score', ylim=(0, 1), grid=True, color='royalblue') + \
            bias.hvplot.line(x='time', title='Bias', grid=True, color='mediumaquamarine') +\
            normalized_bias.hvplot.line(x='time', title='Normalized Bias', grid=True, color='salmon')
           ).cols(2)
    