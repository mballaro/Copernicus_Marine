import xarray as xr
import motuclient
import os
import pandas as pd
import numpy as np

### https://help.marine.copernicus.eu/en/articles/5182598-how-to-consume-the-opendap-api-and-cas-sso-using-python

#! /usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Copernicus Marine User Support Team"
__copyright__ = "(C) 2021 E.U. Copernicus Marine Service Information"
__credits__ = ["E.U. Copernicus Marine Service Information"]
__license__ = "MIT License - You must cite this source"
__version__ = "202104"
__maintainer__ = "D. Bazin, E. DiMedio, C. Giordan"
__email__ = "servicedesk dot cmems at mercator hyphen ocean dot eu"

def copernicusmarine_datastore(dataset, username, password):
    from pydap.client import open_url
    from pydap.cas.get_cookies import setup_session
    
    cas_url = 'https://cmems-cas.cls.fr/cas/login'
    session = setup_session(cas_url, username, password)
    session.cookies.set("CASTGC", session.cookies.get_dict()['CASTGC'])
    
    database = ['my', 'nrt']
    url = f'https://{database[0]}.cmems-du.eu/thredds/dodsC/{dataset}'
    
    try:
        data_store = xr.backends.PydapDataStore(open_url(url, session=session))
    except:
        url = f'https://{database[1]}.cmems-du.eu/thredds/dodsC/{dataset}'
        data_store = xr.backends.PydapDataStore(open_url(url, session=session))
    
    return data_store




def get_cmems_duacs_alongtrack(mission, time_min, time_max, USERNAME, PASSWORD):
    
    import sys
    # can be imporved with lon/lat selcetion, choice of service -d etc ....
    cmd = f'{sys.executable} -m motuclient --motu https://my.cmems-du.eu/motu-web/Motu --service-id SEALEVEL_GLO_PHY_L3_MY_008_062-DGF --product-id cmems_obs-sl_glo_phy-ssh_my_{mission}-l3-duacs_PT1S --date-min "{time_min}" --date-max "{time_max}" --out-dir ./ --out-name data_{mission}.zip --user {USERNAME} --pwd {PASSWORD}'
        
    os.system(cmd)
    
    os.system(f'mkdir -p data_{mission}')
    cmd = f'unzip data_{mission}.zip -d ./data_{mission}'
    os.system(cmd)
    os.system(f'rm -rf data_{mission}.zip')
    
    ds = xr.open_mfdataset(f'./data_{mission}/*.nc', combine='nested', concat_dim='time')
    
    return ds
    
    
def wget_cmems(product, dataset, time_min, time_max, username, password):
    
    product = 'SEALEVEL_GLO_PHY_L3_MY_008_062'
    dataset = 'cmems_obs-sl_glo_phy-ssh_my_alg-l3-duacs_PT1S'
    cmd = f'wget -r --user={username} --password={password} --no-parent -A*.nc "ftp://{username}@my.cmems-du.eu/Core/SEALEVEL_GLO_PHY_L3_MY_008_062/cmems_obs-sl_glo_phy-ssh_my_alg-l3-duacs_PT1S/2017/01/*.nc"' 
    
    
def list_drifters_of_interest(file_index_history, time_min, time_max, repo_directory):
    
    # Open list of dataset info
    df = pd.read_csv(file_index_history, skiprows=5)
    
    # Read start and end date of each sensor
    start_date = np.asarray([np.datetime64(tt) for tt in df[' time_coverage_start'].values])
    end_date = np.asarray([np.datetime64(tt) for tt in df[' time_coverage_end'].values])
    
    interval_study = pd.Interval(pd.Timestamp(np.datetime64(time_min)), pd.Timestamp(np.datetime64(time_max)))
    
    list_of_drifters = []
    for date_index in range(start_date.size):
    
        interval_drifter = pd.Interval(pd.Timestamp(start_date[date_index]), pd.Timestamp(end_date[date_index]))
    
        if interval_drifter.overlaps(interval_study):
             list_of_drifters.append(date_index)  

    index_drifter_of_interest = np.asarray(list_of_drifters)
    
    list_filename_drifter_of_interest = [repo_directory + os.path.basename(filename) for filename in df[' file_name'][index_drifter_of_interest].values]
    
    return list_filename_drifter_of_interest