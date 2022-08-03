import xarray as xr
import motuclient
import os

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
    
    
 