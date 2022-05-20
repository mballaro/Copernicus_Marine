import xarray as xr
import motuclient

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





class MotuOptions:
    def __init__(self, attrs: dict):
        super(MotuOptions, self).__setattr__("attrs", attrs)

    def __setattr__(self, k, v):
        self.attrs[k] = v

    def __getattr__(self, k):
        try:
            return self.attrs[k]
        except KeyError:
            return None
        
        
script_template = 'python -m motuclient --motu https://nrt.cmems-du.eu/motu-web/Motu --service-id GLOBAL_ANALYSIS_FORECAST_PHY_001_024-TDS --product-id global-analysis-forecast-phy-001-024 --longitude-min 100 --longitude-max 155 --latitude-min 10 --latitude-max 40 --date-min "2021-02-21 12:00:00" --date-max "2021-02-27 12:00:00" --depth-min 0.493 --depth-max 0.4942 --variable thetao --out-dir <OUTPUT_DIRECTORY> --out-name <OUTPUT_FILENAME> --user <USERNAME> --pwd <PASSWORD>'

#! /usr/bin/env python3
# -*- coding: utf-8 -*-
__author__ = "Copernicus Marine User Support Team"
__copyright__ = "(C) 2021 E.U. Copernicus Marine Service Information"
__credits__ = ["E.U. Copernicus Marine Service Information"]
__license__ = "MIT License - You must cite this source"
__version__ = "202105"
__maintainer__ = "D. Bazin, E. DiMedio, J. Cedillovalcarce, C. Giordan"
__email__ = "servicedesk dot cmems at mercator hyphen ocean dot eu"

def motu_option_parser(script_template, usr, pwd, output_filename):
    dictionary = dict(
        [e.strip().partition(" ")[::2] for e in script_template.split('--')])
    dictionary['variable'] = [value for (var, value) in [e.strip().partition(" ")[::2] for e in script_template.split('--')] if var == 'variable']  # pylint: disable=line-too-long
    for k, v in list(dictionary.items()):
        if v == '<OUTPUT_DIRECTORY>':
            dictionary[k] = '.'
        if v == '<OUTPUT_FILENAME>':
            dictionary[k] = output_filename
        if v == '<USERNAME>':
            dictionary[k] = usr
        if v == '<PASSWORD>':
            dictionary[k] = pwd
        if k in ['longitude-min', 'longitude-max', 'latitude-min', 
                 'latitude-max', 'depth-min', 'depth-max']:
            dictionary[k] = float(v)
        if k in ['date-min', 'date-max']:
            dictionary[k] = v[1:-1]
        dictionary[k.replace('-','_')] = dictionary.pop(k)
    dictionary.pop('python')
    dictionary['auth_mode'] = 'cas'
    return dictionary


def get_cmems_duacs_alongtrack(mission, time_min, time_max, OUTPUT_FILENAME, USERNAME, PASSWORD):
    
    data_request_options_dict_manual = {
    "service_id": "SEALEVEL_GLO_PHY_L3_MY_008_062-DGF",
    "product_id": f"cmems_obs-sl_glo_phy-ssh_my_{mission}-l3-duacs_PT1S",
    "date_min": f"{time_min}",
    "date_max": f"{time_max}",
#     "longitude_min": 0.,
#     "longitude_max": 360.,
#     "latitude_min": -90.,
#     "latitude_max": 90.,
#     "depth_min": 0.493,
#     "depth_max": 0.4942,
#     "variable": ["sla_unfiltered"],
    "motu": "https://my.cmems-du.eu/motu-web/Motu",
    "out_dir": "./",
    "out_name": OUTPUT_FILENAME,
    "auth_mode": "cas",
    "user": USERNAME,
    "pwd": PASSWORD
    }
    
    #data_request_options_dict_automated = motu_option_parser(script_template, USERNAME, PASSWORD, OUTPUT_FILENAME)
    
    motuclient.motu_api.execute_request(MotuOptions(data_request_options_dict_manual))