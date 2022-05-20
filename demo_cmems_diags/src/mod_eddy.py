from datetime import datetime
from py_eddy_tracker.dataset.grid import RegularGridDataset

import tempfile
import numpy as np
import os



def eddy_detection(dataset, pixel_limit=(5, 2000), shape_error=70, sampling=20):
        
#     # Read Dataset
#     coordinates = ('longitude', 'latitude')
#     regular_grid = RegularGridDataset.with_array(coordinates=coordinates,
#                                                  datas={
#                                                         'adt': np.ma.masked_invalid(dataset['adt'][:, :].values.T),
#                                                         coordinates[0]: dataset['longitude'].values,
#                                                         coordinates[1]: dataset['latitude'].values,
#                                                         },
#                                                 )
    
#     regular_grid.variables_description['adt']['attrs']['units'] = 'm'
    
    tmp_filename = tempfile.NamedTemporaryFile(dir='./').name + '.nc'
    
    encoding = {"longitude": {"dtype": "float32"},
                "latitude": {"dtype": "float32"},
                "adt": {"dtype": "int32", "scale_factor": 0.0001, "zlib": True, "_FillValue": -2147483647},
               }
    
    dataset.to_netcdf(tmp_filename)#, encoding=encoding)
    regular_grid = RegularGridDataset(tmp_filename, "longitude", "latitude",)
    try:
        date = datetime.strptime(np.datetime_as_string(dataset['time'].values[0]), '%Y-%m-%dT%H:%M:%S.000000000')   # For glorys
    except IndexError:
        date = datetime.strptime(np.datetime_as_string(dataset['time'].values), '%Y-%m-%dT%H:%M:%S.000000000')     # for duacs
    
    # Compute geostrophic currents from adt
    regular_grid.add_uv("adt", "ugos_pyeddy", "vgos_pyeddy")
    
    # filter large scale
    wavelength = 400
    regular_grid.copy("adt", "adt_raw")
    regular_grid.bessel_high_filter("adt", wavelength, order=1)
    
    # Detect anticyclonic and cyclonic eddies
    a_adt, c_adt = regular_grid.eddy_identification("adt", 
                                                    "ugos_pyeddy", 
                                                    "vgos_pyeddy", 
                                                    date, 
                                                    0.002, 
                                                    pixel_limit=pixel_limit,
                                                    shape_error=shape_error,
                                                    sampling=sampling )
    
    os.system(f'rm -rf {tmp_filename}')
    
    return regular_grid, a_adt, c_adt