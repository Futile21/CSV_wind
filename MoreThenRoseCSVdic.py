from netCDF4 import Dataset  # dataset allows us to read the netcdf4 files
import bisect
import os
import pandas as pd
import numpy as np
import time
import xarray as xr


da=xr.open_mfdataset(r"C:\Users\futil\OneDrive\GIZ\Run2\*1.nc")
print(da)