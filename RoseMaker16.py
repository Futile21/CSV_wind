from netCDF4 import Dataset  # dataset allows us to read the netcdf4 files
import bisect
import os
import pandas as pd
import numpy as np
import time
import xarray as xr

DATA_DIR = r"C:\AllRun"
filelist = os.listdir(DATA_DIR)