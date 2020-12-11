# =============================================================================
# # Import packages
# =============================================================================
import atlite
import xarray as xr
import pandas as pd
import scipy.sparse as sp
import numpy as np
import pgeocode
from collections import OrderedDict

import matplotlib.pyplot as plt
%matplotlib inline

import seaborn as sns
sns.set_style('whitegrid')

# =============================================================================
# # Download Necessary Data Files
# =============================================================================
import requests
import os
import zipfile

def download_file(url, local_filename):
    # variant of http://stackoverflow.com/a/16696317
    if not os.path.exists(local_filename):
        r = requests.get(url, stream=True)
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)
    return local_filename

# =============================================================================
# # Reference time-series
# =============================================================================
opsd_fn = download_file('https://data.open-power-system-data.org/index.php?package=time_series&version=2019-06-05&action=customDownload&resource=3&filter%5B_contentfilter_cet_cest_timestamp%5D%5Bfrom%5D=2012-01-01&filter%5B_contentfilter_cet_cest_timestamp%5D%5Bto%5D=2013-05-01&filter%5BRegion%5D%5B%5D=DE&filter%5BVariable%5D%5B%5D=solar_generation_actual&filter%5BVariable%5D%5B%5D=wind_generation_actual&downloadCSV=Download+CSV',
                        'time_series_60min_singleindex_filtered.csv')

opsd = pd.read_csv(opsd_fn, parse_dates=True, index_col=0)

# we later use the (in current version) timezone unaware datetime64
# to work together with this format, we have to remove the timezone
# timezone information. We are working with UTC everywhere.

opsd.index = opsd.index.tz_convert(None)

# We are only interested in the 2012 data
opsd = opsd[("2011" < opsd.index) & (opsd.index < "2013")]


# =============================================================================
# PV locations (“Anlagenregister”)
# =============================================================================
eeg_fn = download_file('http://www.energymap.info/download/eeg_anlagenregister_2015.08.utf8.csv.zip',
                        'eeg_anlagenregister_2015.08.utf8.csv.zip')

with zipfile.ZipFile(eeg_fn, "r") as zip_ref:
    zip_ref.extract("eeg_anlagenregister_2015.08.utf8.csv")
    
# =============================================================================
# Create a Cutout from ERA5
# =============================================================================
import cartopy.io.shapereader as shpreader
import geopandas as gpd
import userdefined_atlite_functions
shp = shpreader.Reader(shpreader.natural_earth(resolution='10m', category='cultural', name='admin_0_countries'))
de_record = list(filter(lambda c: c.attributes['ISO_A2'] == 'DE', shp.records()))[0]
de = gpd.GeoSeries({**de_record.attributes, 'geometry':de_record.geometry})
x1, y1, x2, y2 = de['geometry'].bounds

cutout_directory = "C:\MasalaAtlite\ATLITE\cutouts"
cutout_name = 'test_cutout'
cutout = atlite.Cutout(name = cutout_name,
                        cutout_dir=cutout_directory,
                        module='era5',
                        xs=slice(x1-.2,x2+.2),
                        ys=slice(y1-.2, y2+.2),
                        #chunks={'time':100},
                        time="2012")
# cutout.prepare()