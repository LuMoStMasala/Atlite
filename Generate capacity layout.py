# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 17:48:15 2020

@author: Monisha
"""

import atlite
import xarray as xr
import pandas as pd
import scipy.sparse as sp
import numpy as np

import pgeocode #pip install pgeocode
from collections import OrderedDict

import matplotlib.pyplot as plt  #pip innstall matplotlib
%matplotlib inline

import seaborn as sns  #pip install seaborn
sns.set_style('whitegrid')

import requests
import zipfile

import cartopy.io.shapereader as shpreader
import geopandas as gpd
import logging
import os

shp = shpreader.Reader(shpreader.natural_earth(resolution='10m', category='cultural', name='admin_0_countries'))
de_record = list(filter(lambda c: c.attributes['ISO_A2'] == 'DE', shp.records()))[0]
de = gpd.GeoSeries({**de_record.attributes, 'geometry':de_record.geometry})
x1, y1, x2, y2 = de['geometry'].bounds

cwd = os.getcwd() #current working directory
cfp = os.path.realpath(__file__) 

logging.basicConfig(level=logging.INFO)


cutout = atlite.Cutout(name='germany-2012',
                       module="era5",
                       cutout_dir="C:/Users/Monisha/MasterProjectAtlite",
                       xs=slice(x1-.2,x2+.2), 
                       ys=slice(y1-.2, y2+.2),
                       years=slice(2012, 2012),
                       months=slice(1,1)
                       )
def capacity_layout(cutout, typ, cap_range=None, until=None):
    """Aggregate selected capacities to the cutouts grid into a capacity layout.

    Parameters
    ----------
        cutout : atlite.cutout
            The cutout for which the capacity layout is contructed.
        typ : str
            Type of energy source, e.g. "Solarstrom" (PV), "Windenergie" (wind).
        cap_range : (optional) list-like
            Two entries, limiting the lower and upper range of capacities (in kW)
            to include. Left-inclusive, right-exclusive.
        until : str
            String representation of a datetime object understood by pandas.to_datetime()
            for limiting to installations existing until this datetime.

    """
    
    # Load locations of installed capacities and remove incomplete entries
    cols = OrderedDict((('installation_date', 0),
                        ('plz', 2), ('city', 3),
                        ('type', 6),
                        ('capacity', 8), ('level', 9),
                        ('lat', 19), ('lon', 20),
                        ('validation', 22)))
    database = pd.read_csv('eeg_anlagenregister_2015.08.utf8.csv',
                       sep=';', decimal=',', thousands='.',
                       comment='#', header=None,
                       usecols=list(cols.values()),
                       names=list(cols.keys()),
                       # German postal codes can start with '0' so we need to treat them as str
                       dtype={'plz':str},
                       parse_dates=['installation_date'],
                       na_values=('O04WF', 'keine'))

    database = database[(database['validation'] == 'OK') & (database['plz'].notna())]

    # Query postal codes <-> coordinates mapping
    de_nomi = pgeocode.Nominatim('de')
    plz_coords = de_nomi.query_postal_code(database['plz'].unique())
    plz_coords = plz_coords.set_index('postal_code')

    # Fill missing lat / lon using postal codes entries
    database.loc[database['lat'].isna(), 'lat'] = database['plz'].map(plz_coords['latitude'])
    database.loc[database['lon'].isna(), 'lon'] = database['plz'].map(plz_coords['longitude'])

    # Ignore all locations which have not be determined yet
    database = database[database['lat'].notna() & database['lon'].notna()]

    # Select data based on type (i.e. solar/PV, wind, ...)
    data = database[database['type'] == typ].copy()

    # Optional: Select based on installation day
    if until is not None:
        data = data[data['installation_date'] < pd.to_datetime(until)]

    # Optional: Only installations within this caprange (left inclusive, right exclusive)
    if cap_range is not None:
        data = data[(cap_range[0] <= data['capacity']) & (data['capacity'] < cap_range[1])]

    # Determine nearest cells from cutout
    cells = gpd.GeoDataFrame({'geometry': cutout.grid_cells,
                              'lon': cutout.grid_coordinates()[:,0],
                              'lat': cutout.grid_coordinates()[:,1]})

    nearest_cell = cutout.data.sel({'x': data.lon.values,
                                    'y': data.lat.values},
                                   'nearest').coords

    # Map capacities to closest cell coordinate
    data['lon'] = nearest_cell.get('lon').values
    data['lat'] = nearest_cell.get('lat').values

    new_data = data.merge(cells, how='inner')

    # Sum capacities for each grid cell (lat, lon)
    # then: restore lat lon as coumns
    # then: rename and reindex to match cutout coordinates
    new_data = new_data.groupby(['lat','lon']).sum()

    layout = new_data.reset_index().rename(columns={'lat':'y','lon':'x'})\
                        .set_index(['y','x']).capacity\
                        .to_xarray().reindex_like(cutout.data)

    layout = (layout/1e3).fillna(.0).rename('Installed Capacity [MW]')

    return layout



solar_layout = capacity_layout(cutout, 'Solarstrom',cap_range= None until="2012")

solar_layout.plot(cmap="inferno_r", size=8, aspect=1)
plt.title("Installed PV in Germany until 2012")
plt.tight_layout()