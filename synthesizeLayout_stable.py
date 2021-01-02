#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  1 12:08:22 2021

@author: ludwig
"""
# Import modules
import numpy as np
import pandas as pd
import xarray as xr
import cartopy.io.shapereader as shpreader
from shapely.geometry import Point
import matplotlib.pyplot as plt

# =============================================================================
# Set CONSTANTS
# =============================================================================
# Examined country
COUNTRY = 'DE'

# Coordinate boundaries of examined region (rectangular)
BOUNDARIES = {
    'lat_min': 47.2,
    'lat_max': 55.1,
    'lon_min': 5.5,
    'lon_max': 15.5,
    }

# Grid resolution (min RES = 0.001)
RES = 0.1

# =============================================================================
# Prepare empty layout and get shape of country borders
# =============================================================================
# Create empty layout as xarray-DataArray with dimensions & coordinates
# The additional '+ 1e-5' is required to include upper coordinate limit
layout = xr.DataArray(dims = ('lat', 'lon'),
                      coords={'lat': np.arange(BOUNDARIES["lat_min"],
                                               BOUNDARIES["lat_max"] + 1e-5,
                                               RES).round(3),
                              'lon': np.arange(BOUNDARIES["lon_min"],
                                               BOUNDARIES["lon_max"] + 1e-5,
                                               RES).round(3)
                               }
                      )

# Assign zeros to layout
# layout.data = np.zeros((len(layout.coords['lat']), len(layout.coords['lon'])))
print('Layout framework created with ' + str(len(layout.to_series())) + \
      ' cells (lat=' + str(len(layout.coords['lat'].values)) + ', lon=' + \
      str(len(layout.coords['lon'].values)) + ').')

# Get shapes
shp = shpreader.Reader(shpreader.natural_earth(resolution='10m',
                                               category='cultural',
                                               name='admin_0_countries'))
# Extract shape of certain country border
border = list(filter(lambda c: c.attributes['ISO_A2'] == COUNTRY,
                        shp.records()))[0].geometry

# Check if country border does not exceed initial defined map coordinates
if border.bounds[0] < BOUNDARIES["lon_min"] or \
   border.bounds[1] < BOUNDARIES["lat_min"] or \
   border.bounds[2] > BOUNDARIES["lon_max"] or \
   border.bounds[3] > BOUNDARIES["lat_max"]:
    print('The examined country ' + COUNTRY + ' is exceeding the initial' +
          ' defined layout boundaries. This may cause errors.')

# =============================================================================
# Create a layout mask
# =============================================================================
# Iteratre through layout cells and compare each with country border limits
# to see if capacity should be assigned to it later.
def maskLayout(baseLayout, shape):
    """
    Create a mask on layout in regards to given shape. 

    Parameters
    ----------
    layout : xarray.core.dataarray.DataArray
        In rectangular shape.
    shape : shapely.geometry.multipolygon.MultiPolygon
        Gives information on shape in which capacity should be assigned.

    Returns
    -------
    layout : xarray.core.dataarray.DataArray
        In rectangular shape with Ones (1) for coordinates within the shape/
        border and Zeros (0) outside of the shape.

    """
    count = 0
    count_ones = 0
    n_cells = len(baseLayout.coords['lat'].values) * \
              len(baseLayout.coords['lon'].values)
    n_prints = 10
    for x in baseLayout.coords['lat'].values:
        for y in baseLayout.coords['lon'].values:
            if shape.contains(Point(y, x)): # ! Reversed lat and lon ;-)
                # This is where the magic happens
                baseLayout.loc[x, y] = 1
                count_ones += 1
            else:
                baseLayout.loc[x, y] = 0
            count += 1
            if count % int(n_cells / n_prints) == 0 or count == n_cells:
                print(str(count).zfill(len(str(n_cells))) + ' of ' + \
                      str(n_cells) + ' cells evaluated (' + \
                      str(int(100 * count/n_cells)) + '%).')
    print('Done. ' + str(count_ones) + ' cells were filled with ones.')
    return baseLayout

layout = maskLayout(layout, border)
layout.rename(lon='x', lat='y')
layout.plot()
