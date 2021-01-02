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
    'lon_min': 5.5,
    'lon_max': 15.5,
    'lat_min': 47.2,
    'lat_max': 55.1,
    }

# Grid resolution (min RES = 0.001)
RES = 0.1

#Layout name
lname = 'Installed (artificial) capacity MW'

# =============================================================================
# Prepare empty layout and get shape of country borders
# =============================================================================
# Create empty layout as xarray-DataArray with dimensions & coordinates
# The additional '+ 1e-5' is required to include upper coordinate limit
layout = xr.DataArray(name=lname,
                      dims = ('lon', 'lat'),
                      coords={'lon': np.arange(BOUNDARIES["lon_min"],
                                               BOUNDARIES["lon_max"] + 1e-5,
                                               RES).round(3),
                              'lat': np.arange(BOUNDARIES["lat_min"],
                                               BOUNDARIES["lat_max"] + 1e-5,
                                               RES).round(3)
                               }
                      )

# Assign zeros to layout
# layout.data = np.zeros((len(layout.coords['lat']), len(layout.coords['lon'])))
print('Layout framework created with ' + str(len(layout.to_series())) + \
      ' cells (lon(x)=' + str(len(layout.coords['lon'].values)) + \
      ', lat(y)=' + str(len(layout.coords['lat'].values)) + ').')

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
    #Rename dimensions from x/y to lon/lat if necessary
    if 'y' in baseLayout.dims:
        baseLayout = baseLayout.rename(y='lat')
    if 'x' in baseLayout.dims:
        baseLayout = baseLayout.rename(x='lon')
        
    count = 0
    count_ones = 0
    n_cells = len(baseLayout.coords['lon'].values) * \
              len(baseLayout.coords['lat'].values)
    n_prints = 10
    
    # Iterate over cells
    for x in baseLayout.coords['lon'].values:
        for y in baseLayout.coords['lat'].values:
            if shape.contains(Point(x, y)):
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
    
    # Rename dimensions to lon=x and lat=y to be aligned with atlite cutouts
    baseLayout = baseLayout.rename(lon='x', lat='y')
    
    return baseLayout

# Apply function i.e. assign Ones in border shape and Zeros outside of borders
layout = maskLayout(layout, border)

#Rename dimensions from x/y to lon/lat if necessary
if 'lon' in layout.dims:
    layout = layout.rename(lon='x')
if 'lat' in layout.dims:
    layout = layout.rename(lat='y')

layout.plot()

# <codecell>
# =============================================================================
# New approach by using cutout coords and solar layout as basis
# =============================================================================

# Load cutout
CUTOUT_NAME = 'Cutout_DE_01'
CUTOUT_DIR = '/home/ludwig/AtliteMasala/cutouts'
cutout = atlite.Cutout(name = CUTOUT_NAME,
                       cutout_dir = CUTOUT_DIR)
print("Cutout '" + CUTOUT_NAME + "' exist already. Loaded successfully.")

# Run atlite pv generation
# atlite_op = cutout.pv(panel="CSi",orientation='latitude_optimal',layout=layout)
# generation = atlite_op.to_series()


import createSolarLayout
layout_s = createSolarLayout.createSolarLayout()

# List of capacities with coords
capList = layout_s.to_dataframe()

# Again, check coordinates if located in Border Shape (Germany)
count = 0
count_ones = 0
printInfo = False
for idx in capList.index:
    count += 1
    x = idx[1]
    y = idx[0]
    # cap = capList.loc[(capList.index.get_level_values('x') == x) & \
    #                   (capList.index.get_level_values('y') == y)][0]
    cap = capList.loc[idx, capList.columns[0]] # Installed Capacity [MW]
    if printInfo == True:
        print('\nComputing index nr: ' + str(count))
        print('x = ' + str(round(x, 4)))
        print('y = ' + str(round(y, 4)))
        print('Capacity = ' + str(round(cap, 4)))
    
    if shape.contains(Point(x, y)):
        capList.loc[idx, capList.columns[0]] = 100
        count_ones += 1
    else:
        capList.loc[idx, capList.columns[0]] = 0
    # if count >= 2:
    #     break

# new capacity distribution
newDist = capList.to_xarray().to_array().values[0]

# Assign new capacity distribution to layout
layout_syn = layout_s.copy()
layout_syn.data = newDist

# Run atlite again
atlite_op = cutout.pv(panel="CSi",
                      orientation='latitude_optimal',
                      layout=layout_syn)

# Extract generation time-series
generation = atlite_op.to_pandas().T
generation.rename(columns={0: "Generation [MW]"}, inplace=True)

generation.plot()