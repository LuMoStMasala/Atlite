# -*- coding: utf-8 -*-
"""
Created on Mon Dec 21 09:52:23 2020

LEARNING XARRAY

@author: Lui
"""
# Import necessary modules
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib as plt

# =============================================================================
# Create xArray using constructors
# =============================================================================
testArr = xr.DataArray(range(5), name='myName', dims=['dimensions'])
testArr
solar_layout


# Create equal distribution array
dimx = 39
dimy = 32

# capacity value for even distribution
evenCap = 2.6

data = np.full((dimy, dimx), evenCap, dtype=float)

evenLayout = xr.DataArray(data, name="Even wind capacity layout in GW", dims=("x", "y"))


example = xr.DataArray(np.random.randn(2, 3), dims=("x", "y"), coords={"x": [10, 20]})



myLayout = xr.DataArray(range(0), name='wind capacity layout', dims=['x', 'y'])
