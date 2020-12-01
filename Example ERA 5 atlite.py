# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 10:55:07 2020

This sheet reproduces an example of renewable power prediction using atlite to
improve my understanding of this package

@author: Stefan Burbach
"""
# =============================================================================
# Import packages
# =============================================================================
import atlite
import logging

#define logging
logging.basicConfig(level = logging.DEBUG)

# =============================================================================
# Defining the cutout extent
# =============================================================================
#cutout_dir = 'C:/Benutzer/stefa/data/cutouts'
# cutout = atlite.Cutout(path="western-europe-2011-01.nc",
#                        module="era5",
#                        x=slice(-13.6913, 1.7712),
#                        y=slice(49.9096, 60.8479),
#                        time="2011-01"
#                       )
# cutout = atlite.Cutout(name ="western-europe-2011-01.nc",
#                        module = "era5",
#                        xs = slice(-13.6913, 1.7712),
#                        ys = slice(49.9096, 60.8479),
#                        years = slice(2011,2011),
#                        months = slice(1,1)
#                        )

#%%
cutout = atlite.Cutout(name="europe-2011-01",
                       cutout_dir = "C:/Users/stefa/Studium/Master Sustainable Systems Engineering/3.Semester/Masters Project/Code/Example ERA 5",
                       module="era5",
                       xs=slice(-13.6913, 1.7712),
                       ys=slice(49.9096, 60.8479),
                       years=slice(2011, 2011),
                       months=slice(1,1)
                       )

#%%
#run the following line seperately
cutout.prepare()