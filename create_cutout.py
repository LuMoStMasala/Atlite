# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 19:33:30 2020

@author: Monisha, Ludwig, Stefan
"""


import atlite
import logging

#get current working directory (cwd)
import os
cwd = os.getcwd() #current working directory
cfp = os.path.realpath(__file__) #current file path


logging.basicConfig(level=logging.INFO)

# Create cutout
cutout = atlite.Cutout(name="europe-2011-01-v2",
                       cutout_dir="C:\MasalaAtlite\ATLITE\cutouts",
                       module="era5",
                       xs=slice(-13.6913, 1.7712),
                       ys=slice(49.9096, 60.8479),
                       years=slice(2011, 2011),
                       months=slice(1,1)
                       )
#4: Prepare cutout

#cutout.prepare()
