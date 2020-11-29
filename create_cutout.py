# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 19:33:30 2020

@author: Monisha
"""


import atlite
import logging

logging.basicConfig(level=logging.INFO)

cutout = atlite.Cutout(name="europe-2011-01",
                       cutout_dir="C:/Users/Monisha/ProjectLumost",
                       module="era5",
                       xs=slice(-13.6913, 1.7712),
                       ys=slice(49.9096, 60.8479),
                       years=slice(2011, 2011),
                       months=slice(1,1)
                       )
#4: Prepare cutout

