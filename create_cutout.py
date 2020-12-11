# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 19:33:30 2020

Main author: Ludwig
Other authors: Monisha, Stefan (a.k.a. "LuMoSt Masala")
"""


import atlite
import logging
import userdefined_atlite_functions

# Get current working directory (cwd)
import os
cwd = os.getcwd() #current working directory
cfp = os.path.realpath(__file__) #current file path

# Set logging type
logging.basicConfig(level=logging.INFO)

#introduce variable that flips every time the script runs
if 'flipVar' in globals():
    if flipVar == True:
        flipVar = False
    else:
        flipVar = True
else:
    flipVar = True
    
# Define whether new cutout should be created
newCutout = False

if newCutout == False:
    loadCutout = True
else:
    loadCutout = False

# Define cutout directory and name
cutout_directory = 'C:\MasalaAtlite\ATLITE\cutouts'
cutout_name = 'newCutout_name_example23'

# Check cutout directory (not working yet)
#check_cutout_dir(cutout_std_directory, cutout_name)




# =============================================================================
# Create/Prepare cutout 
# =============================================================================
if flipVar == True and newCutout == True:
    cutout = atlite.Cutout(name = cutout_name,
                            cutout_dir = cutout_directory,
                            module = "era5",
                            xs = slice(-13.6913, 1.7712),
                            ys = slice(49.9096, 60.8479),
                            years = slice(2011, 2011),
                            months = slice(1,1)
                            )

elif flipVar == False and 'cutout' in globals():
    # only exectued if cutout variable was already created
    cutout.prepare()
    
# =============================================================================
# Load cutout
# =============================================================================
if loadCutout == True:
    cutout_Loading_test = atlite.Cutout(name = cutout_name, cutout_dir = cutout_directory)
