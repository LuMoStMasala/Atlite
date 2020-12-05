# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 19:33:30 2020

@author: Monisha, Ludwig, Stefan (a.k.a. "LuMoSt Masala")
"""


import atlite
import logging
import time
import subprocess

#get current working directory (cwd)
import os
cwd = os.getcwd() #current working directory
cfp = os.path.realpath(__file__) #current file path


logging.basicConfig(level=logging.INFO)

#from define_cutout import cutout
# =============================================================================
# Define cutout
# =============================================================================
# import subprocess
# print("Defining cutout information")
# cu = subprocess.Popen("define_cutout.py", shell=True)

# from subprocess import check_output
# out = check_output(["ntpq", "-p"])

#import subprocess
#p = subprocess.Popen("define_cutout.py", shell=True)
#a = subprocess.check_output('define_cutout.py', shell=True)
# p = subprocess.Popen("define_cutout.py", stdout=subprocess.PIPE)
# result = p.communicate()[0]
# print(result)

cutout = atlite.Cutout(name="germany_2011_v01",
                        cutout_dir="C:\MasalaAtlite\ATLITE\cutouts",
                        module="era5",
                        xs=slice(5.5, 15.5),
                        ys=slice(47.2, 55.1),
                        years=slice(2011, 2011),
                        months=slice(1,1)
                        )

# =============================================================================
# Create/Prepare cutout 
# =============================================================================
#cutout.prepare()

# =============================================================================
# Output
# =============================================================================
cutout
cutout.coords #get coordinates
cutout.data
cutout.wind()
