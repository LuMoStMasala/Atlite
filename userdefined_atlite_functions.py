# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 19:34:32 2020

User defined functions for Atlite Master's Project'
@author: Ludwig
"""
def check_cutout_dir(cutout_std_directory, cutout_name):
    cutout_path = cutout_std_directory + "\\" + cutout_name
    import os
    #from os import path
    if path.isdir(cutout_path):
        answer = input("The cutout " + cutout_name + " already exists. Delete it? (y/n):")
        if answer in ['y','Y']:
            #delete this cutout folder
            try:
                os.rmdir(cutout_path)
            except OSError as e:
                print("Error: %s : %s" % (cutout_path, e.strerror))
        else:
            print('The existing cutout was not deleted. Errors may occur :(')
    else:
        print('Everything is fine.... :)')
        
def test_if_udfFun_Imported():
    print('The file "userdefined_atlite_functions.py" is already loaded.')


def open_atlite_directory():
# Opens the directory of atlite module in file explorer
    import os
    import sys
    if 'atlite' in sys.modules:
        package_dir = os.path.dirname(os.path.realpath(atlite.__file__))
        os.startfile(package_dir)
    else:
        print('Atlite module was not imported.')