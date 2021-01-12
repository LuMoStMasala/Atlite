# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 19:33:30 2020

Authors: Ludwig, Monisha, Stefan (a.k.a. "LuMoSt Masala")

Version 2.0.
Changelog:
    12.01.2021:
        Different allocation methods for pv capacity added and compared
    
"""
# =============================================================================
# # Import packages
# =============================================================================
import os
from collections import OrderedDict
import zipfile
import atlite
import pandas as pd
import pgeocode
import matplotlib.pyplot as plt
import seaborn as sns
import requests
import cartopy.io.shapereader as shpreader
import geopandas as gpd
import xarray as xr
import numpy as np
import copy
from operator import itemgetter


# Set seaborn style
sns.set_style('whitegrid')



#Identify whether code is run from script or console
try:
    __file__
except NameError:
    print('This code is executed from console.')
    # This is an IPython command (doesn't work in script)
    # --> It could raise error later when its functionality is used
    #%matplotlib inline # commented out for Pylint code analysis
else:
    print('This code is executed from script.')

#Change working directory (in case it was not set automatically)
if os.path.dirname(os.path.realpath(__file__)) != os.getcwd():
    os.chdir(os.path.dirname(os.path.realpath(__file__)))
    print('Working directory was changed to ' + \
          os.path.dirname(os.path.realpath(__file__)))

# Define download directory
DOWNLOAD_DIR = 'downloads'

# Define cutout name and directory
import platform
if platform.node() == 'T520': # identify Ludwigs Linux Computer
    CUTOUT_DIRECTORY = r'/home/ludwig/AtliteMasala/cutouts'
elif platform.node() == 'stefan-VirtualBox':
    CUTOUT_DIRECTORY = r'/home/stefan/Mastersproject/Cutouts'
elif platform.node() == 'Computername_User_02':
    CUTOUT_DIRECTORY = r'Cutout_direcory for User 02'
else:
    CUTOUT_DIRECTORY = r"C:\MasalaAtlite\ATLITE\cutouts"

if not os.path.isdir(CUTOUT_DIRECTORY):
    os.mkdir(CUTOUT_DIRECTORY)
    print('Local cutout-directory "' + CUTOUT_DIRECTORY + '" created.')
if platform.system()=='Linux':
    DIR_SEP = r"/"
else:
    DIR_SEP = r"\\"
    
# =============================================================================
# CONSTANTS
# =============================================================================
P_TOTAL_SOLAR = 23.75 * 1e3 #MW

# <codecell>
# =============================================================================
# Download Necessary Data Files
# =============================================================================
# Define function to download files
def download_file(url, local_directory, local_filename):
    """Downloads a file if it doesn't exist already"""
    # variant of http://stackoverflow.com/a/16696317
    local_filepath = local_directory + DIR_SEP + local_filename
    if not os.path.isdir(local_directory):
        os.mkdir(local_directory)
        print('Local download directory "' + local_directory + '" created.')
    if not os.path.exists(local_filepath):
        req = requests.get(url, stream=True)
        with open(local_filepath, 'wb') as file:
            for chunk in req.iter_content(chunk_size=1024):
                if chunk:
                    file.write(chunk)
    return local_filepath


# <codecell>
# Installed Capacities in Germany (Anlagenregister) - Zip file
EEG_FN = download_file('http://www.energymap.info/download/eeg_'
                       'anlagenregister_2015.08.utf8.csv.zip',
                       DOWNLOAD_DIR,
                       'eeg_anlagenregister_2015.08.utf8.csv.zip')

# Create ZipFile object 'zip_ref' with generator data
with zipfile.ZipFile(EEG_FN, "r") as zip_ref:
    zip_ref.extract("eeg_anlagenregister_2015.08.utf8.csv",DOWNLOAD_DIR)
    
# <codecell>
# Reference time-series from Open-Power-System-Data (OPSD)
OPSD_FN = download_file('https://data.open-power-system-data.org/index.php?' \
                        'package=time_series&version=2019-06-05&action=' \
                        'customDownload&resource=3&filter%5B_content' \
                        'filter_cet_cest_timestamp%5D%5Bfrom%5D=2012-01' \
                        '-01&filter%5B_contentfilter_cet_cest_timestamp' \
                        '%5D%5Bto%5D=2013-05-01&filter%5BRegion%5D%5B%5' \
                        'D=DE&filter%5BVariable%5D%5B%5D=solar_generati' \
                        'on_actual&filter%5BVariable%5D%5B%5D=wind_gene' \
                        'ration_actual&downloadCSV=Download+CSV',
                        DOWNLOAD_DIR,
                        'time_series_60min_singleindex_filtered.csv')


opsd = pd.read_csv(OPSD_FN, parse_dates=True, index_col=0)

# we later use the (in current version) timezone unaware datetime64
# to work together with this format, we have to remove the timezone
# timezone information. We are working with UTC everywhere.

opsd.index = opsd.index.tz_convert(None)

# We are only interested in the 2012 data
opsd = opsd[(opsd.index > "2011") & (opsd.index < "2013")]

# <codecell>
# =============================================================================
# Create a Cutout from ERA5
# =============================================================================
shp = shpreader.Reader(shpreader.natural_earth(resolution='10m',
                                               category='cultural',
                                               name='admin_0_countries'))
de_record = list(filter(lambda c: c.attributes['ISO_A2'] == 'DE',
                        shp.records()))[0]

# Create vector (Geopandas Series) 'de' with several information about Germany
de = gpd.GeoSeries({**de_record.attributes, 'geometry':de_record.geometry})

# Determining the bounds of the MultiPolygon-geometry of Germany
# (de['geometry']), i.e. convert shape to rectangle coordinates
x1, y1, x2, y2 = de['geometry'].bounds

CUTOUT_NAME = 'Cutout_DE_01'

# Create cutout if it doesn't exist already
if not os.path.isdir(CUTOUT_DIRECTORY + DIR_SEP + CUTOUT_NAME):
    cutout = atlite.Cutout(name = CUTOUT_NAME,
                            cutout_dir=CUTOUT_DIRECTORY,
                            module='era5',
                            xs=slice(x1-.2, x2+.2),
                            ys=slice(y1-.2, y2+.2),
                            #chunks={'time':100},
                            years = slice(2012, 2012))
    print("Cutout '" + CUTOUT_NAME + "' was newly created.")
    cutout.prepare()
else:
    cutout = atlite.Cutout(name = CUTOUT_NAME, cutout_dir = CUTOUT_DIRECTORY)
    print("Cutout '" + CUTOUT_NAME + "' exist already. Loaded successfully.")

# <codecell>
# =============================================================================
# Capacity layout
# =============================================================================
def capacity_layout(cutout_input, typ, cap_range=None, until=None):
    """Aggregate selected capacities to the cutouts grid into a capacity
    layout.

    Parameters
    ----------
    cutout_input : atlite.cutout
        The cutout for which the capacity layout is contructed.
    typ : str
        Type of energy source, e.g. "Solarstrom" (PV), "Windenergie" (wind).
    cap_range : (optional) list-like
        Two entries, limiting the lower and upper range of capacities (in kW)
        to include. Left-inclusive, right-exclusive.
    until : str
        String representation of a datetime object understood by
        pandas.to_datetime()
        for limiting to installations existing until this datetime.

    """

    # Load locations of installed capacities and remove incomplete entries
    cols = OrderedDict((('installation_date', 0),
                        ('plz', 2), ('city', 3),
                        ('type', 6),
                        ('capacity', 8), ('level', 9),
                        ('lat', 19), ('lon', 20),
                        ('validation', 22)))
    database = pd.read_csv(DOWNLOAD_DIR + DIR_SEP + \
                           'eeg_anlagenregister_2015.08.utf8.csv',
                       sep=';', decimal=',', thousands='.',
                       comment='#', header=None,
                       usecols=list(cols.values()),
                       names=list(cols.keys()),
                       # German postal codes can start with '0' so we need to
                       # treat them as str
                       dtype={'plz':str},
                       parse_dates=['installation_date'],
                       na_values=('O04WF', 'keine'))

    # kick out data with without 'OK' or with no PLZ
    database = database[(database['validation'] == 'OK') &
                        (database['plz'].notna())]

    # Query postal codes <-> coordinates mapping
    de_nomi = pgeocode.Nominatim('de')
    plz_coords = de_nomi.query_postal_code(database['plz'].unique())
    plz_coords = plz_coords.set_index('postal_code')

    # Fill missing lat / lon using postal codes entries
    database.loc[database['lat'].isna(), 'lat'] = \
        database['plz'].map(plz_coords['latitude'])
    database.loc[database['lon'].isna(), 'lon'] = \
        database['plz'].map(plz_coords['longitude'])

    # Ignore all locations which have not be determined yet
    database = database[database['lat'].notna() & database['lon'].notna()]

    # Select data based on type (i.e. solar/PV, wind, ...)
    data = database[database['type'] == typ].copy()

    # Optional: Select based on installation day
    if until is not None:
        data = data[data['installation_date'] < pd.to_datetime(until)]

    # Optional: Only installations within this caprange
    # (left inclusive, right exclusive)
    if cap_range is not None:
        data = data[(cap_range[0] <= data['capacity']) &
                    (data['capacity'] < cap_range[1])]

    # Determine nearest cells from cutout_input
    cells = gpd.GeoDataFrame({'geometry': cutout_input.grid_cells,
                              'lon': cutout_input.grid_coordinates()[:,0],
                              'lat': cutout_input.grid_coordinates()[:,1]})

    nearest_cell = cutout_input.meta.sel({'x': data.lon.values,
                                    'y': data.lat.values},
                                   'nearest').coords

    # Map capacities to closest cell coordinate
    data['lon'] = nearest_cell.get('lon').values
    data['lat'] = nearest_cell.get('lat').values

    new_data = data.merge(cells, how='inner')

    # Sum capacities for each grid cell (lat, lon)
    # then: restore lat lon as coumns
    # then: rename and reindex to match cutout_input coordinates
    new_data = new_data.groupby(['lat','lon']).sum()

    # layout = new_data.reset_index().rename(columns={'lat':'y','lon':'x'})\
    #                     .set_index(['y','x']).capacity\
    #                     .to_xarray().reindex_like(cutout_input.data)
    layout = new_data.reset_index().rename(columns={'lat':'y','lon':'x'})\
                    .set_index(['y','x']).capacity\
                    .to_xarray().reindex_like(cutout_input.meta)

    layout = (layout/1e3).fillna(.0).rename('Installed Capacity [MW]')

    return layout
# <codecell>

# Create solar layout with capacities from "Anlagenregister" (2012)
solar_layout = capacity_layout(cutout, 'Solarstrom', until="2012")

# Plot solar layout
solar_layout.plot(cmap="inferno_r", size=8, aspect=1)
plt.title("Installed PV in Germany until 2012")
plt.tight_layout()
# <codecell>

# Create DataFrame with generation data
generation = pd.DataFrame(index=cutout.meta.time.values)

# Add generation time series from real installed capacities from
# Anlagenregister in 2012
generation['Anlagenregister'] = \
    cutout.pv(panel="CSi",
              orientation='latitude_optimal',
              layout=solar_layout
              ).to_pandas().iloc[0]

# Compute capacity factors for PV in 2012
capacityFactors_PV = cutout.pv(panel="CSi",
                               orientation='latitude_optimal',
                               capacity_factor = True
                               )

#%%
# =============================================================================
# Create capacity layout
# =============================================================================
def synthesize_layout(baseLayout_org, capacityFactors, totalCapacity, mode, cellLimit=1e20):
    """
    

    Parameters
    ----------
    baseLayout_org : DataArray
        Existing capacity layout that fits to cutout, for example solar_layout
    capacityFactors : DataArray
        Atlite output by appling capacity_factors=True to pv.
    totalCapacity : int
        Total allocated capacity in MW.
    mode : str
        Mode/Method of allocation. Options:
            'CF highest'
            'CF lowest'
            'even'
            'CF percentage'
            'CF median'
    cellLimit : int
        Optional. Maximum capacity to be allocated to one cell.
            

    Returns
    -------
    Capacity layout.

    """
    baseLayout = copy.deepcopy(baseLayout_org)

    # Create mask of only 0 outside Germany and 1 inside Germany
    baseLayout.data[baseLayout.data > 0] = 1
    
    newLayout = copy.deepcopy(baseLayout)
    
    # Create new layout by using capacity factors and set all values outside
    # of Germany to zero with mask
    layout = capacityFactors * baseLayout.data
    
    # Create array of zeros with shape of layout
    capacityArr = np.zeros(layout.data.shape)
    
    # Allocation according to highest Capacity Factor (CF)
    if mode == 'CF highest':
        cf_sorted = sorted(np.ndenumerate(layout.data),
                           key=itemgetter(1),
                           reverse=True)
        # delete elements with CF = 0
        cf_sorted = [x for x in cf_sorted if x[1] != 0]
        # Iterate over CFs
        for cell in cf_sorted:
            if totalCapacity - capacityArr.sum() >= cellLimit:
                capacityArr[cell[0][0], cell[0][1]] = cellLimit
            else:
                capacityArr[cell[0][0], cell[0][1]] = \
                    totalCapacity - capacityArr.sum()
                break
        
    # Allocation according to lowest Capacity Factor (CF)
    elif mode == 'CF lowest':
        cf_sorted = sorted(np.ndenumerate(layout.data),
                           key=itemgetter(1),
                           reverse=True)
        # delete elements with CF = 0
        cf_sorted = [x for x in cf_sorted if x[1] != 0]
        # Iterate over CFs
        for cell in cf_sorted:
            if totalCapacity - capacityArr.sum() >= cellLimit:
                capacityArr[cell[0][0], cell[0][1]] = cellLimit
            else:
                capacityArr[cell[0][0], cell[0][1]] = \
                    totalCapacity - capacityArr.sum()
                break
        
    # Allocation according to an even distribution over all cells
    elif mode == 'even':
        # number of cells within Germany
        n_cells = (baseLayout.data > 0).sum()
        
        # Assign same capacities to each cell
        capacityArr[baseLayout.data > 0] = totalCapacity / n_cells

    
    # Allocation according to each cell's share of CF sum
    elif mode == 'CF percentage':
        capacityArr = totalCapacity * layout.data / layout.data.sum()
        
    
    # Allocation according to cells higher than median
    elif mode == 'CF median':
        medianBool = layout.data >= np.median(layout.data)
        capacityArr[medianBool] = totalCapacity * \
            (layout.data[medianBool] / layout.data[medianBool].sum())        
        
    # layout.data = capacityArr
    newLayout.data = capacityArr
    print(str((newLayout.data > 0).sum()) + ' out of ' + \
          str((baseLayout.data > 0).sum()) + ' cells were filled with ' + \
              'capacities according to mode: ' + mode)
    
    if newLayout.data.sum() < totalCapacity:
        print('Caution! Allocated capacity is less than total capacity' + \
              ' specified. Please check cell limit.')

    return newLayout


# <codecell>
# =============================================================================
# Evaluate various capacity layouts
# =============================================================================
std_limit = 200
layout_variation = {
    'highest CF, 50 MW':
        synthesize_layout(solar_layout, capacityFactors_PV,
                          P_TOTAL_SOLAR,
                          'CF highest',
                          cellLimit=50),
    'highest CF, 100 MW':
        synthesize_layout(solar_layout, capacityFactors_PV,
                          P_TOTAL_SOLAR,
                          'CF highest',
                          cellLimit=100),
    'highest CF, 250 MW':
        synthesize_layout(solar_layout, capacityFactors_PV,
                          P_TOTAL_SOLAR,
                          'CF highest',
                          cellLimit=250),
    'highest CF, 500 MW':
        synthesize_layout(solar_layout, capacityFactors_PV,
                          P_TOTAL_SOLAR,
                          'CF highest',
                          cellLimit=500),
    'highest CF, 1000 MW':
        synthesize_layout(solar_layout, capacityFactors_PV,
                          P_TOTAL_SOLAR,
                          'CF highest',
                          cellLimit=1000),
    'highest CF, no limit':
        synthesize_layout(solar_layout, capacityFactors_PV,
                          P_TOTAL_SOLAR,
                          'CF highest'),
    'lowest CF, 200 MW':
        synthesize_layout(solar_layout, capacityFactors_PV,
                          P_TOTAL_SOLAR,
                          'CF lowest',
                          cellLimit=std_limit),
    'CF Median':
        synthesize_layout(solar_layout, capacityFactors_PV,
                          P_TOTAL_SOLAR,
                          'CF median'),
    'CF percentage':
        synthesize_layout(solar_layout, capacityFactors_PV,
                          P_TOTAL_SOLAR,
                          'CF percentage'),
    'Even distribution':
        synthesize_layout(solar_layout, capacityFactors_PV,
                          P_TOTAL_SOLAR,
                          'even')
    }

# Compute generation data for all capacity layout variations
for key in layout_variation:
    generation[key] = cutout.pv(panel="CSi",
                                orientation='latitude_optimal',
                                layout=layout_variation[key]).to_pandas().iloc[0]
# Add reference (OPSD) time-series
generation['OPSD'] = opsd['DE_solar_generation_actual']

# Plot weekly averages
generation.divide(1000).resample('1W').mean().plot(figsize=(10,6))
plt.ylabel("Feed-In [GW]")
plt.title('PV time-series (weekly avrg.) in Germany. Total PV Capacity ' + \
          'P_tot_pv: ' + str(round(P_TOTAL_SOLAR/1e3,1)) + ' GW')
plt.tight_layout()
