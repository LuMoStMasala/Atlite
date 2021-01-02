# -*- coding: utf-8 -*-
"""
Created on Fri Nov  6 19:33:30 2020

Main author: Ludwig
Other authors: Monisha, Stefan (a.k.a. "LuMoSt Masala")

Version 1.0.1
Changelog:
    31.12.2020:
        function download_file() was modified to specify download folder
        Therefore, "DOWNLOAD_DIR" was created to define local download dir
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


# Modules maybe needed later
#import xarray as xr
#import scipy.sparse as sp
#import numpy as np

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
elif platform.node() == 'DESKTOP-580IP6U': # identify Ludwigs Win Computer
    CUTOUT_DIRECTORY = r"C:\MasalaAtlite\ATLITE\cutouts"
elif platform.node() == 'Computername_User_01':
    CUTOUT_DIRECTORY = r'Cutout_direcory for User 01'
elif platform.node() == 'Computername_User_02':
    CUTOUT_DIRECTORY = r'Cutout_direcory for User 02'
elif platform.node() == 'Computername_User_03':
    CUTOUT_DIRECTORY = r'Cutout_direcory for User 03'
else:
    CUTOUT_DIRECTORY = r"C:\MasalaAtlite\ATLITE\cutouts"
if not os.path.isdir(CUTOUT_DIRECTORY):
    os.mkdir(CUTOUT_DIRECTORY)
    print('Local cutout-directory "' + CUTOUT_DIRECTORY + '" created.')
    
# Identify separator for paths to make code run on Windows and Linux
if platform.system()=='Linux':
    DIR_SEP = r"/"
else:
    DIR_SEP = r"\\"

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
# Installed Capacities in Germany (Anlagenregister) - Zip file
EEG_FN = download_file('http://www.energymap.info/download/eeg_'
                       'anlagenregister_2015.08.utf8.csv.zip',
                       DOWNLOAD_DIR,
                       'eeg_anlagenregister_2015.08.utf8.csv.zip')

# Create ZipFile object 'zip_ref' with generator data
with zipfile.ZipFile(EEG_FN, "r") as zip_ref:
    zip_ref.extract("eeg_anlagenregister_2015.08.utf8.csv",DOWNLOAD_DIR)

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

solar_layout = capacity_layout(cutout, 'Solarstrom', until="2012")

solar_layout.plot(cmap="inferno_r", size=8, aspect=1)
plt.title("Installed PV in Germany until 2012")
plt.tight_layout()
# <codecell>

# pv = cutout.pv(panel="CSi", orientation={'slope': 30., 'azimuth': 0.},
#                layout=solar_layout)
pv = cutout.pv(panel="CSi",
               orientation='latitude_optimal',
               layout=solar_layout)

# compare = pd.DataFrame(dict(atlite=pv.to_pandas(),
#                             opsd=opsd['DE_solar_generation_actual'])) /1e3
compare = pd.DataFrame(dict(atlite=pv.to_pandas().iloc[0],
                            opsd=opsd['DE_solar_generation_actual'])) /1e3 #GW

compare.resample('1W').mean().plot(figsize=(8,5))
plt.ylabel("Feed-In [GW]")
plt.title('PV time-series Germany 2012')
plt.tight_layout()


# <codecell>
# =============================================================================
# Create layout with even capacity distribution
# =============================================================================

# Required total installed capacity in GW
# installedCap = 49
installedCap = solar_layout.data.sum() / 1e3

# New layout
evenDistribution_GER = solar_layout.data.copy()

# Assign dummy value of 1 to each cell with installed capacity in GER
evenDistribution_GER[evenDistribution_GER > 0] = 1

# Number of cells withing German boarder
n_cells = int(evenDistribution_GER.sum())


# Average capacity for each cell
evenCap = installedCap * 1000 / n_cells


# Assign new even capacity to only the cells that contain solar capacity
evenDistribution_GER[evenDistribution_GER > 0] = evenCap

even_GER_layout = solar_layout.copy()
even_GER_layout.rename('Even Capacity Distribution [MW]')
even_GER_layout.data = evenDistribution_GER

even_GER_layout.plot(cmap="cool", size=8, aspect=1)
plt.title("Even Capacity Distribution within German boarder: " +
          str(round(evenCap,2)) + " MW")
plt.tight_layout()

# Compute electricity generation with atlite
# pv_even = cutout.pv(panel="CSi", orientation={'slope': 30., 'azimuth': 0.},
#                     layout=even_GER_layout)
pv_even = cutout.pv(panel="CSi", orientation='latitude_optimal',
                    layout=even_GER_layout)
pv_even_output = pv_even.to_pandas().iloc[0]

# Determine time step length in hours (h)
dt = pv_even_output.index[1]-pv_even_output.index[0]
dt = dt.days * 24 + dt.seconds // 3600

# Determine total electricity generation
E_tot = pv_even_output.sum() * dt / 1e6 # in TWh

# Determine total installed capacity in GW
P_tot = even_GER_layout.data.sum() / 1e3


# Plot 1-week averaged electricity generation in GW
plt.figure()
pv_even_output.divide(1000).resample('1W').mean().plot(figsize=(8,5))
plt.ylabel("Feed-In [GW]")
plt.title('PV time-series (weekly avrg.) in Germany with equal capacity ' \
          'distribution (' + str(round(evenCap,1)) + \
              " MW).\n Total electricity generation E_tot = " + \
              str(round(E_tot,1)) + " GWh")
plt.tight_layout()

# Output
print('==============================================')
print('Findings:')
print('Number of cells in Germany: ' + str(n_cells))
print('Installed capacity in each cell: ' + str(round(evenCap,2)) + ' MW')
print('Total installed capacity: ' + str(round(P_tot,2)) + ' GW')
print('Total annual electricity generation: ' + str(round(E_tot,3)) + ' TWh/a')

# =============================================================================
# Comparison
# =============================================================================
curve1 = opsd['DE_solar_generation_actual'] # in MW
curve2 = pv.to_pandas().iloc[0] # in MW
curve3 = pv_even.to_pandas().iloc[0] # in MW
compare = pd.DataFrame(dict(opsd=curve1,
                            atlite_real=curve2,
                            atlite_even=curve3)) /1e3 # in GW

plt.figure()
compare.resample('1W').mean().plot(figsize=(8,5))
plt.ylabel("Feed-In [GW]")
plt.title('Comparison of PV time-series in Germany 2012\n' +\
          'Total installed capacity in Atlite P_tot = ' + \
              str(round(P_tot,2)) + ' GW')
plt.tight_layout()

#plt.close('all')

# =============================================================================
# Wind capacity layout
# =============================================================================
wind_layout = capacity_layout(cutout, 'Windkraft', until='2012')
