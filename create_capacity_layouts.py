# <codecell>
# =============================================================================
# # Import packages
# =============================================================================
import atlite
import xarray as xr
import pandas as pd
import scipy.sparse as sp
import numpy as np
import pgeocode
from collections import OrderedDict

import matplotlib.pyplot as plt
#Identify whether code is run from script or console
try:
    __file__
except:
    print('This code is executed from console.')
    # This is an IPython command (doesn't work in script)
    # --> It could raise error later when its functionality is used
    %matplotlib inline 
else:
    print('This code is executed from script.')


import seaborn as sns
sns.set_style('whitegrid')

#Change working directory (in case it was not set automatically)
import os
if os.path.dirname(os.path.realpath(__file__)) != os.getcwd():
    os.chdir(os.path.dirname(os.path.realpath(__file__)))
    print('Working directory was changed to ' + os.path.dirname(os.path.realpath(__file__)))


# <codecell>
# =============================================================================
# Download Necessary Data Files
# =============================================================================
import requests
import os
import zipfile
# Standard download folder
download_dir = "C:\MasalaAtlite\ATLITE\downloads"
# Define function to download files
def download_file(url, local_filename):
    # variant of http://stackoverflow.com/a/16696317
    if not os.path.exists(local_filename):
        r = requests.get(url, stream=True)
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)
    return local_filename

# <codecell>
# =============================================================================
# # Reference time-series
# =============================================================================
opsd_fn = download_file('https://data.open-power-system-data.org/index.php?package=time_series&version=2019-06-05&action=customDownload&resource=3&filter%5B_contentfilter_cet_cest_timestamp%5D%5Bfrom%5D=2012-01-01&filter%5B_contentfilter_cet_cest_timestamp%5D%5Bto%5D=2013-05-01&filter%5BRegion%5D%5B%5D=DE&filter%5BVariable%5D%5B%5D=solar_generation_actual&filter%5BVariable%5D%5B%5D=wind_generation_actual&downloadCSV=Download+CSV',
                        'time_series_60min_singleindex_filtered.csv')

opsd = pd.read_csv(opsd_fn, parse_dates=True, index_col=0)

# we later use the (in current version) timezone unaware datetime64
# to work together with this format, we have to remove the timezone
# timezone information. We are working with UTC everywhere.

opsd.index = opsd.index.tz_convert(None)

# We are only interested in the 2012 data
opsd = opsd[("2011" < opsd.index) & (opsd.index < "2013")]

# <codecell>
# =============================================================================
# PV locations (“Anlagenregister”)
# =============================================================================
eeg_fn = download_file('http://www.energymap.info/download/eeg_anlagenregister_2015.08.utf8.csv.zip',
                        'eeg_anlagenregister_2015.08.utf8.csv.zip')

# Create ZipFile object 'zip_ref' with generator data
with zipfile.ZipFile(eeg_fn, "r") as zip_ref:
    zip_ref.extract("eeg_anlagenregister_2015.08.utf8.csv")
 
# <codecell>
# =============================================================================
# Create a Cutout from ERA5
# =============================================================================
import cartopy.io.shapereader as shpreader
import geopandas as gpd
import userdefined_atlite_functions
shp = shpreader.Reader(shpreader.natural_earth(resolution='10m', category='cultural', name='admin_0_countries'))
de_record = list(filter(lambda c: c.attributes['ISO_A2'] == 'DE', shp.records()))[0]

# Create vector (Geopandas Series) 'de' with several information about Germany
de = gpd.GeoSeries({**de_record.attributes, 'geometry':de_record.geometry})
# Determining the bounds of the MultiPolygon-geometry of Germany (de['geometry']), i.e. convert shape to rectangle coordinates
x1, y1, x2, y2 = de['geometry'].bounds

# Define cutout name and directory
cutout_directory = "C:\MasalaAtlite\ATLITE\cutouts"
cutout_name = 'Cutout_DE_01'
# Create cutout if it doesn't exist already
if not os.path.isdir(cutout_directory + "\\" + cutout_name):
    cutout = atlite.Cutout(name = cutout_name,
                            cutout_dir=cutout_directory,
                            module='era5',
                            xs=slice(x1-.2, x2+.2),
                            ys=slice(y1-.2, y2+.2),
                            #chunks={'time':100},
                            years = slice(2012, 2012))
    print("Cutout '" + cutout_name + "' was newly created.")
    cutout.prepare()
else:
    cutout = atlite.Cutout(name = cutout_name, cutout_dir = cutout_directory)
    print("Cutout '" + cutout_name + "' exist already. Loaded successfully.")
    
# <codecell>


# <codecell>
# =============================================================================
# Capacity layout
# =============================================================================
def capacity_layout(cutout, typ, cap_range=None, until=None):
    """Aggregate selected capacities to the cutouts grid into a capacity layout.

    Parameters
    ----------
        cutout : atlite.cutout
            The cutout for which the capacity layout is contructed.
        typ : str
            Type of energy source, e.g. "Solarstrom" (PV), "Windenergie" (wind).
        cap_range : (optional) list-like
            Two entries, limiting the lower and upper range of capacities (in kW)
            to include. Left-inclusive, right-exclusive.
        until : str
            String representation of a datetime object understood by pandas.to_datetime()
            for limiting to installations existing until this datetime.

    """

    # Load locations of installed capacities and remove incomplete entries
    cols = OrderedDict((('installation_date', 0),
                        ('plz', 2), ('city', 3),
                        ('type', 6),
                        ('capacity', 8), ('level', 9),
                        ('lat', 19), ('lon', 20),
                        ('validation', 22)))
    database = pd.read_csv('eeg_anlagenregister_2015.08.utf8.csv',
                       sep=';', decimal=',', thousands='.',
                       comment='#', header=None,
                       usecols=list(cols.values()),
                       names=list(cols.keys()),
                       # German postal codes can start with '0' so we need to treat them as str
                       dtype={'plz':str},
                       parse_dates=['installation_date'],
                       na_values=('O04WF', 'keine'))
    
    # kick out data with without 'OK' or with no PLZ
    database = database[(database['validation'] == 'OK') & (database['plz'].notna())]

    # Query postal codes <-> coordinates mapping
    de_nomi = pgeocode.Nominatim('de')
    plz_coords = de_nomi.query_postal_code(database['plz'].unique())
    plz_coords = plz_coords.set_index('postal_code')

    # Fill missing lat / lon using postal codes entries
    database.loc[database['lat'].isna(), 'lat'] = database['plz'].map(plz_coords['latitude'])
    database.loc[database['lon'].isna(), 'lon'] = database['plz'].map(plz_coords['longitude'])

    # Ignore all locations which have not be determined yet
    database = database[database['lat'].notna() & database['lon'].notna()]

    # Select data based on type (i.e. solar/PV, wind, ...)
    data = database[database['type'] == typ].copy()

    # Optional: Select based on installation day
    if until is not None:
        data = data[data['installation_date'] < pd.to_datetime(until)]

    # Optional: Only installations within this caprange (left inclusive, right exclusive)
    if cap_range is not None:
        data = data[(cap_range[0] <= data['capacity']) & (data['capacity'] < cap_range[1])]

    # Determine nearest cells from cutout
    cells = gpd.GeoDataFrame({'geometry': cutout.grid_cells,
                              'lon': cutout.grid_coordinates()[:,0],
                              'lat': cutout.grid_coordinates()[:,1]})

    nearest_cell = cutout.meta.sel({'x': data.lon.values,
                                    'y': data.lat.values},
                                   'nearest').coords

    # Map capacities to closest cell coordinate
    data['lon'] = nearest_cell.get('lon').values
    data['lat'] = nearest_cell.get('lat').values

    new_data = data.merge(cells, how='inner')

    # Sum capacities for each grid cell (lat, lon)
    # then: restore lat lon as coumns
    # then: rename and reindex to match cutout coordinates
    new_data = new_data.groupby(['lat','lon']).sum()

    # layout = new_data.reset_index().rename(columns={'lat':'y','lon':'x'})\
    #                     .set_index(['y','x']).capacity\
    #                     .to_xarray().reindex_like(cutout.data)
    layout = new_data.reset_index().rename(columns={'lat':'y','lon':'x'})\
                    .set_index(['y','x']).capacity\
                    .to_xarray().reindex_like(cutout.meta)

    layout = (layout/1e3).fillna(.0).rename('Installed Capacity [MW]')

    return layout

solar_layout = capacity_layout(cutout, 'Solarstrom', until="2012")

solar_layout.plot(cmap="inferno_r", size=8, aspect=1)
plt.title("Installed PV in Germany until 2012")
plt.tight_layout()
# <codecell>

pv = cutout.pv(panel="CSi", orientation={'slope': 30., 'azimuth': 0.}, layout=solar_layout)

#compare = pd.DataFrame(dict(atlite=pv.to_pandas(), opsd=opsd['DE_solar_generation_actual'])) /1e3 # in GW
compare = pd.DataFrame(dict(atlite=pv.to_pandas().iloc[0], opsd=opsd['DE_solar_generation_actual'])) /1e3 # in GW

compare.resample('1W').mean().plot(figsize=(8,5))
plt.ylabel("Feed-In [GW]")
plt.title('PV time-series Germany 2012')
plt.tight_layout()


# <codecell>
# =============================================================================
# Create layout with even capacity distribution
# =============================================================================
# Create equal distribution array
dimx = 39
dimy = 32

# capacity value for even distribution
evenCap = 2.6

# # new Data: even capacity distribution
# evenDistribution = np.full((dimy, dimx), evenCap, dtype=float)

# # new layout (from copy of solar layout)
# even_layout = solar_layout
# even_layout.rename('Even Capacity Distribution [MW]')
# even_layout.data = evenDistribution

# # plot new layout
# even_layout.plot(cmap="cool", size=8, aspect=1)
# plt.title("Even Capacity Distribution [MW]")
# plt.tight_layout()

# New layout 2
evenDistribution_GER = solar_layout.data

# Assign new even capacity to only the cells that contain solar capacity
evenDistribution_GER[evenDistribution_GER > 0] = evenCap

even_GER_layout = solar_layout
even_GER_layout.rename('Even Capacity Distribution [MW]')
even_GER_layout.data = evenDistribution_GER

even_GER_layout.plot(cmap="cool", size=8, aspect=1)
plt.title("Even Capacity Distribution within German boarder: " + str(evenCap) + " MW")
plt.tight_layout()

# Compute electricity generation with atlite
pv_even = cutout.pv(panel="CSi", orientation={'slope': 30., 'azimuth': 0.}, layout=even_GER_layout)
pv_even_output = pv_even.to_pandas().iloc[0]

# Determine time step length in hours (h)
dt = pv_even_output.index[1]-pv_even_output.index[0]
dt = dt.days * 24 + dt.seconds // 3600

# Determine total electricity generation 
E_tot = pv_even_output.sum() * dt / 1000 # in GWh


plt.figure()
pv_even_output.resample('1W').mean().plot(figsize=(8,5))
plt.ylabel("Feed-In [GW]")
plt.title('PV time-series Germany with equal capacity distribution (' + str(evenCap) + " MW).\n Total electricity generation E_tot = " + str(round(E_tot,1)) + " GWh")
plt.tight_layout()
