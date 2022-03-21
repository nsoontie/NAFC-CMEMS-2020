# A script to demonstrate plotting GLORYS12 surface fields.
# Nancy Soontiens Feb 2022

# package imports
import cartopy.crs as ccrs # a package for plotting maps
import matplotlib.pyplot as plt # a package for creating plots
import xarray as xr # a package for opening netCDF files

# Define path to a GLORYS12 file
data_file='/data/cmems/my.cmems-du.eu/GLOBAL_REANALYSIS_PHY_001_030/global-reanalysis-phy-001-030-daily/2018/07/mercatorglorys12v1_gl12_mean_20180715_R20180718.nc'

# Define region for plot
lon_min=-70
lat_min=40
lon_max=-30
lat_max=70

# Load data file and truncate to area of interest
d = xr.open_dataset(data_file)
dsel = d.sel(latitude=slice(lat_min, lat_max),
             longitude=slice(lon_min, lon_max))
# Extract variables of interest
thetao_surface = dsel['thetao'].isel(depth=0)
lon = dsel['longitude']
lat = dsel['latitude']
thetao_surface.plot()
plt.show()

# More detailed plotting with cartopy (include coastlines, gridlines, projections etc) 
# Create a figure and axis object
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
# Create the pcolormesh
mesh = ax.pcolormesh(lon,
                     lat,
                     thetao_surface.squeeze(),
                     transform=ccrs.PlateCarree(),
                     vmin=-2,
                     vmax=25,
                     cmap='plasma')
# Add colorbar
cbar = plt.colorbar(mesh, ax=ax)
cbar.set_label('{} [{}]'.format(thetao_surface.long_name,
                                thetao_surface.units))
# Add coastlines
ax.coastlines()
# Add gridlines 
gl = ax.gridlines(draw_labels=True)
gl.xlabels_top = False
gl.ylabels_right = False
plt.show()

# More cartopy resources
# https://scitools.org.uk/cartopy/docs/v0.13/matplotlib/advanced_plotting.html
# https://scitools.org.uk/cartopy/docs/v0.13/matplotlib/gridliner.html
