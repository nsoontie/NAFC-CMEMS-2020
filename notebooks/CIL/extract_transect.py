# A script to process model data along transects

import glob
import os
import datetime
import geopy.distance
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.interpolate as interp
import sys
import xarray as xr

# Some global constants
TRANSECTS_DICT = {
    'SealIsland': {'sub_domain': [-56,53,-52,56]},
    'WhiteBay' : {'sub_domain': [-57, 50, -49, 53]},
    'FlemishCap': {'sub_domain': [-53, 46.5, -41.5, 47.5]},
    'Bonavista': {'sub_domain': [-53, 48.5, -47.5, 50.5]},
    'BonneBay': {'sub_domain': [-60, 49, -58, 50.5]},
    'HalifaxLine': {'sub_domain': [-64, 42.5, -61, 44.5]},
    'Louisbourg': {'sub_domain': [-60, 43, -57, 46]}
}
START_YEAR=1993
END_YEAR=2020
# Year to switch from GLORYS12 to PSY4
SWITCH_YEAR=2020
SEASONS={'summer': {'months': [6,7,8],
                    'label' : 'JJA'},
         'annual': {'months': [1,2,3,4,5,6,7,8,9,10,11,12],
                    'label': 'annual'},
         'fall': {'months': [9,10,11],
                  'label': 'SON'},
         'winter': {'months': [12, 1,2],
                    'label': 'DJF'},
         'spring': {'months': [3,4,5],
                    'label': 'MAM'}
        }

#setting up paths monthly averages
hank_dir='/home/soontiensn/remote2/hank/'
csv_outputs = '../../data/csv/transects/'

MODEL_PATHS={
    'GLORYS12': {
        'base': os.path.join(hank_dir,'cmems/GLOBAL_REANALYSIS_PHY_001_030'),
        'subdir': 'global-reanalysis-phy-001-030-monthly',
        'filepart': 'mercatorglorys12v1_gl12_mean'
    },
    'PSY4':{
        'base': os.path.join(
            hank_dir,
            'nrt.cmems-du.eu/GLOBAL_ANALYSIS_FORECAST_PHY_001_024'
        ),
        'subdir': 'global-analysis-forecast-phy-001-024-monthly',
        'filepart': 'mercatorpsy4v3r1_gl12_mean'
    },
}

def extract_transect_timeseries(transect,
                                season,
                                start_year,
                                end_year
):
    years = np.arange(start_year, end_year+1)
    months = SEASONS[season]['months']
    label = SEASONS[season]['label']
    print(       os.path.join(
            csv_outputs,
            transect,
            '{}.txt'.format(transect.lower())))
    stations = pd.read_table(
        os.path.join(
            csv_outputs,
            transect,
            '{}.txt'.format(transect.lower())),
        delimiter=' ')
    sub_domain = TRANSECTS_DICT[transect]['sub_domain']
    df=pd.DataFrame()
    for year in years:
        print('Processing year: {}'.format(year))
        if year >= SWITCH_YEAR:
            # Use Psy4
            model='PSY4'
        else: 
            # Use GLORYS
            model='GLORYS12'
        stations = pd.read_table(
            os.path.join(
                csv_outputs,
                transect,
                '{}.txt'.format(transect.lower())),
            delimiter=' ')
        weight_sum=0
        var_avg=0
        if season =='winter' and year==start_year:
            continue
        for month in months:
            if month==12 and season =='winter':
                print('Processing year {} for December winter'.format(year-1))
                distance, depths, var, num_days = \
                    load_model_data(model,
                                    year-1,
                                    month,
                                    csv_outputs,
                                    stations,
                                    sub_domain,
                                    transect)
            else:
                print('Processing month: {}'.format(month))
                distance, depths, var, num_days = \
                load_model_data(model,
                                year,
                                month,
                                csv_outputs,
                                stations,
                                sub_domain,
                                transect)
            var_avg += var*num_days
            weight_sum += num_days
        var_avg = var_avg/weight_sum
        var_avg = np.ma.masked_invalid(var_avg)
        # plot
        fig,ax=plt.subplots(1,1, figsize=(15,5))
        xx,yy,var_regrid = regrid(distance, depths, var_avg)
        plot_transect(xx, yy, var_regrid, ax)
        ax.set_title('{} {} - {} {}'.format(model, transect, label, str(year)))
        out_file = os.path.join('figures', str(year),
                                '{}_{}_{}_{}.png'.format(model, transect, label, str(year)))
        fig.savefig(out_file, bbox_inches='tight')
        # area
        area = calculate_area(xx,yy,var_regrid)
        dnew = pd.DataFrame({'Year': year, 'Model': model,
                             'Area (km^2)': area, 'season': [label]})
        df = pd.concat([df, dnew])
    # Save area dataframe
    df.to_csv(os.path.join(csv_outputs, transect,
                           '{}_{}_CIL_Area.csv'.format(transect, label)))

    
def load_model_data(model,
                    year,
                    month,
                    csv_outputs,
                    stations,
                    sub_domain,
                    transect,
                    variable='Temperature [degrees_C]'):
    """Load model data for given year and season.
    Subset to a sub_domain
    """
    # model directories etc
    basedir = MODEL_PATHS[model]['base']
    subdir = MODEL_PATHS[model]['subdir']
    filepart = MODEL_PATHS[model]['filepart']
    #subset file lookup first
    csv_outdir = os.path.join(csv_outputs, transect, subdir, str(year))
    if not os.path.exists(csv_outdir):
        os.makedirs(csv_outdir)
    # Load netcdf and create file to load later
    csv_outfile = '{}_{}{:02d}_{}.csv'.format(transect,
                                              str(year),
                                              month,
                                              filepart)
    csv = os.path.join(csv_outdir, csv_outfile)
    if not os.path.exists(csv):
        pth=os.path.join(basedir, subdir,
                         str(year), '*{:02d}.nc'.format(month))
        fname = glob.glob(pth)[0]
        print("Loading netcdf file: {}".format(fname))
        d = xr.open_dataset(fname)
        # subset
        sd=sub_domain
        d = d.sel(longitude=slice(sd[0],sd[2]), latitude=slice(sd[1], sd[3]))
        df = interp_to_stations(d, stations)
        df.to_csv(csv,  index=False, na_rep='nan')
        # Load csv
    distance, depths, var, num_days = load_csv(csv, stations)
    return distance, depths, var, num_days


def regrid(xo, yo, tempo, nx=400, ny=400, depthmax=700):
    xnew = np.linspace(xo[0,0],xo[-1,0], num=nx )
    ynew = np.linspace(yo[0,0],depthmax, num=ny )
    yy, xx = np.meshgrid(ynew, xnew)
    tempm = tempo.filled(np.nan)
    tempnew = interp.griddata((xo.flatten(),yo.flatten()),tempm.flatten(), (xx,yy),
                              fill_value=np.nan)
    tempnew = np.ma.masked_invalid(tempnew)
    return xx, yy, tempnew


def plot_transect(x,y,t,ax, contour=True):
    if contour:
        mesh=ax.contourf(x,y,t,cmap='plasma',levels=np.arange(-2,11))
    else:
        below_zero = np.ma.masked_greater(t,0)
        mesh=ax.pcolormesh(x,y,below_zero,cmap='plasma',vmin=-2,vmax=10)
    contour = ax.contour(x, y, t, [0,], colors='k')
    cbar = plt.colorbar(mesh,ax=ax)
    cbar.set_label('Temperature [deg C]')
    ax.set_ylim([500,0])
    ax.set_ylabel('Depth [m]')
    ax.set_xlabel('Distance [km]')
    ax.grid()


def calculate_area(x,y,temp):
    below_zero = np.ma.masked_greater(temp, 0)
    diffs_x = np.diff(x, axis=0)
    diffs_y = np.diff(y, axis=1)
    area = np.nansum(diffs_x[:, :-1]*diffs_y[:-1, :]/1000*(1-below_zero.mask[:-1,:-1]))
    return area


def load_csv(fname, stations, variable='Temperature [degrees_C]'):
    print("Loading csv: {}".format(fname))
    d = pd.read_csv(fname)
    numstations=stations.shape[0]
    # distance along transect in km
    distance = [geopy.distance.distance((stations.iloc[0]['LAT'],
                                         stations.iloc[0]['LON']),
                                        (stations.iloc[i]['LAT'],
                                         stations.iloc[i]['LON'])).km \
                for i in range(numstations)]
    # Create arrays
    var = d[variable].values
    depths = d['Depth [m]'].values
    size=len(var)
    # Reshape data
    numdepths = int(size/numstations)
    var = var.reshape((numstations, numdepths ))
    var = np.ma.masked_invalid(var)
    depths = depths.reshape((numstations, numdepths))
    _, distance = np.meshgrid(depths[0, :], distance)
    num_days = d['num_days'].mean()
    return distance, depths, var, num_days


def interp_to_stations(d, stations, variable='thetao'):
    month_length = d.time.dt.days_in_month
    drop=[]
    for v in d.variables:
        if v not in [variable, 'longitude', 'time', 'latitude', 'depth']:
            drop.append(v)
    print("Interpolating data to stations")
    dnew = pd.DataFrame()
    numrows=stations.shape[0]
    numdepths=d.depth.values.shape[0]
    for i in range(stations.shape[0]):
        print('Station {} of {}'.format(i, stations.shape[0]))
        lat=stations.iloc[i]['LAT']
        lon=stations.iloc[i]['LON']
        sid=stations.iloc[i]['STATION']
        dstation = d.interp(latitude=lat, longitude=lon)
        dstation = dstation.drop(drop)
        df=dstation.to_dataframe()
        df['station'] = sid
        df['num_days'] = month_length.values[0]
        df = df.reset_index()
        dnew = pd.concat([dnew,df])
    # Rename columns
    rename={var: '{} [{}]'.format(d[var].attrs['long_name'],
                                  d[var].attrs['units'])
        for var in dnew.columns if var not in ['station', 'time', 'num_days']}
    rename['station'] = 'STATION'
    rename['time'] = 'TIME'
    rename['longitude'] = 'LON'
    rename['latitude'] = 'LAT'
    rename['num_days'] = 'num_days'
    dnew = dnew.rename(columns=rename)
    dnew = dnew.reset_index(drop=True)
    return dnew


def main():
    # Settings for this transect
    transect = sys.argv[1]
    season = sys.argv[2]
    extract_transect_timeseries(transect, season,
                                START_YEAR,
                                END_YEAR)

if __name__ == '__main__':
    main()
