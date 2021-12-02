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
from matplotlib.colors import from_levels_and_colors

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
# For climatology
START_YEAR=1993
END_YEAR=2020
#setting up paths monthly averages
csv_outputs = '../../data/csv/transects/'


def load_transect_timeseries(transect,
                             season):
    label = SEASONS[season]['label']
    csv_file = os.path.join(csv_outputs, transect, 
                            '{}_{}_CIL_Area.csv'.format(transect, label))
    df = pd.read_csv(csv_file)
    area = df['Area (km^2)']
    dclim = df[(df['Year'] >=START_YEAR) & (df['Year'] <= END_YEAR)]
    clim_area = dclim['Area (km^2)']
    df['anomaly'] = (area - clim_area.mean())/clim_area.std()
    df['Clim mean'] = clim_area.mean()
    df['Clim std'] = clim_area.std()
    df['transect'] = transect
    return df
    

def main():
    # Settings for this transect
    season = sys.argv[1]
    label=SEASONS[season]['label']
    transects = ['SealIsland', 'WhiteBay', 'Bonavista', 'FlemishCap']
    dnew = pd.DataFrame()
    for transect in transects:
        df = load_transect_timeseries(transect, season)
        data=list(df['anomaly'].values)
        data.append(df['Clim mean'].values[0])
        data.append(df['Clim std'].values[0])
        cols=list(df['Year'].values)
        cols.append('MEAN')
        cols.append('SD')
        subdf = pd.DataFrame(data,index=cols,columns=[transect,] )
        #subdf['transect']=transect
        #subdf=subdf.set_index('transect')
        dnew = pd.concat([dnew,subdf], axis=1)
    dnew = dnew.transpose()
    dnew =dnew.rename(columns={
    'MEAN' : r'$\rm \overline{x}$',
    'SD': 'sd'})
    # Start the plotting - plotting set up from Fred
    # Build the colormap
    vmin = -3.49
    vmax = 3.49
    midpoint = 0
    levels = np.linspace(vmin, vmax, 15)
    midp = np.mean(np.c_[levels[:-1], levels[1:]], axis=1)
    colvals = np.interp(midp, [vmin, midpoint, vmax], [-1, 0., 1])
    normal = plt.Normalize(-3.49, 3.49)
    reds = plt.cm.Reds(np.linspace(0,1, num=7))
    blues = plt.cm.Blues_r(np.linspace(0,1, num=7))
    whites = [(1,1,1,1)]*2
    colors = np.vstack((blues[0:-1,:], whites, reds[1:,:]))
    colors = np.concatenate([[colors[0,:]], colors, [colors[-1,:]]], 0)
    cmap, norm = from_levels_and_colors(levels, colors, extend='both')
    # Cell parameters
    hcell, wcell = 0.6, 1
    hpad, wpad = 1, .5
    nrows, ncols = dnew.index.size+1, dnew.columns.size 
    fig, ax=plt.subplots(1,1,figsize=(ncols*wcell+wpad, nrows*hcell+hpad))
    ax.axis('off')
    #do the table
    header = ax.table(cellText=[['']],
                      colLabels=['-- NL Shelf transects - {} --'.format(label)],
                      loc='center'
                      )
    header.set_fontsize(13)
    vals_color= -dnew.values #reversed for area
    vals_color [:,-2:] = 0 # zero color for last two columsn (mean/std)
    vals = np.around(dnew.values,1)
    years = np.array(dnew.columns.astype(str))
    for i,yr in enumerate(years[:-2]):
        years[i] =yr[-2:]
    the_table=ax.table(cellText=vals,
                       rowLabels=dnew.index,
                       colLabels=years,
                       loc='center',
                       cellColours=cmap(norm(vals_color)),
                       cellLoc='center',
                       bbox=[0, 0, 1, 0.5]
                    )
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(11)
    plt.show()
    fig.savefig('figures/CIL_scorecard_{}.png'.format(label),
                bbox_inches='tight')

if __name__ == '__main__':
    main()
