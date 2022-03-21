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

START_YEAR=1993
END_YEAR=2020
#setting up paths monthly averages
csv_outputs = '../../data/csv/transects/'


def plot_transect_timeseries(transect,
                             season, ax, anomaly=True):
    label = SEASONS[season]['label']
    csv_file = os.path.join(csv_outputs, transect, 
                            '{}_{}_CIL_Area.csv'.format(transect, label))
    df = pd.read_csv(csv_file)
    area = df['Area (km^2)']
    if anomaly:
        data = (area - area.mean())/area.std()
    else:
        data = area
    ax.plot(df['Year'], data, label='GLORYS12')
    

def main():
    # Settings for this transect
    season = sys.argv[1]
    label=SEASONS[season]['label']
    anomaly=False
    transects = ['BonneBay', 'Louisbourg', 'HalifaxLine']
    fig, axs = plt.subplots(3, 1, figsize=(8,10), sharex=True)
    for ax, transect in zip(axs, transects):
        plot_transect_timeseries(transect, season, ax, anomaly=anomaly)
        ax.text(0.05, 0.8, transect, transform=ax.transAxes)
        ax.set_xlim([START_YEAR, END_YEAR])
        if anomaly:
            ylabel='Standardized \nAnomaly'
            ax.set_ylim([-3,3])
            fname = os.path.join('figures', 'other_CIL_Area_{}_anomaly.png'.format(label))
        else:
            ylabel='CIL Area (km^2)'
            fname = os.path.join('figures', 'other_CIL_Area_{}.png'.format(label))
        ax.grid()
    ax=axs[0]
    ax.set_ylabel(ylabel)
    ax.legend(bbox_to_anchor=(1.04,1), loc="upper left")
    ax.set_title('CIL Area along \n transects \n {} average'.format(label),
                 fontweight='bold',
                 fontsize=14)
    ax=axs[-1]
    ax.set_xlabel('Year')
    plt.show()
    fig.savefig(fname,bbox_inches='tight')
    

if __name__ == '__main__':
    main()
