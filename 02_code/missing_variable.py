# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 07:40:59 2022

@author: rr70wedu
"""

from argparse import ArgumentParser
from os.path import join, dirname
from loess.loess_1d import loess_1d
import pandas as pd
import numpy as np
import fiona as fn
import yaml


parser = ArgumentParser(description = 'locate savanna occurrences')
parser.add_argument('config', help='configuration file', type=str)

options = parser.parse_args()
config = options.config

# extract base directory
wdir = dirname(config)

# load parameters
config = yaml.safe_load(open(config, "r"))

# determine number of features
nr_features = len(fn.open(join(wdir, config['mining_areas'])))

# minimum number of time steps
min_length = config['change_detection']['nr_steps']

# string format to get input file name
sformat = '{0:0' + str(len(str(nr_features))) + 'd}'

change = np.zeros(nr_features, dtype='float32')

for i in range(0,nr_features):
    
    #========================================================================#
    # load and format data
    #========================================================================#
    
    try:
        reflectances = pd.read_csv(join(wdir, '01_analysis','sr', sformat.format(i) + '_sr.csv'))
        reflectances = reflectances.groupby(['timestamp','place']).mean()
        reflectances['timestamp'] = [i[0] for i in reflectances.index]
        reflectances['place'] = [i[1] for i in reflectances.index]
        length_inside = np.sum(reflectances['place'] == 'inside')
        length_outside = np.sum(reflectances['place'] == 'outside')
    except:
        length_inside = 0
        length_outside = 0
    
    if (length_inside >= min_length) & (length_outside >= min_length):
        
        #====================================================================#
        # find common dates
        #====================================================================#
        
        all_dates = np.unique(reflectances['timestamp'].values)
        all_inside = np.unique(reflectances.loc[
            reflectances['place'] == 'inside','timestamp'])
        
        all_outside = np.unique(reflectances.loc[
            reflectances['place'] == 'outside','timestamp'])
        
        all_dates = all_dates[np.where(np.isin(all_dates, all_inside) & 
                                       np.isin(all_dates, all_outside))]
        
        #====================================================================#
        # smooth ndvi for inside the mine
        #====================================================================#
        
        vii = reflectances['ndvi'].values[np.where((reflectances['place'] == 'inside') & np.isin(reflectances['timestamp'], all_dates))].copy()
        vio = reflectances['ndvi'].values[np.where((reflectances['place'] == 'outside') & np.isin(reflectances['timestamp'], all_dates))].copy()
        
        xo, yo, wo = loess_1d(np.array(range(1,len(vii)+1)), vii-vio, frac=0.3)
        
        change[i] = np.abs(yo[0]-yo[len(yo)-1])
    
    print(i)


odf = pd.DataFrame({'change':change})
odf.to_csv(join(wdir, 'change.csv'), index=False)

