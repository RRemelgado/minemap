# -*- coding: utf-8 -*-
"""extract
------------------------------------------------------------------------------
The algorithm retrieves Landsat multispectral data from Google Earth Engine 
for a set of of samples. These samples are distinguished between those 
"inside" and "outside" a minining area and are use to extract data from 
a list of image collections defined in the required configuration file. 
Before exporting the extracted data, these are aggregated for each timestamp
and place (i.e. "inside"/"outside") through averaging.
------------------------------------------------------------------------------
------------------------------------------------------------------------------
Created on Sun Aug 15 15:48:13 2021
@author: Ruben Remelgado
"""

#----------------------------------------------------------------------------#
# load required modules and arguments
#----------------------------------------------------------------------------#

from progress.bar import Bar
from argparse import ArgumentParser
from os.path import dirname, join
import pandas as pd
import numpy as np
import fiona as fn
import yaml
import ee

parser = ArgumentParser(description='extract landsat SR time series')
parser.add_argument("config", help="path to configuration yaml", type=str)
parser.add_argument("index", help="target file index", type=int)
options = parser.parse_args()
config = options.config
index = options.index

# infer project path from input file
wdir = dirname(config)

# load configuration file
config = yaml.safe_load(open(config, 'r'))

# determine number of features
nr_features = len(fn.open(join(wdir, config['mining_areas'])))

# string format to get input file name
fid = ('{0:0' + str(len(str(nr_features))) + 'd}').format(index)

#----------------------------------------------------------------------------#
# construct masking function used to remove clouds/shadows from imagery
#----------------------------------------------------------------------------#

def mask(image):
    
    # build mask of clouds and cloud shadows
    qa = image.select('QA_PIXEL')
    cloud = qa.bitwiseAnd(1 << 3).Or(qa.bitwiseAnd(1 << 4))
    
    # Remove edge pixels that don't occur in all bands
    edges = image.mask().reduce(ee.Reducer.min())
    
    return image.updateMask(cloud.Not()).updateMask(edges)

#----------------------------------------------------------------------------#
# configure EarthEngine call
#----------------------------------------------------------------------------#

# start API
ee.Initialize()

# load parameters used to search for data
collections = [c[0] for c in config['reflectances']['collections']]
col_bands =  [c[1] for c in config['reflectances']['collections']]
band_names = config['reflectances']['band_names'] + ['timestamp']
sdate = [c[2][0] for c in config['reflectances']['collections']]
edate = [c[2][1] for c in config['reflectances']['collections']]

#----------------------------------------------------------------------------#
# load data and extract random subset to respect GEE's query limits
#----------------------------------------------------------------------------#

# read samples
samples = pd.read_csv(join(wdir, '01_analysis/xy/' + fid + '_xy.csv'))

# exract samples from inside the mining site
px = np.where(samples['place'] == 'inside')[0]

ns = 50
if len(px) < 50: 
    ns = len(px)

np.random.shuffle(px)
xy_inside = ee.Geometry.MultiPoint(coords=[
    ee.Geometry.Point(samples['x'].values[px[i]], 
                      samples['y'].values[px[i]]) for i in range(0,ns)])

# exract samples from outside the mining site
px = np.where(samples['place'] == 'outside')[0]

ns = 50
if len(px) < 50: 
    ns = len(px)

np.random.shuffle(px)
xy_outside = ee.Geometry.MultiPoint(coords=[
    ee.Geometry.Point(samples['x'].values[px[i]], 
                      samples['y'].values[px[i]]) for i in range(0,ns)])

#----------------------------------------------------------------------------#
bar = Bar('extract surface reflectances to samples', max=len(collections))
#----------------------------------------------------------------------------#

output_df = [] # file where to write samples (will be concatenated)

# iterate through each collection
for c in range(0,len(collections)):
    
    col = collections[c]
    bands = col_bands[c] + ['timestamp']
    start = str(sdate[c]) + '-01-01'
    end = str(edate[c]) + '-12-31'
    
    try:
        
        # access landsat collection (subset to desired date range)
        landsat = ee.ImageCollection(col).filterDate(start, end).map(mask)
        
        # extract data from inside the mine
        data_1 = landsat.getRegion(xy_inside, 30).getInfo()
        data_1 = pd.DataFrame(data_1[1:], columns=data_1[0])
        data_1['timestamp'] = pd.to_datetime(data_1['time'], unit='ms')
        data_1 = data_1.dropna()
        data_1 = data_1[bands]
        data_1.columns = band_names
        data_1['ndvi'] = (data_1['nir']-data_1['red']) / (data_1['nir']+data_1['red'])
        data_1['place'] = 'inside'
        
        # extract data from outside the mine
        data_2 = landsat.getRegion(xy_outside, 30).getInfo()
        data_2 = pd.DataFrame(data_2[1:], columns=data_2[0])
        data_2['timestamp'] = pd.to_datetime(data_2['time'], unit='ms')
        data_2 = data_2.dropna()
        data_2 = data_2[bands]
        data_2.columns = band_names
        data_2['ndvi'] = (data_2['nir']-data_2['red']) / (data_2['nir']+data_2['red'])
        data_2['place'] = 'outside'
        
        # extract data from GEE for given dates and points
        output_df.append(pd.concat((data_1, data_2)))
    
    except:
        
        pass
    
    bar.next()

bar.finish()

#----------------------------------------------------------------------------#
# summarize and write extracted data
#----------------------------------------------------------------------------#

# concatenate extracted data and derive mean values
# (derives means per each timestamp, both "inside" and "outside" the polygon)
output_df = pd.concat(output_df)

if output_df.shape[0] > 0:    
    output_df = output_df.groupby(['timestamp','place']).mean()
    output_df['timestamp'] = [i[0] for i in output_df.index]
    output_df['place'] = [i[1] for i in output_df.index]

# write output
output_df.to_csv(join(wdir, '01_analysis', 'sr', fid + '_sr.csv'), index=False)
