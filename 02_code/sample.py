# -*- coding: utf-8 -*-
"""  sample
------------------------------------------------------------------------------
The algorithm collects samples on mine occurrences based on pre-computed data 
on mine extents. For every polygon in a reference vector dataset specified in 
a yaml configuration file, the algorithm will first rasterize it. The pixel 
resolution of the rasterized polygon is specified in the yaml. Then, the 
algorithm will find all pixels inside the polygon and return their spatial 
coordinates. Then, through covolution, the algorithm will find the pixels 
that are along the borders of the reference polygon, returning their spatial 
coordinates. The output is a set a data.frame, writen as a CSV file, reporting 
on the coordinates of each selected pixel and on their respective "place" 
in relation to the polygon (i.e. "inside" or "outside" the polygon). The 
samples for each polygon are writen in a directory pre-defined in the 
configuration yaml file, and the corresponding sample file is named according 
to the processing order (e.g. polygon 1 in a vector file with 21060 entries 
leads to a sample file named "00001.csv").
------------------------------------------------------------------------------
------------------------------------------------------------------------------
Created on Sun Aug 15 17:18:32 2021
@author: Ruben Remelgado
"""

#----------------------------------------------------------------------------#
# load modules and required arguments
#----------------------------------------------------------------------------#

from progress.bar import Bar
from argparse import ArgumentParser
from rasterio.features import geometry_mask
from scipy.ndimage import convolve
from shapely.geometry import shape
import rasterio as rt
import pandas as pd
import numpy as np
import fiona as fn
import yaml
import os

parser = ArgumentParser(description = 'Collect samples')
parser.add_argument('config', 'path to configuration yaml')

options = parser.parse_args()
config = options.config
wdir = os.path.dirname(config)
config = yaml.safe_load(open(config, 'r'))

# read polygons with mining areas
polygons = fn.open(os.path.join(wdir, config['mining_areas']))

# pixel resolution
pr = config['pixel_resolution']

#----------------------------------------------------------------------------#
nr_mines = len(polygons)
bar = Bar('sample mining areas', max=nr_mines)
#----------------------------------------------------------------------------#

for i in range(0,nr_mines):
    
    # translate geometry of polygon into shape (needed to extract bounds)
    s = shape(polygons[i]['geometry']).buffer(config['pixel_resolution']*2)
    
    #========================================================================#
    # rasterize polygon
    #========================================================================#
    
    # spatial transform (i.e. starting x/y and pixel resolution)
    t = rt.transform.from_origin(s.bounds[0]-pr, s.bounds[3]+pr,pr,pr)
    
    # array dimensions
    nr = int((s.bounds[3]-s.bounds[1])/pr)+1 # number of rows
    nc = int((s.bounds[2]-s.bounds[0])/pr)+1 # number of columns
    
    # rasterization
    ia = geometry_mask(geometries=[polygons[i]['geometry']], \
                       out_shape=(nr,nc), transform=t, \
                           invert=True, all_touched=True)
    
    #========================================================================#
    # extract samples
    #========================================================================#
    
    x_coord = [] # x coordinates
    y_coord = [] # y coordinates
    place = [] # location of samples (i.e. "inside" or "outside" polygon)
    
    # sample inside polygon
    px = np.where(ia == True)
    xy = rt.transform.xy(t, px[0],px[1])
    x_coord = x_coord + xy[0]
    y_coord = y_coord + xy[1]
    place = place + ['inside']*len(xy[0])
    
    # sample outside polygon (only along the border)
    ma = convolve(ia.astype('float32'), weights=np.ones((3,3))) # find borders
    px = np.where((ia == False) & (ma > 0))
    xy = rt.transform.xy(t, px[0],px[1])
    x_coord = x_coord + xy[0]
    y_coord = y_coord + xy[1]
    place = place + ['outside']*len(xy[0])
    
    #========================================================================#
    # combine and write samples
    #========================================================================#
    
    # combine samples into one dta.frame
    odf = pd.DataFrame({'x':x_coord, 'y':y_coord, 'place':place})
    odf['id'] = i # add polygon unique identifier
    
    # write samples
    oname = os.path.join(wdir, config['directories']['coordinates'] + '{0:05d}'.format(i) + '.csv')
    odf.to_csv(oname, index=False)
    
    bar.next()

bar.finish()

