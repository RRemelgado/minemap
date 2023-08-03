# -*- coding: utf-8 -*-
""" detect_change
------------------------------------------------------------------------------
The algorithm uses elbow and knee detection to find the year when a mining 
operation started. First, the algorithm will sort the observations by date 
and derive a cummulative RMSE. This step helps reduce multi-temporal noise 
and highlight abrupt changes in time, which are likely related to start of 
mining operations. To detect this event, the algorithm will use the elbow 
and knee detection method implemented in the kneenow module. Kneebow rotates 
the cummulative RMSE values so that a curve looks down. Then, the elbow will 
correspond t the minimum value of the rotated vector, while the knee will 
correspond to the the maximum. To know whether to find the elbow or the knee, 
the change detection algorithm uses the Area Under the Curve (AUC). If the AUC 
is lower than half of the maximum, the cummulative RMSE time-series becomes 
convex, which is known as an "elbow". Alternatively, if the AUC is greater 
than half of the maximum, the time-series forms a concave curve, know as a 
"knee". The algorithm will report on the likely change year, corresponding 
R2 and number of time-steps used. Note that the algorithm will not do such 
reporting if the minimum R2 and minimum number of time steps don't respect 
the values specified in the configuration yaml.
------------------------------------------------------------------------------
------------------------------------------------------------------------------
Created on Sat Sep 11 12:16:22 2021
@author: Ruben Remelgado
"""

#----------------------------------------------------------------------------#
# load required modules and arguments
#----------------------------------------------------------------------------#

from progress.bar import Bar
from os.path import join, dirname
from sklearn.metrics import r2_score
from argparse import ArgumentParser
from loess.loess_1d import loess_1d
from kneebow.rotor import Rotor
from sklearn.metrics import auc
import pandas as pd
import numpy as np
import fiona as fn
import yaml

parser = ArgumentParser(description='extract landsat SR time series')
parser.add_argument('config', help='path to configuration yaml', type=str)

options = parser.parse_args()
config_file = options.config

# extract working directory
wdir = dirname(config_file)

# access parameters
config = yaml.safe_load(open(config_file, 'r'))

# determine number of features
nr_features = len(fn.open(join(wdir, config['mining_areas'])))

# minimum number of time steps
min_length = config['change_detection']['nr_steps']

# string format to get input file name
sformat = '{0:0' + str(len(str(nr_features))) + 'd}'

#----------------------------------------------------------------------------#
bar = Bar('estimate starting year of each mining site', max=nr_features)
#----------------------------------------------------------------------------#

year = [] # inferred starting year
dates = [] # number of dates used in inferrence
r2 = [] # R2 between time and RMSE
change = [] # NDVI change

for i in range(0,nr_features):
    
    #========================================================================#
    # load and format data
    #========================================================================#
    
    try:
        
        iname = join(wdir, '01_analysis','sr', sformat.format(i) + '_sr.csv')
        reflectances = pd.read_csv(iname)
        reflectances = reflectances.groupby(['timestamp','place']).mean()
        reflectances['timestamp'] = [i[0] for i in reflectances.index]
        reflectances['place'] = [i[1] for i in reflectances.index]
        length_inside = np.sum(reflectances['place'] == 'inside')
        length_outside = np.sum(reflectances['place'] == 'outside')
        
    except:
        
        length_inside = 0
        length_outside = 0
    
    if (length_inside >= min_length) & (length_outside >= min_length):
        
        all_dates = np.unique(reflectances['timestamp'].values)
        all_inside = np.unique(reflectances.loc[
            reflectances['place'] == 'inside','timestamp'])
        
        all_outside = np.unique(reflectances.loc[
            reflectances['place'] == 'outside','timestamp'])
        
        all_dates = all_dates[np.where(np.isin(all_dates, all_inside) & 
                                       np.isin(all_dates, all_outside))]
        
        # register number of time-steps
        dates.append(len(all_dates))
        
        #====================================================================#
        # build NDVI time-series of inside/outside absolute differences
        # proceed to smooth the time series using a loess curve
        #====================================================================#
        
        vii = reflectances['ndvi'].values[
                np.where((reflectances['place'] == 'inside') & 
                np.isin(reflectances['timestamp'], all_dates))]
        
        vio = reflectances['ndvi'].values[
                np.where((reflectances['place'] == 'outside') & 
                np.isin(reflectances['timestamp'], all_dates))]
        
        xo, yo, wo = loess_1d(np.array(range(1,len(vii)+1)), 
                                    np.abs(vii-vio), frac=0.3)
        
        #====================================================================#
        # format data to be used with rotor module
        #====================================================================#
        
        # restructure the data.frame
        # NOTE: needed step to use the knee/elbow detection method
        data = np.zeros((len(xo),2),dtype='float32')
        data[:,0] = xo
        data[:,1] = np.cumsum(yo)
        
        # calculate absolute ndvi change across time-series
        change.append(np.sum(np.diff(yo)))
        
        #====================================================================#
        # estimate starting year
        #====================================================================#
        
        # correlation between time and RMSE
        r2.append(r2_score(data[:,0]/np.max(data[:,0]), data[:,1]/np.max(data[:,1])))
        
        # calculate area under the curve: used to take elbow or knee
        area = auc((data[:,0]/np.max(data[:,0])),(data[:,1]/np.max(data[:,1])))
        
        # fit rotated time-series
        rotor = Rotor()
        rotor.fit_rotate(data)
        
        # estimate elbow (if area < 0.5) or knee (if area > 0.5 )
        if area < 0.5:
            di = rotor.get_elbow_index()
        else:
            di = rotor.get_knee_index()
        
        # estract year of change
        year.append(int(all_dates[di][0:4]))
        
    else:
        
        r2.append(np.nan)
        dates.append(0)
        year.append(config['change_detection']['starting_year'])
        change.append(np.nan)
    
    bar.next()

bar.finish()

#----------------------------------------------------------------------------#
# combine and write output
#----------------------------------------------------------------------------#

odf = pd.DataFrame({'r2':r2, 'change_year':year, 
                'nr_dates':dates, 'change':change})

odf.to_csv(join(wdir, '01_analysis', 'mine_start.csv'), index=False)
