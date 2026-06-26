
'''
This is the sub-routine for GFAS correction processes including:
- attenuation correction
- qualifier/sizer filter
- scale factor input
- sizing regression

It uses the following datasets: 
- scale factors in PbP_processing/detector_scaling_factors

author: Lea Haberstock, Stockholm University, Department of Environmental Science, Atmospheric Unit 
developed toether with: Darrel Baumgardner and Paul Zieger

contact: lea.haberstock@aces.su.se 
used in publication: Haberstock et al. 2026 (submitted to AMT) 


Last modified June 26 2026

'''
#%% GFAS processing function script
# load packages
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import sys  
#import Data_excluded
import numpy as np
import pandas as pd
import glob2
import glob, os
from io import BytesIO
import datetime as dt
import calendar
from pandas.errors import EmptyDataError
from scipy.ndimage import shift
from zipfile import ZipFile
import dask
import dask.dataframe as dd

#%% Attenuation correction
def CF(fwhm):
    CF = 1 + 0.00404079 * fwhm**1.07336 
    return CF
def attenuation(df):    
    cf_sizer = CF(df['Sizer FWHM']) # calculate correction factor for every particle
    cf_qual = CF(df['Qual FWHM'])
    cf_s = CF(df['S FWHM'])
    cf_p = CF(df['P FWHM'])

    df['Sizer Peak'] = df['Sizer Peak']*cf_sizer # apply correction factor to the peak values
    df['Qual Peak'] = df['Qual Peak']*cf_qual
    df['S Peak'] = df['S Peak']*cf_s
    df['P Peak'] = df['P Peak']*cf_p
    return df

#%% remove any particle where Sizer Peak and Qual Peak have a ratio of Qual/Sizer < 0.6
def qual_filter(df):
    df = df[df['Qual Peak']/df['Sizer Peak']>= 0.6]
    return df

#%% sizing regression

scale_factors = pd.read_csv('/home/leha7253/instrumentation-GFAS/PbP_processing/detector_scaling_factors', index_col = 0)

def sizing_model(df):
    scale_forward = 1/scale_factors.iloc[0, 0] #1.8*10**10 #1.96*10**10#
    m = [2.436249909162449, -0.1961175657829281, -0.06641428293828897]
    small_m = [3.54411956337398, 0.643890961346859, 0.02465975028767153] # 3.54411956337398,#3.55911956337398, 
    df['scaled Sizer Peak'] = df['Sizer Peak']/scale_forward
    df['Diameter'] = np.where(df['scaled Sizer Peak'] >= 2.37e-08, # 1.86e-08, # 
                            10**(m[0] + m[1]*np.log10(df['scaled Sizer Peak']) + m[2]*np.log10(df['scaled Sizer Peak'])**2), 
                            10**(small_m[0] + small_m[1]*np.log10(df['scaled Sizer Peak']) + small_m[2]*np.log10(df['scaled Sizer Peak'])**2))
    return df


