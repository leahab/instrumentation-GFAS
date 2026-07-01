
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
    ''' 
    calculates the correction factor of the attenuation correction for every individual particle
    fwhm: Series of the full width half maximum (FWHM) of the GFAS pbp data 
    output: Series of correction factors
    '''
    CF = 1 + 0.00404079 * fwhm**1.07336 
    return CF

def attenuation(df):   
    '''
    calculates the attenuation correction of the gfas pbp data for each detector.
    dependencies: function CF 
    df: GFAS PbP as dataframe
    output: corrected GFAS PbP data as full dataset
    '''
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
    ''' 
    Filters the attenuation corrected dataset to only keep values where the signal of Qual/Sizer >= 0.6
    df: attenuation corrected GFAS PbP as dataframe 
    output: filtered GFAS PbP as dataframe 
    '''    
    df = df[df['Qual Peak']/df['Sizer Peak']>= 0.6]
    return df

#%% sizing regression

scale_factors = pd.read_csv('/home/leha7253/instrumentation-GFAS/PbP_processing/detector_scaling_factors', index_col = 0)

def sizing_model(df):
    '''
    calculates particle diameter from forward scattering "Sizer Peak" column of the GFAS PbP dataframe. 
    creates two new columns within the dataframe: "Diameter" and "Scaled Sizer Peak" 
    df: GFAS PbP as dataframe
    ouput: GFAS PbP as dataframe with added columns "Diameter" and "Scaled Sizer Peak"  
     '''
    scale_forward = 1/scale_factors.iloc[0, 0] 
    m = [2.436249909162449, -0.1961175657829281, -0.06641428293828897]
    small_m = [3.54411956337398, 0.643890961346859, 0.02465975028767153] 
    df['scaled Sizer Peak'] = df['Sizer Peak']/scale_forward
    df['Diameter'] = np.where(df['scaled Sizer Peak'] >= 2.37e-08, # 1.86e-08, # 
                            10**(m[0] + m[1]*np.log10(df['scaled Sizer Peak']) + m[2]*np.log10(df['scaled Sizer Peak'])**2), 
                            10**(small_m[0] + small_m[1]*np.log10(df['scaled Sizer Peak']) + small_m[2]*np.log10(df['scaled Sizer Peak'])**2))
    return df

#%% Flags time series for a parameter and threshold value 
def ultimate_flag(vis=None, threshold=1000, min_percentage=0.9, min_duration='30min', plot=False, color='r', alpha=0.5, ax=None):
    '''
    This function flags time intervals in a visibility time series where the visibility is below a specified threshold for a minimum percentage of a specified duration. Can be used for other variables as well. the threshold value just has to be adjustend. It can also plot the flagged intervals. 
    vis: pandas Series of the parameter that is to be flagged (e.g., visibility)
    threshold: float, the threshold value below which the parameter is considered flagged (default: 1000)
    min_percentage: float, the minimum percentage of time below the threshold required to flag the interval (default: 0.9)
    min_duration: str, the minimum duration of the interval to be flagged (default: '30min')
    plot: bool, whether to plot the flagged intervals (default: False)
    color: str, the color of the flagged intervals in the plot (default: 'r')
    alpha: float, the transparency of the flagged intervals in the plot (default: 0.5)
    ax: matplotlib Axes object, if provided, the flagged intervals will be plotted on this axes (default: None)
    output: pandas DataFrame with columns 'Start', 'End', and 'Timespan' indicating the flagged intervals, or a plot of the flagged intervals if plot=True   
    '''
    # Convert the visibility data to float
    param = vis.astype('float')
    
    # Create a binary flag where 1 indicates below the threshold and 0 indicates above
    flag_below_threshold = (param < threshold).astype('int')
    
    # Convert the minimum duration to minutes and calculate the number of required points
    min_duration_minutes = pd.to_timedelta(min_duration).total_seconds() / 60
    required_points = int(min_percentage * min_duration_minutes)
    
    # Use a centered rolling window to calculate the sum of the flag values
    rolling_sum = flag_below_threshold.rolling(window=int(min_duration_minutes), center=True, min_periods=1).sum()
    
    # Identify valid intervals where the rolling sum meets or exceeds the required points
    valid_intervals = rolling_sum >= required_points

    # Find the start and end times of these intervals
    param_flags = pd.DataFrame(columns=['Start', 'End'])
    param_diff = valid_intervals.astype('int').diff().fillna(0)
    start_times = param_diff[param_diff == 1].index
    end_times = param_diff[param_diff == -1].index
    
    # Handle cases where the last interval doesn't have an end time
    if len(end_times) < len(start_times):
        end_times = end_times.append(pd.Index([param.index[-1]]))
    
    param_flags['Start'] = start_times
    param_flags['End'] = end_times
    
    # Filter intervals to ensure they meet the minimum duration
    param_flags['Timespan'] = param_flags['End'] - param_flags['Start']
    param_flags = param_flags[param_flags['Timespan'] >= pd.to_timedelta(min_duration)]
    
    # Merge overlapping intervals if necessary
    merged_intervals = []
    if not param_flags.empty:
        current_start = param_flags.iloc[0]['Start']
        current_end = param_flags.iloc[0]['End']
    
        for i in range(1, len(param_flags)):
            row = param_flags.iloc[i]
            if row['Start'] <= current_end:
                current_end = max(current_end, row['End'])
            else:
                merged_intervals.append([current_start, current_end])
                current_start = row['Start']
                current_end = row['End']
        merged_intervals.append([current_start, current_end])
    
    param_flags = pd.DataFrame(merged_intervals, columns=['Start', 'End'])
    param_flags['Timespan'] = param_flags['End'] - param_flags['Start']
    
    # Plotting
    if plot:
        for _, flag in param_flags.iterrows():
            if ax is None:
                plt.axvspan(pd.to_datetime(flag['Start']),
                            pd.to_datetime(flag['End']),
                            color=color, alpha=alpha)
            else:
                ax.axvspan(pd.to_datetime(flag['Start']),
                            pd.to_datetime(flag['End']),
                            color=color, alpha=alpha)
    else:
        return param_flags