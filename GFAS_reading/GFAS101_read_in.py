
'''
This is the sub-routine for GFAS read-in

author: Lea Haberstock, Stockholm University, Department of Environmental Science, Atmospheric Unit 
developed toether with: Darrel Baumgardner and Paul Zieger

contact: lea.haberstock@aces.su.se 
used in publication: Haberstock et al. 2026 (submitted to AMT) 


Last modified June 26 2026
'''

#%% Read in functions for the GFAS101 version (manufactured after 2023, with only one gain)

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

#%% read in pbp data
# check if given string is a date format
def is_date_format(name, datefmt):
    '''checks if given string is actually a date
    name: given string
    datefmt: the format the date is meant to be in 
    '''
    try:
        dt.datetime.strptime(name, datefmt)
        return True
    except ValueError:
        return False
    

# this function is used for selecting time periods within the folder structure of the GFAS
def timeselect(f,date,datefmt,starttime,endtime,timefmt):    
    ''' this function is used for selecting time periods within the folder structure of the GFAS
    f: directory
    date: part of the directory that is the date
    datefmt: the format the date is meant to be in 
    starttime: beginning of timespan in fmt datefmt
    endtime: end of timespan in fmt datefmt
    timefmt: the format the date is meant to be in 
    '''        
    start = calendar.timegm(dt.datetime.strptime(starttime, timefmt).date().timetuple())
    end = calendar.timegm(dt.datetime.strptime(endtime, timefmt).date().timetuple())

    ctime=calendar.timegm(dt.datetime.strptime(date, datefmt).date().timetuple())

    return start<=ctime and end>=ctime

# read in for zipped files
def read_csv_from_zip(zip_path, filetype):
    ''' read in routine if folders are zipped
    zip_path: filepath
    filetype: targeted GFAS filetype: 'GFAS_PbP_Data' or  'GFAS_User_Data'
    '''
    try:
        with ZipFile(zip_path) as zip_file:
            for text_file in zip_file.infolist():
                if text_file.filename.endswith('.csv') and text_file.filename.startswith(filetype):
                    df = pd.read_csv(BytesIO(zip_file.open(text_file.filename).read()), encoding='latin1')
                    df.rename(columns=lambda x: x[0:6] if x.startswith('Bin') else x.strip(), inplace=True)
                    df['datetime'] = (pd.to_datetime(dt.datetime(1904, 1, 1)) + 
                                      pd.to_timedelta(df['Computer Time (sec)'].values, 's'))
                    
                    return df
    except pd.errors.EmptyDataError:
        print(f"No columns to parse from file {zip_path}")
    return None

#read in for unzipped files
def read_csv_file(file_path, filetype):
    try:
        df = pd.read_csv(file_path, encoding='latin1')
        df.rename(columns=lambda x: x[0:6] if x.startswith('Bin') else x.strip(), inplace=True)
        df['datetime'] = (pd.to_datetime(dt.datetime(1904, 1, 1)) + 
                          pd.to_timedelta(df['Computer Time (sec)'].values, 's'))
       
        return df
    except pd.errors.EmptyDataError:
        print(f"No columns to parse from file {file_path}")
    return None

#this function calls the functions above and reads in the GFAS data either as a dask product or as a dataframe
def readPbP101(path, starttime, endtime, timefmt,
               filetype='GFAS_PbP_Data', compute=True):
    '''
    this function calls the functions above and reads in the GFAS data either as a dask product or as a dataframe
    path: path in which all the different GFAS data are stored. This path should lead to the folders that are called by the day of the sampling
    starttime: beginning of timespan in fmt datefmt
    endtime: end of timespan in fmt datefmt
    timefmt: the format the date is meant to be in
    filetype: targeted GFAS filetype: 'GFAS_PbP_Data' or  'GFAS_User_Data'
    compute: if True, computes the data to pandas dataframe. If False, keeps it as a dask thingie
    '''

    if path is None:
        raise ValueError("Path to data folder must be provided")

    # ----------------------------
    # find folders in time window
    # ----------------------------
    flist = [
        f for f in glob.glob(os.path.join(path, "*"))
        if timeselect(f, f[-8:], '%Y%m%d', starttime, endtime, timefmt)
    ]

    if not flist:
        print("nothing found in flist")
        return pd.DataFrame() if compute else None

    # ----------------------------
    # collect files
    # ----------------------------
    zip_files = [
        os.path.join(folder, f)
        for folder in flist
        for f in glob.glob(os.path.join(folder, "*.zip"))
    ]

    csv_files = [
        csv for folder in flist
        for csv in glob.glob(os.path.join(folder, "**", "*.csv"), recursive=True)
        if filetype in csv
    ]

    all_files_exist = bool(csv_files or zip_files)
    if not all_files_exist:
        return pd.DataFrame() if compute else None

    # ----------------------------
    # infer schema from one file
    # ----------------------------
    sample = None

    for f in csv_files:
        sample = read_csv_file(f, filetype)
        if sample is not None and not sample.empty:
            break

    if sample is None:
        for z in zip_files:
            sample = read_csv_from_zip(z, filetype)
            if sample is not None and not sample.empty:
                break

    if sample is None:
        return pd.DataFrame() if compute else None

    meta = sample.iloc[0:0]  # empty dataframe with correct schema

    
    # ----------------------------
    # safe reader (schema enforced)
    # ----------------------------
    def safe_read(func, filepath, filetype):
        df = func(filepath, filetype)
        if df is None or df.empty:
            return meta.copy()
        return df.reindex(columns=meta.columns)

    # ----------------------------
    # delayed loading
    # ----------------------------
    delayed_csv = [
        dask.delayed(safe_read)(read_csv_file, f, filetype)
        for f in csv_files
    ]

    delayed_zip = [
        dask.delayed(safe_read)(read_csv_from_zip, z, filetype)
        for z in zip_files
    ]

    delayed_all = delayed_csv + delayed_zip

    print(f"Scheduled {len(csv_files)} CSV and {len(zip_files)} ZIP files.")

    # ----------------------------
    # build dask dataframe
    # ----------------------------
    ddf = dd.from_delayed(delayed_all, meta=meta)

    if not compute:
        return ddf

    # ----------------------------
    # compute to pandas
    # ----------------------------
    GFAS = ddf.compute()
    GFAS.reset_index(drop=True, inplace=True)
    return GFAS

#%% read in 1Hz data
def readGFAS101(path=None, starttime=None, endtime=None, timefmt=None, filetype='GFAS_User_Data'):
    '''
    this function reads in the GFAS 1Hz data either as a dataframe
    path: path in which all the different GFAS data are stored. This path should lead to the folders that are called by the day of the sampling
    starttime: beginning of timespan in fmt datefmt
    endtime: end of timespan in fmt datefmt
    timefmt: the format the date is meant to be in
    filetype: targeted GFAS filetype: should be 'GFAS_User_Data'
    '''

    if (path is None) or (starttime is None) or (endtime is None) or (timefmt is None):
        print('Path, starttime, endtime, and timefmt as input needed...')
        #return None

    if filetype not in ['GFAS_User_Data', 'GFAS_Diagnostic_Data']:
        print('Filetype has to be GFAS_User_Data or GFAS_Diagnostic_Data')
        #return None

    flist = [f for f in glob.glob(path + "/*") if timeselect(f, f[-8:], '%Y%m%d', starttime, endtime, timefmt)]

    # Collect information about files that failed to load
    failed_files = []

    # Process zipped folders
    dfs = []
    flists2 = []
    for file in flist:
        flist2 = [f for f in glob.glob(file + "/*") if f.endswith('.zip')]
        flists2.extend(flist2)

    for file in flists2:
        try:
            zip_file = ZipFile(file)
            for text_file in zip_file.infolist():
                if text_file.filename.endswith('.csv') and text_file.filename.startswith(filetype):
                    df = pd.read_csv(BytesIO(zip_file.open(text_file.filename).read()))
                    df.rename(columns=lambda x: x[:6] if x.startswith('Bin') else x, inplace=True)
                    df.rename(columns=lambda x: x.strip(), inplace=True)
                    df.set_index('Time Stamp (UTC sec)', inplace=True)
                    df['datetime'] = pd.to_datetime(df.index, unit='s', origin='1904-01-01')
                    df['Computer Time (sec)'] = pd.to_datetime(df['Computer Time (sec)'], unit='s', origin='1904-01-01')
                    df.rename(columns={'Computer Time (sec)': 'Computer Time (UTC)'}, inplace=True)
                    df.set_index('datetime', inplace=True)
                    dfs.append(df)
        except Exception as e:
            print(f"Error reading file {file}: {e}")
            failed_files.append(file)

    GFAS_zipped = pd.concat(dfs) if dfs else pd.DataFrame()

    # Process non-zipped folders
    flists3 = []
    for file in flist:
        flist3 = [f for f in glob.glob(file + "/*") if (f[-4:] != '.zip' and f[-4:] != '.log')]
        flists3.extend(flist3)

    flists4, flists5 = [], []

    for file in range(len(flists3)):
        #User data:
        flist4 = [f for f in glob.glob(flists3[file] + "/*") if (f.endswith('.csv') and f[-38:-24]==filetype)]
        flists4=flists4+flist4
        #Diagnostic data:
        flist5 = [f for f in glob.glob(flists3[file] + "/*") if (f.endswith('.csv') and f[-44:-24]==filetype)]
        flists5=flists5+flist5


    li = []
    for file in flists4 if filetype == 'GFAS_User_Data' else flists5:
        try:
            df3 = pd.read_csv(file, encoding='iso-8859-1')
            df3.rename(columns=lambda x: x[:6] if x.startswith('Bin') else x, inplace=True)
            df3.rename(columns=lambda x: x.strip(), inplace=True)
            df3.set_index('Time Stamp (UTC sec)', inplace=True)
            df3['datetime'] = pd.to_datetime(df3.index.values, unit='s', origin='1904-01-01')
            df3['Computer Time (sec)'] = pd.to_datetime(df3['Computer Time (sec)'], unit='s', origin='1904-01-01')
            df3.rename(columns={'Computer Time (sec)': 'Computer Time (UTC)'}, inplace=True)
            df3.set_index('datetime', inplace=True)
            li.append(df3)
        except Exception as e:
            print(f"Error reading file {file}: {e}")
            failed_files.append(file)

    GFAS_unzipped = pd.concat(li) if li else pd.DataFrame()

    GFAS = pd.concat([GFAS_zipped, GFAS_unzipped], axis=0).sort_index() if not GFAS_zipped.empty or not GFAS_unzipped.empty else pd.DataFrame()

    if failed_files:
        print(f"The following files could not be read:\n{failed_files}")

    if filetype == 'GFAS_User_Data' and not GFAS.empty:
        bins = np.array([0.5, 0.7, 0.9, 1.1, 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18.,
                            19., 20., 21., 22., 23., 24., 25., 26., 27., 28., 29., 30., 31., 32., 33., 34., 35., 38., 40.]) * 10**(-6)
        for i in range(1, 41):
            GFAS['%.9f' % bins[int(i-1)]] = GFAS[f'Bin {int(i)}'] / GFAS['Sample Volume Flow Rate (cm^3/s)']

        Dpmin = 0.4499999880791
        dlogDp = np.log10(GFAS.columns[-40:].astype(float).values) - np.roll(np.log10(GFAS.columns[-40:].astype(float).values), 1)
        dlogDp[0] = np.log10(GFAS.columns[-40:].astype(float).values)[0] - np.log10(Dpmin * 10**(-6))
        GFAS.iloc[:, -40:] = GFAS.iloc[:, -40:] / dlogDp

    return GFAS
