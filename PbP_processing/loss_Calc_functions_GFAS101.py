import sys  
import numpy as np
import pandas as pd
import xarray as xr
import glob, os
import datetime as dt
import calendar
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib import pyplot
from pandas.errors import EmptyDataError
from scipy.ndimage.interpolation import shift
import metpy
from metpy.calc import wind_components,wind_speed,wind_direction
from metpy.units import units

# calculate angle between wind direction and sampling direction
def angle_diff(a,b):
    ''' 
    calculates the steep angle between two wind directions in degree
    a: angle 1
    b: angle 2
    returns: steep angle between the two angles
    '''
    head_diff = a-b
    head_diff = head_diff.where(head_diff < 180, (360+b)-a)
    head_diff = head_diff.where(head_diff > -180, (360+a)-b)
    head_diff = head_diff.where(head_diff > 0, b-a)
    return head_diff

############### function for calculating sample efficiencies ################
def loss_calc_GFAS101(GFAS_df, Wind, wind_GFAS = True, min_efficiency=0.01, return_sep = False):
    ''' GFAS_df: pandas dataframe with GFAS data, index has to be called 'DateTime', 
                 data should only be dN/dlogDp and the columns should be the diameter bins (in m), 
                 and the following columns: 'Sample Velocity (m/s)','GFAS Heading (deg)',('Wind Speed (m/s)','Wind Direction (deg)')
        wind_GFAS: if True, uses the wind speed and direction from the GFAS data, if False, uses the wind speed and direction from the Wind dataframe
        Wind: pandas dataframe with ambient wind speed and direction data in it must be >= 1 min resolution. needs to be called df['Wind Speed (m/s)']  needs to be called df['Wind Direction (deg)'] or df['GFAS Heading (deg)']
        return_sep: if True, returns the separate efficiencies (nasp, n_trm, n_tsp) as xarray DataArrays

        returns: total efficiency as xarray DataArray

        needs: angle_diff, xarray, numpy, pandas
    '''
    
    gfas = GFAS_df.resample('1min').mean() # make sure that the data is in 1 min resolution

    if wind_GFAS == True: # if no wind speed data is given, use the wind speed from the GFAS data
        GFAS_df['u'], GFAS_df['v']= wind_components(GFAS_df['Wind Speed (m/s)'].values*units('m/s'), GFAS_df['Wind Direction (deg)'].values*units.deg)
        gfas.loc[:,['u','v']] = GFAS_df.loc[:,['u','v']] .resample('1min').mean()
        gfas['Wind Direction (deg)'] = wind_direction(gfas['u'].values* units("m/s"), gfas['v'].values* units("m/s"), convention='from')
        WindDirection = gfas['Wind Direction (deg)']
    else:
        U, V = wind_components(Wind['Wind Speed (m/s)'].values*units('m/s'), Wind['Wind Direction (deg)'].values*units.deg)
        U = U.resample('1min').mean()
        V = V.resample('1min').mean()
        WindDirection = wind_direction(U, V)

    WindSpeed = Wind['Wind Speed (m/s)'].resample('1min').mean()
    idx = np.intersect1d(gfas.index, WindSpeed.index)
    gfas = gfas.loc[idx]
    WindSpeed = WindSpeed.loc[idx] # make sure that the wind speed is in the same time frame as the gfas data

    WindDirection = WindDirection.loc[idx] # calculate wind direction from U and V components


    ##field data
    U0 =  xr.DataArray(WindSpeed) # wind speed in m/s
    index = U0.DateTime

    def is_number(x):
        try:
            float(x)
            return True
        except:
            return False

    numeric_name_cols = [c for c in gfas.columns if is_number(c)] 
    D = gfas[numeric_name_cols].columns.astype(float).values # diameter in m (dataframe of all different diameter bins)
    # see if any value in D is larger than 0.0001, if so, convert to m (since some GFAS data is in µm)
    if (D > 0.1).any():
        D = D/10**6 # convert to m if in µm

    TAS = xr.DataArray(gfas['Sample Velocity (m/s)'], dims = 'DateTime', coords={'DateTime':index}) # true air speed; mean value of flow velocity in the wind tunnel im m/s ; laser cross section area is 0.22 mm^2

    theta_s = xr.DataArray(angle_diff(gfas['GFAS Heading (deg)'],WindDirection))
    theta_sr = np.deg2rad(theta_s)      #sampling angle in rad

        ## assumptions atmosphere (Constants)
    rho_g = 1.20 # air density[kg/m^3]; (at 293 K, Hinds, W.C., Aerosol Technology, Wiley, 1999)
    rho_p =   997   # paricle density (kg/m3, water as assumption)
    g = 9.81 # gravitation constant (m/(s^2))
    phi = np.radians(90) # zenith angle of inlet in rad 
    lamda = 66*10**(-9)  # mean free path (m) in air usually between 64 - 68 *10^(-9) m (Hinds, W.C., Aerosol Technology, Wiley, 1999)
    mu = 1.81*10**(-5) # dynamic viscosity of air (Pa*s) (Hinds, W.C., Aerosol Technology, Wiley, 1999)

        ## device specific parameters (GFAS 3)!!! ask for changes 
    # Constants 
    Lw = 0.14354 # length of wind tunnel in m
    d0 = 0.02504 # diameter of wind tunnel in m
    di = 0.05007   # diameter of inlet in m
    Li = 0.12005  # length of inlet in m
    theta = 0  # angle of inlet inclination (0° for horizontal flows) 
    theta_cont = np.arctan(((di-d0)/2)/Li) #contraction, half angle of contraction part in rad; (5.95 degrees)
    A0 = np.pi*(0.5*d0)**2 # cross section of area of wind tunnel in m^2
    Ai = np.pi*(0.5*di)**2# cross section area of the inlet in m^2  

        ## sampling Parameters
    U = (d0**2)*TAS/(di**2)  # sampling speed in m/s 
    Rv = U0/U
    Re = rho_g*U*di/mu  # reynoldsnumber

    # size dependent parameters
    Cc = 1+lamda/D*(2.34+1.05*np.exp(-0.39*D/lamda)) #dimensionless Cunningham correction factor (Hinds, W.C., Aerosol Technology, Wiley, 1999)
    rel_t = xr.DataArray(rho_p*D**2*Cc/(18*mu), dims = 'D', coords = {'D': D}) # relaxation time in s (Hinds, W.C., Aerosol Technology, Wiley, 1999)
    V_ts =xr.DataArray(g*rel_t, dims='D', coords= {'D': D})# terminal settling velocity in m/s (Hinds, W.C., Aerosol Technology, Wiley, 1999)

        ## dimensionless parameters dependent on wind speed
    Stk = xr.DataArray((U0*rel_t)/d0, dims = ['DateTime','D',], coords = {'DateTime' : index, 'D': D}) # stokes number (Hinds, W.C., Aerosol Technology, Wiley, 1999)
    Stk_ = xr.DataArray(Stk*np.exp(0.022*np.rad2deg(theta_sr)), dims = ['DateTime', 'D' ], coords = {'DateTime' : index, 'D': D}) # corrected stokes number for measurement angel (Spiegel et al., 2012)

    ########## Aspiration efficiency calculation

        ### moving reigme U0 > 2.1 m/s (Hangan and Willeke, 1990)
        # for theta_s 0-60 degrees and stokes number 0.01 < Stk < 6, 0.5 < Rv < 2
    f1 = (1 - (1 + Stk_*(2+0.617/Rv))**(-1))/(1 - (1 + 2.617*Stk_)**(-1)) # dim = U0xDxtheta_s
    f2 = 1-(1 + 0.55*Stk_*np.exp(0.25*Stk_))**(-1)# dim = U0xDxtheta_s
    f= f1*f2
    itheta_sr = theta_sr.where(theta_sr < np.deg2rad(60), float('nan'))
    nasp_move = 1+(Rv*np.cos(itheta_sr)-1)*f

        # for theta_s 60-90 degrees and 0.003 < Stk < 0.2, 0.5 < Rv < 2
    iStk = Stk.where((Stk <= 0.2), float('nan')) #&(Stk >= 0.003)
    itheta_sr = theta_sr.where((theta_sr > np.deg2rad(60))&(theta_sr < np.deg2rad(90)), float('nan'))
    nasp_move1 = 1 + (Rv*np.cos(itheta_sr) - 1)*(3*iStk**(Rv**(-0.5)))


    # merge datasets for different theta_s
    nasp_move = nasp_move.where(theta_sr < np.deg2rad(60), nasp_move1)

            ### calm regime Rv < 0.5 (Grinshpun et al., 1993)
    iStk = Stk.where(Stk <= 100) # dim: U0xD
    nasp_calm = np.exp( -(4*iStk**( 1+np.sqrt(V_ts/U) )/(1+2*iStk)) ) + (V_ts/U)*np.cos(phi) # dim: U0xD # only needed for sampling that's not horizontal: + VU*np.cos(phi) 

        #For Vts/U > 1 we set ηasp(D) = NaN
    nasp_calm = nasp_calm.where(V_ts/U < 1, float('nan'))

        ### slow moving regime 0.5 < U0 < 2.1 (Grinshpun et al., 1993)

    VU0 = (V_ts/U0).T  #Grinshpun et al. (1993, 1994)   
        # correction factors for slow moving regime
    f_move = np.exp(-VU0) # dim = U0xD
    f_calm = 1-np.exp(-VU0) # dim = U0xD
    d = xr.DataArray(2*np.cos(theta_sr+phi)) # dim = theta_s
    delta=xr.DataArray((VU0+d)*VU0, dims = ['DateTime', 'D'], coords = {'DateTime' : index, 'D': D}) 

        # apply correction factors
    calm = (nasp_calm*f_calm) 
    move = (nasp_move*f_move*(1+delta)**0.5)

    slow = calm+move #np.where(np.isnan(calm), move, calm + np.nan_to_num(move)) #add calm and move with awareness to nan data 
    slow = xr.DataArray(slow,dims = ['DateTime', 'D'], coords = {'DateTime' : index, 'D': D}) #np.where(slow == 0, float('nan'),slow)

        # bring together all air regimes into one data array for aspiration efficiency nasp
    nasp = xr.DataArray( dims = ['DateTime', 'D'], coords = {'DateTime' : index, 'D': D})

    # for 10−3 ≤ Vts/U ≤ 1 and 0.001 ≤ Stk ≤ 100
    nasp.loc[U0[U0<=0.5].DateTime,:]= nasp_calm.loc[U0[U0<=0.5].DateTime,: ]  # calm regime
    nasp.loc[U0[(U0>0.5)&(U0 < 2)].DateTime,:]= slow.loc[U0[(U0>0.5)&(U0 < 2)].DateTime,:] # slow moving regime

    # moving regime conditions for theta_s < 60: 0.5 < Rv <2, 0.01 < Stk <6 ;  for theta_s >60:  0.003 < Stk < 0.2, 0.5 <Rv < 2 ## the lower thresholds for Stk were ignored since the efficiency is very close to 1 anyways
    nasp.loc[Rv[(Rv>=0.5)&(Rv <= 2)].DateTime,:]= nasp_move.loc[Rv[(Rv>=0.5)&(Rv <= 2)].DateTime,:]

    ###### Transmission Efficiency ######
        #(von der Weiden et al. (2009))
    #for 0.25 ≤ Rv ≤ 1 otherwise Iv = 0
        # losses in the vena contracta (I_v)
    I_v = 0.09 * (Stk*((U-U0)/U0 * np.cos(theta_sr)))**0.3
    I_v.loc[Rv[Rv >1].DateTime,:] = 0

    #for 0.02 ≤ Stk ≤ 4, and 0.25 ≤ Rv ≤ 4, where
        # losses by wall impact (I_w)
    wRv = Rv[(Rv < 0.25) | (Rv > 4)].DateTime
    iStk = Stk.where((Stk >= 0.02) & (Stk <= 4), float('nan'))

    # alpha only valid when calculated in deg not rad
    alpha =  np.deg2rad(12*((1-np.rad2deg(theta_sr)/90)-np.exp(-np.rad2deg(theta_sr))))
    y = np.sqrt(Rv)*np.sin(theta_sr-alpha)* np.sin((theta_sr-alpha)/2) #used for downward sampling (sampling probe faces upward),
    z = np.sqrt(Rv)*np.sin(theta_sr+alpha)* np.sin((theta_sr+alpha)/2) #for upward sampling (sampling probe faces downward) and
    I_w = iStk*y
    #for 0.25 ≤ Rv ≤ 4
    I_w.loc[wRv,:] = float('nan')

    #0.02 ≤ Stk ≤ 4, 0.25 ≤ Rv ≤ 4, and 0◦ ≤ θs ≤ 90◦:
        # combine to transmission efficiency (n_trm)
    tot = np.where(np.isnan(I_v), I_w, I_v+np.nan_to_num(I_w)) # add I_v and I_w with awareness to nan data
    n_trm = xr.DataArray(np.exp(-75*(tot)**2), dims = ['DateTime', 'D'], coords = {'DateTime' : index, 'D': D})
    iStk = Stk.where( (Stk <= 4), float('nan')) 
    n_trm = n_trm*iStk/iStk # *iStk/iStk makes everything nan that is not in the range of 0.02 ≤ Stk ≤ 4 
    n_trm.loc[wRv,:] = float('nan')

    ###### Transport Efficiency ######
        ## wind tunnel with constant diameter
    # sedimentation n_grav valid for V_ts*sin(theta)/TAS << 1
    a = 4*V_ts*Lw*np.cos(theta)
    b = TAS*d0*np.pi
    n_grav = np.exp(-(a/b))#xr.DataArray( np.exp(-(a/b)).T, dims = [ 'D'], coords = { 'D': D})

    # Turbulent inertial deposition n_turb_I
    tao_plus =  0.0395*Stk*Re**(1/8) # dimensionless particle relaxation time
    V_plus= 0.0006*tao_plus**2 # dimensionless turbulent velocity
    V_plus = V_plus.where(V_plus <= 0.0006*12.9**2, 0.1) # only for tao_plus <= 12.9, else V_plus = 0.1
    V_t = (V_plus*TAS)/(5.03*Re**(1/8)) #deposition velocity for turbulent inertial deposition

    n_turb_I =  xr.DataArray(np.exp(-(4*V_t*Lw)/(TAS*d0*np.pi)))

        ## in the contraction part with changing diameter
    # Inertial impaction in contraction n_cont 
    c = Stk*(1-A0/Ai)/(3.14*np.exp(-0.0185*np.rad2deg(theta_cont)))
    n_cont =xr.DataArray( 1- 1/(1+2*((c)**(-1.24)))) # it's -1.24 (see Muyshondt et al 1996), and not 1.24 like in spiegel et al. 2012

    # sedimentation in contraction n_cont_grav
    c = 4*V_ts*Li*np.cos(theta)
    d_ = TAS*di*np.pi
    n_cont_grav =np.exp(-(c/d_))

    # turbulent inertial deposition in contraction n_cont_turb_I
    n_cont_turb_I = xr.DataArray( np.exp(-(4*V_t*Li)/(TAS*di*np.pi)))

    # Combine all transport efficiencies
    n_tsp = n_grav*n_turb_I*n_cont*n_cont_grav*n_cont_turb_I
    # in nasp, n_trm and n_tsp make everything below 0.01 nan
    nasp = nasp.where(nasp > min_efficiency, float('nan'))
    n_trm = n_trm.where(n_trm > min_efficiency, float('nan'))
    n_tsp = n_tsp.where(n_tsp > min_efficiency, float('nan'))
    
    if return_sep == True:
        return nasp, n_trm, n_tsp
    else:
        ###### combine all efficiencies; where eff == nan, put 1 so we can still multiply 
        pnasp = xr.DataArray(np.nan_to_num(nasp, nan=1), dims=nasp.dims, coords=nasp.coords)
        pn_trm = xr.DataArray(np.nan_to_num(n_trm, nan=1), dims=n_trm.dims, coords=n_trm.coords)
        pn_tsp = xr.DataArray(np.nan_to_num(n_tsp, nan=1), dims=n_tsp.dims, coords=n_tsp.coords)
        total_efficiency = pnasp*pn_trm*pn_tsp
        return total_efficiency

    
    

############## function for returning corrected data #################
def correction(df, efficiency, min_efficiency = 0.01, pandas = True):
    ''' This function applies the calculated efficiencies on the size distribution data.
        input:
        df: GFAS size distribution data (only). has to be the same size bins as used in loss_calc_GFAS101
        efficiency: the overall sampling efficiency calculated in loss_calc_GFAS101 with return_sep = False (total_efficiency)
        min_efficiency: cutoff for when to apply efficiency value. if too low sample loss corrections might skyrocket
        pandas: if True returns pandas df, else xarray dataset
        
        returns: corrected size distribution

    '''    
    df = df.resample('1min').mean() # make sure that the data is in 1 min resolution
    # df: pandas dataframe with only size data!
    idx = np.intersect1d(df.index, efficiency.DateTime)
    df = df.loc[idx]
    efficiency = efficiency.loc[idx]
    D = df.columns.astype(float).values
    if (D > 0.1).any():
        D = D/10**6 # convert to m if in µm

    index = df.index
    corrected = xr.DataArray(df, dims=['DateTime', 'D'], coords={'DateTime': index, 'D':D})/(efficiency.where(efficiency > min_efficiency, 1))
    if pandas == True:
        corrected = pd.DataFrame(corrected, columns = D, index = index)
        return corrected
    else: return corrected

################ function for calculating N, LWC, and ED ################
def correct_parameter (corrected_data, min_val_series,  Dmin =0.4499999880791, min_valED =1 , log = False):
    ''' corrected_data: pandas dataframe with corrected data (only number concentration in bins)
        min_val_series: df from which the minimum value (min_valED) for the mask is supposed to be taken 
        Dmin: lower bin limit in µm
        min_valED: minimum value for the number concentration to calculate ED
        log: if True: size distribution was in log space, if False: SD was in natural space
        returns: pandas dataframe with dN/dlogDp, N, ED, and LWC
    '''

    df = corrected_data.resample('1min').mean() # make sure that the data is in 1 min resolution
    min_val_series = min_val_series.resample('1min').mean() # make sure that the data is in 1 min resolution
    
    d = df.columns.astype(float).values
    if (d > 0.1).any():
        d = d/10**6 # convert to m if in µm
    D = d-(d-shift(d, 1, cval = np.nan))/2 # mid diameter
    D[0] = d[0]-Dmin*10**(-6)/2
    
    if log == False:
        dlogDp= D-shift( D,1,cval=np.nan)
        dlogDp[0]= D[0]- (Dmin*10**(-6))
        dlogDp = dlogDp*10**6
    else:

        dlogDp=np.log10(D)-shift(np.log10(D),1,cval=np.nan)
        dlogDp[0]=np.log10(D)[0]-np.log10(Dmin*10**(-6)) 

    # make mask where minimum value in the min_val_series is larger than min_valED, and get the index of that mask to use for calculating ED; 
    # this is to make sure that we only calculate ED for size bins where there are enough particles to make a meaningful calculation
     
    mask = min_val_series[min_val_series.values>min_valED].index # minimum criteria to calculate ED 

    N =(df*dlogDp).sum(axis = 1)
    ED = (df.loc[mask]*(D*10**(6))**3).sum(axis = 1)/(df.loc[mask]*(D*10**(6))**2).sum(axis=1)
    LWC = (((df*dlogDp)*(D*10**6)**3).sum(axis = 1)*np.pi/6)/10**6

    new = pd.concat([N,ED,LWC], axis = 1, keys = ['N','ED','LWC'])
    NEW = pd.concat([df,new], axis = 1)
    return NEW
#%%
##################### flagging functions ############################

def efficiency_flags(n_asp, n_trm, n_tsp):
    ''' n_asp: aspiration efficiency
        n_trm: transmission efficiency
        n_tsp: transport efficiency

        For interpreation of flagging:
    0 = all three efficiencies (nasp, n_trm, ntsp) were present
    1 = nasp was missing (NaN)
    2 = n_trm was missing (NaN)
    4 = n_tsp was missing (NaN)
        Combine flags to represent cases where more than one value was missing, e.g.,
    3 = nasp and n_trm were missing (1+2), etc.
    5 = nasp and n_tsp were missing (1+4)
    6 = n_trm and n_tsp were missing (2+4)
    7 = nasp, n_trm, and n_tsp were all missing (1+2+4)

    returns: flagging matrix with dim = time x D
    '''
    flag = np.zeros_like(n_asp)
    # Create the flagging matrix
    flag += np.isnan(n_asp).astype(int) * 1  # Flag as 1 where nasp is NaN
    flag += np.isnan(n_trm).astype(int) * 2  # Flag as 2 where n_trm is NaN
    flag += np.isnan(n_tsp).astype(int) * 4  # Flag as 4 where n_tsp is NaN
    return flag

def classify_quality(flag):
    # Apply a reduction over the 'D' dimension to classify each timestep
    # Case 0: Fully Assured -> If all flags in the 'D' dimension are 0
    fully_assured = (flag == 0).all(dim='D')
    
    # Case 1: Partially Assured -> If there is any flag that is 1, 2, or 4 (only one efficiency missing)
    partially_assured = ((flag == 1) | (flag == 2) | (flag == 4)).any(dim='D')
    
    # Case 2: Poorly Assured -> If there is any flag that is 3, 5, or 6 (two efficiencies missing)
    poorly_assured = ((flag == 3) | (flag == 5) | (flag == 6)).any(dim='D')
    
    # Case 3: Not Assured -> If all flags in the 'D' dimension are 7 (all efficiencies missing)
    not_assured = (flag == 7).all(dim='D')
    
    # Initialize an array for quality levels (default to fully assured)
    quality_assurance = xr.DataArray(np.zeros_like(fully_assured, dtype=int), 
                                     coords={'DateTime': flag['DateTime']},
                                     dims=['DateTime'])
    
    # Assign values for each classification
    quality_assurance = xr.where(not_assured, 3, quality_assurance)  # First prioritize "Not Assured"
    quality_assurance = xr.where(poorly_assured & ~not_assured, 2, quality_assurance)  # Poorly Assured (exclude Not Assured)
    quality_assurance = xr.where(partially_assured & ~poorly_assured & ~not_assured, 1, quality_assurance)  # Partially Assured (exclude above)
    
    # The remaining values are 0, which corresponds to "Fully Assured"
    
    return quality_assurance

import matplotlib.colors as colors
def plot_flag_parameters (flag, cbar):
    ''' defining cmap, norm, midpoints and tick labels for flagging plot'''
    # Plot the DataArray with a discrete colormap
    cmap = plt.get_cmap("Set3", int(flag.max())+1)  # Discrete colormap with 4 colors
    norm = colors.BoundaryNorm(boundaries= np.linspace(-0.5, int(flag.max()) + 0.5, int(flag.max()) + 2), ncolors= int(flag.max())+1)

    # Set tick positions to the middle of each color (midpoints of the boundaries)
    midpoints = np.linspace(0, int(flag.max()),int(flag.max())+1)
    flag_names = ['All available', 
                  'Transmission \nTransport', 
                  'Aspiration \nTransport', 
                  'Transport', 
                  'Aspiration \nTransmission', 
                  'Transmission', 
                  'Aspiration', 
                  'All missing']

    # Set the tick locations and custom labels
    cbar.set_ticks(midpoints)
    cbar.set_ticklabels(flag_names[0:int(flag.max())+1])
    