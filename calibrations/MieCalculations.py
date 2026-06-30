'''
This is the sub-routine for Mie calculations icluding:
- calculating the scattering cross section over a given opening angle and for a range of diameters


author: Lea Haberstock, Stockholm University, Department of Environmental Science, Atmospheric Unit 
developed toether with: Darrel Baumgardner and Paul Zieger

contact: lea.haberstock@aces.su.se 
used in publication: Haberstock et al. 2026 (submitted to AMT) 


Last modified June 26 2026
'''

import PyMieScatt  as pms
import numpy as np
import pandas as pd


#%% Mie calculations

## P pol and Spol scattering cross section vs Dp for different m
def sct_SLSR(m, diam, wvl, ang):
    ''' 
    calculate the scattering cross section for an array of diameters. It calculates the scattering cross section for 30 steps between min and maximum opening angle integrates over the opening angle
    
    m: refractive index as a single value (complex or natural)
    diam: array of diameters
    wvl: wavelength of the incident light
    ang: array of scattering angle with a minimum of two values, the min and the max opening angle
    '''
    sc_L = np.full((len(diam)), np.nan)  # output vector
    sc_R = np.full((len(diam)), np.nan)  # output vector
    sc_TOT = np.full((len(diam)), np.nan)  # output vector
    angle_res = (ang[-1]-ang[0])/30
    for k1, d in enumerate(diam):
        x = d * np.pi / wvl  # size parameter
        theta, SL, SR, SU = pms.ScatteringFunction(m=m, wavelength=wvl, diameter=d, nMedium=1.0,  minAngle=ang[0], maxAngle=ang[-1], angularResolution= angle_res)
        sc_l = x ** (-2) * np.trapz((SL) * np.sin(theta), theta)
        sc_r = x ** (-2) * np.trapz((SR) * np.sin(theta), theta)
        sc_tot =x ** (-2) * np.trapz((SL+SR) * np.sin(theta), theta)
        sc_L[k1] = sc_l * 0.25 * np.pi * (d * 1e-9) ** 2
        sc_R[k1] = sc_r * 0.25 * np.pi * (d * 1e-9) ** 2
        sc_TOT[k1] = sc_tot * 0.25 * np.pi * (d * 1e-9) **2
        sc_L = sc_L 
        sc_R = sc_R 
        sc_TOT = sc_TOT 
    return sc_L, sc_R, sc_TOT