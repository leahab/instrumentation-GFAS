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
        theta, SL, SR, SU = pms.ScatteringFunction(m=m, wavelength=wvl, diameter=d, nMedium=1.0,  
                                                   minAngle=ang[0], maxAngle=ang[-1], angularResolution= angle_res)
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

def mie_calculation(diameters, m=1.33):
     '''
     Process the Mie Calculations for a defined list of diameters and a range of diameters

     parameters:
     diameters : 1D-array    List of diameters
     m : float or complex    Refractive index
     output: scattering cross sections for 
        sigma_fw: forward scattering
        sigma_S: perpendicularly polarized backward scattering
        sigma_P: parallel polarized backward scattering
        sigma_x: measurable parallel polarized backward scattering
        sigma_y: measurable perpendicularly polarized backward scattering
     '''

     n=len(diameters)
     ang_bw = np.linspace(168, 176.5, 30) # backscattering opening angles
     ang_fw = np.linspace(3.8, 20, 30)  # forward scattering opening angles

     sigma_fw=np.zeros(n)
     sigma_S=np.zeros(n)
     sigma_P=np.zeros(n)
     sigma_J=np.zeros(n)
     for i,di in enumerate(diameters):
         x = di * np.pi / 660  # size parameter, unities [nm]/[nm]
         theta_fw, SL1, SR1, _ = pms.ScatteringFunction(m=m, wavelength=660, diameter=di, nMedium=1.0,  
                                                        minAngle=ang_fw[0], maxAngle=ang_fw[-1], angularResolution= (ang_fw[-1]-ang_fw[0])/100 )
         theta_bw, SL, SR, _ = pms.ScatteringFunction(m=m, wavelength=660, diameter=di, nMedium=1.0,  
                                                      minAngle=ang_bw[0], maxAngle=ang_bw[-1], angularResolution= (ang_bw[-1]-ang_bw[0])/100 )
         sc_tot =x ** (-2) * np.trapz((SL1+SR1) * np.sin(theta_fw), theta_fw)
         sc_L =x ** (-2) * np.trapz((SL) * np.sin(theta_bw), theta_bw)
         sc_R =x ** (-2) * np.trapz((SR) * np.sin(theta_bw), theta_bw)
         sigma_fw[i]=1e4*sc_tot*0.25 * np.pi * (di * 1e-9) ** 2 #di converted to meters, sigma converted to cm^2
         sigma_P[i]=1e4*sc_R*0.25 * np.pi * (di * 1e-9) ** 2 #di converted to meters, sigma converted to cm^2
         sigma_S[i]=1e4*sc_L*0.25 * np.pi * (di * 1e-9) ** 2 #di converted to meters, sigma converted to cm^2

         S1,S2=np.zeros(30, dtype=complex), np.zeros(30, dtype=complex)
         theta= ang_bw*np.pi/180 # into rad
         for j in range(30):
             S1[j],S2[j]=pms.MieS1S2(m, x, np.cos(theta[j]))
         sigma_J[i]= x**(-2)/4 *1e4 *np.pi * (di * 1e-9)**2 * np.trapz(np.real(S1 * np.conjugate(S2))* np.sin(theta), theta)

     sigma_x=.75*(sigma_S+sigma_P)-.5*sigma_J
     sigma_y=.25*(sigma_S+sigma_P)+.5*sigma_J

     return sigma_fw, sigma_S, sigma_P, sigma_x, sigma_y