# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 11:31:27 2026

@author: siirias
"""
import numpy as np
from scipy.optimize import brentq
from dataclasses import dataclass
import matplotlib.pyplot as plt
import xarray as xr
import argopy
import gsw
from argopy import DataFetcher, ArgoFloat
from argopy import ArgoIndex  #  This is the class to work with Argo index
from argopy.plot import scatter_map, scatter_plot  # Functions to easily make maps and plots

# Shut down some warning messages for reading
import warnings
warnings.filterwarnings("ignore")
xr.set_options(display_expand_attrs = False)
LOAD_DATA = True  #WHen testing aroudn several times, no need to redownload the data every time.


@dataclass
class SBE41cpCoefficients:
    """Calibration coefficients for SBE 41cp conductivity sensor."""
    ## f = INST_FREQ * sqrt(1.0 + WBOTC * t) / 1000.0
    ## Conductivity = (g + h*f^2 + i*f^3 + j*f^4) / (1 + CTcor*Temp + CPcor*Pres)
    g: float = -1.028007e+000
    h: float =  1.495265e-001
    i: float = -3.432114e-004
    j: float =  4.724081e-005
    CPcor: float = -9.5700e-008   # epsilon (pressure correction)
    CTcor: float =  3.2500e-006   # delta   (temperature correction)
    WBOTC: float = -4.2500e-007   # wire bridge offset temp coefficient


def inst_freq_to_conductivity(inst_freq, temp_C, pres_dbar,
                               coeffs: SBE41cpCoefficients) -> np.ndarray:
    """
    Convert raw instrument frequency to conductivity (Siemens/meter).
    
    Parameters
    ----------
    inst_freq : array-like, Hz
    temp_C    : array-like, degrees C (ITS-90)
    pres_dbar : array-like, decibars
    coeffs    : SBE41cpCoefficients
    
    Returns
    -------
    conductivity in S/m
    """
    inst_freq = np.asarray(inst_freq, dtype=float)
    temp_C    = np.asarray(temp_C,    dtype=float)
    pres_dbar = np.asarray(pres_dbar, dtype=float)

    f = inst_freq * np.sqrt(1.0 + coeffs.WBOTC * temp_C) / 1000.0

    cond = (coeffs.g + coeffs.h * f**2 + coeffs.i * f**3 + coeffs.j * f**4) / \
           (1.0 + coeffs.CTcor * temp_C + coeffs.CPcor * pres_dbar)
    return cond



def conductivity_to_inst_freq(target_cond, temp_C, pres_dbar, coeffs,
                             f_min=1000.0, f_max=10000.0):
    """
    Invert conductivity to instrument frequency using root finding.

    This solves:
        inst_freq_to_conductivity(freq, T, P) - target_cond = 0

    Parameters
    ----------
    target_cond : array-like Conductivity (S/m)
    temp_C : array-like      Temperature (°C, ITS-90)
    pres_dbar : array-like   Pressure (dbar)
    coeffs : SBE41cpCoefficients     Calibration coefficients
    f_min : float, optional  Minimum search frequency (Hz)
    f_max : float, optional  Maximum search frequency (Hz)

    Returns
    -------
    inst_freq : ndarray
        Instrument frequency (Hz), same shape as inputs
    """

    # Convert inputs to arrays and broadcast to same shape
    target_cond = np.asarray(target_cond, dtype=float)
    temp_C = np.asarray(temp_C, dtype=float)
    pres_dbar = np.asarray(pres_dbar, dtype=float)

    target_cond, temp_C, pres_dbar = np.broadcast_arrays(
        target_cond, temp_C, pres_dbar
    )

    # Output array
    inst_freq = np.empty_like(target_cond, dtype=float)

    # Flatten for iteration (brentq is scalar)
    flat_cond = target_cond.ravel()
    flat_temp = temp_C.ravel()
    flat_pres = pres_dbar.ravel()
    flat_out = inst_freq.ravel()

    for i in range(flat_cond.size):
        c = flat_cond[i]
        t = flat_temp[i]
        p = flat_pres[i]

        # Define residual function for this point
        def residual(f):
            cond_val = inst_freq_to_conductivity(f, t, p, coeffs)
            return float(cond_val) - c        # Solve for frequency
        flat_out[i] = brentq(residual, f_min, f_max)

    return inst_freq

def conductivity_to_psal(cond_Sm, temp_C, pres_dbar):
    cond_mScm = np.asarray(cond_Sm, dtype=float) * 10.0
    SP = gsw.SP_from_C(cond_mScm, temp_C, pres_dbar)
    return np.nan_to_num(SP, nan=0.0)

def psal_to_conductivity(psal, temp_C, pres_dbar,
                         cmin_Sm=0.001, cmax_Sm=30.0):
    psal = np.asarray(psal, dtype=float)
    temp_C = np.asarray(temp_C, dtype=float)
    pres_dbar = np.asarray(pres_dbar, dtype=float)

    psal, temp_C, pres_dbar = np.broadcast_arrays(psal, temp_C, pres_dbar)
    out = np.empty_like(psal, dtype=float)

    flat_psal = psal.ravel()
    flat_t = temp_C.ravel()
    flat_p = pres_dbar.ravel()
    flat_out = out.ravel()

    for i in range(flat_psal.size):
        sp = flat_psal[i]
        t = flat_t[i]
        p = flat_p[i]

        if np.isnan(sp) or np.isnan(t) or np.isnan(p):
            flat_out[i] = np.nan
            continue
        if np.abs(sp) <= 1e-6:
            flat_out[i] = 0.0
            continue
        def residual(c_Sm):
            sp_est = gsw.SP_from_C(c_Sm * 10.0, t, p)
            return float(sp_est) - sp
        
        flat_out[i] = brentq(residual, cmin_Sm, cmax_Sm)

    return out

#argopy.set_options(ds='bgc')
argopy.set_options(ds='phy')

float_no = 6902013
#float_no = 6901773

if LOAD_DATA:
    idx = ArgoIndex(index_file='core').load() 
    #idx = ArgoIndex(index_file='bgc-b').load() 
    data = DataFetcher(mode = 'expert', params = 'all').float(float_no).load()


# Plot the profile cloud, to see what we are workign with
plt.figure()
for cycle_num, profile_ds in data.data.groupby("CYCLE_NUMBER"):
    S = profile_ds.PSAL
    P = profile_ds.PRES
    plt.plot(S,-1.0*P)
    plt.grid(True)

coeffs = SBE41cpCoefficients()
coeffs2 = SBE41cpCoefficients(
    g=-1.016771e+000,
    h=1.479133e-001,
    i=-3.481665e-004,
    j=4.656400e-005,
    CPcor=-9.5700e-008,
    CTcor=3.2500e-006,
    WBOTC=-4.2500e-007
)
coeffs3 = SBE41cpCoefficients(
    g=-1.017400e+000,
    h=1.481721e-001,
    i=-4.263132e-004,
    j=5.244764e-005,
    CPcor=-9.5700e-008,
    CTcor=3.2500e-006,
    WBOTC=-4.2500e-007
)

# Plot the difference with two sets of calibration parameters as a cloud
plt.figure()
for cycle_num, test_d in data.data.groupby("CYCLE_NUMBER"):
    S = test_d.PSAL
    T = test_d.TEMP
    P = test_d.PRES
    Cond = psal_to_conductivity(S,T,P)
    f = conductivity_to_inst_freq(Cond, T,P, coeffs)
    Cond_new = inst_freq_to_conductivity(f, T, P, coeffs2)
    S_new = conductivity_to_psal(Cond_new,T,P)
    plt.plot(S_new-S,-1.0*P)
    plt.grid(True)



# Test the round trip of fucntions, to see how much error comes,
# When the data is rotated:
# frequency -> conductivity -> Salinity -> conductivity frequency
Temp = [22.0, 1.0, 4.4999, 15.0, 18.5, 23.994, 29.0, 32.50]
Pres = [0.0,  0.0, 0.0, 0.0,  0.0,  0.0,    0.0,  0.0]
Freq = [2627.11, 5185.88, 5380.60,5959.86,6150.29,6445.61,6710.37,6892.70]
Conduct = [0.0, 2.97956,  3.28699, 4.26986, 4.61540, 5.17339, 5.69650, 6.06937]
cond = inst_freq_to_conductivity(Freq, Temp, Pres, coeffs)
Salinities = []
for i,T,P in zip(cond,Temp,Pres):
    sal = conductivity_to_psal(i,T,P)
    Salinities.append((sal))
    print(f"Conductivity: {i:.5f} S/m = Salinity: {sal:.5f}")  
new_cond = psal_to_conductivity(Salinities,Temp, Pres)
new_freq = conductivity_to_inst_freq(new_cond, Temp, Pres, coeffs)
for n_f,o_f,n_c,o_c in zip(new_freq,Freq,new_cond,cond):
    print(f"New Frequency: {n_f:.2f} old {o_f:.2f} diff:{n_f-o_f:.4f}")
