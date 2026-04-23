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
from argopy import DataFetcher, ArgoFloat
from argopy import ArgoIndex  #  This is the class to work with Argo index
from argopy.plot import scatter_map, scatter_plot  # Functions to easily make maps and plots

# Shut down some warning messages for reading
import warnings
warnings.filterwarnings("ignore")
xr.set_options(display_expand_attrs = False)


@dataclass
class SBE41cpCoefficients:
    """Calibration coefficients for SBE 41cp conductivity sensor."""
    ## f = INST_FREQ * sqrt(1.0 + WBOTC * t) / 1000.0
    ## Conductivity = (g + h*f + i*f² + j*f³) / (1 + δ*t + ε*p)
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


#argopy.set_options(ds='bgc')
argopy.set_options(ds='phy')

LOAD_DATA = True  #WHen testing aroudn several times, no need to redownload the data every time.
float_no = 6902013
#float_no = 6901773

if LOAD_DATA:
    idx = ArgoIndex(index_file='core').load() 
    #idx = ArgoIndex(index_file='bgc-b').load() 
    data = DataFetcher(mode = 'expert', params = 'all').float(float_no).load()


#scatter_map(data.index, set_global=False)
#data.plot('trajectory')
#scatter_plot(data.data, 'TEMP')
plt.figure()
for cycle_num, profile_ds in data.data.groupby("CYCLE_NUMBER"):
    S = profile_ds.PSAL_ADJUSTED
    P = profile_ds.PRES_ADJUSTED
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
plt.figure()
for cycle_num, test_d in data.data.groupby("CYCLE_NUMBER"):
    S = test_d.PSAL_ADJUSTED
    T = test_d.TEMP_ADJUSTED
    P = test_d.PRES_ADJUSTED
    freq = conductivity_to_inst_freq(S, T,P, coeffs)
    S_new = inst_freq_to_conductivity(freq, T, P, coeffs3)
    plt.plot(S_new-S,-1.0*P)
    plt.grid(True)



#Test the adjustment.

# cond = inst_freq_to_conductivity([2627.11, 5185.88, 5380.60], [22.0, 1.0,5.0], 0.0, coeffs)
# for i in cond:
#     print(f"Conductivity: {i:.5f} S/m")  
# freq = conductivity_to_inst_freq([0.0,2.97955,3.28698], [22.0, 1.0,5.0], 0.0, coeffs)
# for i in freq:
#     print(f"Frequency: {i:.5f}")
