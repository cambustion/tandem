"""
==============================================================
Tandem Aerosol Data Fitting and Multi-Lognormal Analysis Tool
==============================================================
Description:
------------
This script performs multi-lognormal fitting of particle size
distribution data obtained from tandem experiments (e.g., CPMA-SMPS, 
AAC-SMPS, CPMA-AAC). It processes data natively using NumPy, allows 
manual selection of peaks via an interactive UI, performs curve fitting, 
and dynamically applies mobility/aerodynamic diameter calculations 
based on the specified classifier metadata.
"""
#%%
import os
import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import curve_fit, fsolve
from datetime import datetime
from typing import Tuple, List, Dict, Any, Callable

import tandem_data_utils_2

# =============================================================
# Constants
# =============================================================
E_CHARGE = 1.60218e-19  # Electron charge [C]


# =============================================================
# Physical Helper Functions
# =============================================================
def Cc(d_m: np.ndarray, P_atm: np.ndarray, T_K: np.ndarray) -> np.ndarray:
    """Compute the Cunningham slip correction factor for arrays (diameter in meters)."""
    alpha, beta, gamma = 1.165 * 2, 0.483 * 2, 0.997 / 2
    lambda_air_std = 67.30e-9

    lambda_air = lambda_air_std * (T_K / 296.15)**2 * (1 / P_atm) * ((110.4 + 296.15) / (T_K + 110.4))
    return 1 + (lambda_air / d_m) * (alpha + beta * np.exp(-gamma * d_m / lambda_air))


def compute_mobility_diameter(d_m_initial: float, P_atm: float, T_K: float, charge_number: int) -> Tuple[float, float]:
    """Compute electrical mobility and the corresponding diameter for a given charge number."""
    mu = 1.81809e-5 * (T_K / 293.15)**1.5 * (293.15 + 110.4) / (T_K + 110.4)

    z = E_CHARGE * Cc(d_m_initial, P_atm, T_K) / \
        (3 * np.pi * mu * d_m_initial)

    def equation(d_m_val):
        return d_m_val - (charge_number * E_CHARGE *
                          Cc(d_m_val, P_atm, T_K)) / (3 * np.pi * mu * z)

    d_m_solution = fsolve(equation, charge_number * d_m_initial)
    return z, d_m_solution[0]

def newton(fun, x0, n=3):
    """
    Solve for a zero using Newtons method.
    """
    n = 3 if n is None else n  # set default
    for _ in range(n):
        f = fun(x0)
        
        # FIX 1: Use a relative step size so large numbers don't lose precision
        h = x0 * 1e-5 
        fprime = (fun(x0 + h) - f) / h
        
        # FIX 2 (Vectorized): Safely handle division by zero across arrays
        with np.errstate(divide='ignore', invalid='ignore'):
            step = f / fprime
            
            # If fprime is 0, step becomes inf or nan. We force those specific 
            # elements to 0.0 so their x0 guess simply stays where it is.
            step = np.where(np.isfinite(step), step, 0.0)
            
        x0 = x0 - step
        
        # BOUNDARY FIX: Prevent negative diameters
        x0 = np.maximum(x0, 1e-9) 
        
    return x0

def fzero(fun, x0, n=3):
    """
    Wrapper to decide which zero solver to use.
    """
    if n == 0:
        return x0  # do not iterate
    elif n == np.inf:
        return fsolve(fun, x0)  # then use imported solver to default tolerances
    else:
        return newton(fun, x0, n)  # else do fast solve, with limited computations


def mp_da2dm(mp: np.ndarray, da: np.ndarray, T: np.ndarray, P: np.ndarray, n_iter: int = 3) -> np.ndarray:
    """
    Compute mobility diameter (m) from particle mass (kg) and aerodynamic diameter (m).
    Uses the custom fzero/newton solver to process arrays efficiently.
    """
    # Ensure inputs are arrays
    mp = np.asarray(mp, dtype=float)
    da = np.asarray(da, dtype=float)
    T = np.asarray(T, dtype=float)
    P = np.asarray(P, dtype=float)

    # 1. Calculate the base guess (dm without slip correction ratio applied)
    dm_base = (da**-2) * mp * (6.0 / (np.pi * 1e3))
    
    if n_iter <= 0:
        return dm_base
        
    # 2. Precompute the slip correction for the aerodynamic diameter (remains constant)
    cc_da = Cc(da, P, T)
    
    # 3. Define the function to solve
    # We want to find the root where: dm_guess = dm_base * (Cc(dm_guess) / Cc(da))
    def fun(dm_guess):
        cc_dm = Cc(dm_guess, P, T)
        return 1e9 * ((dm_base * (cc_dm / cc_da)) - dm_guess)

    # 4. Solve for mobility diameter using the fzero wrapper
    dm_final = fzero(fun, dm_base, n_iter)
        
    return dm_final

def dm_da2mp(dm: np.ndarray, da: np.ndarray, T: np.ndarray, P: np.ndarray) -> np.ndarray:
    """
    Compute mass (kg) from mobility and aerodynamic diameters (m).
    """
    dm = np.asarray(dm, dtype=float)
    da = np.asarray(da, dtype=float)
    T = np.asarray(T, dtype=float)
    P = np.asarray(P, dtype=float)

    # Precompute the slip correction
    cc_da = Cc(da, P, T)
    cc_dm = Cc(dm, P, T)

    mp = da**2 * dm * (np.pi * 1e3 / 6 * cc_da / cc_dm)

    return mp


# =============================================================
# Setpoint Calculation Functions 
# =============================================================
def generate_cpma_setpoints(cpma_spd, cpma_v) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate CPMA setpoint masses (m_star) across all data rows 
    using vectorized NumPy math.
    """
    
    # 1. Perform math on the entire arrays instantly (Vectorization)
    omega = np.asarray(cpma_spd)
    V = np.asarray(cpma_v)
    
    # CPMA properties
    r1 = 0.06
    r2 = 0.061
    r_hat = r1 / r2
    rc = (r1 + r2)/2
    omega_hat = 32 / 33
    
    # Calculate terms for the whole array
    omega1 = omega / (((r_hat**2 - omega_hat) / (r_hat**2 - 1)) + 
                      (r1**2 * (omega_hat - 1) / ((r_hat**2 - 1) * rc**2)))
                      
    alpha = omega1 * (r_hat**2 - omega_hat) / (r_hat**2 - 1)
    beta = omega1 * r1**2 * (omega_hat - 1) / (r_hat**2 - 1)
    
    # Calculate target mass (m_star)
    m_star = V / ((np.log(1 / r_hat) / E_CHARGE) * (alpha * rc + beta / rc)**2)
    
    # 2. Package results into the dictionary list format expected by the main script
    sp_list = []
    n = len(omega)
    
    for i in range(n):       
        sp = {
            'Mp': m_star[i],  # Map 'Mp' key for the main processing loop
            'omega': omega[i],
            'V': V[i],
            'omega1': omega1[i],
            'alpha': alpha[i],
            'beta': beta[i]
        }
        sp_list.append(sp)
        
    return np.asarray(sp_list, dtype=object), m_star


def tau_to_da(tau_star_arr, mu_arr, T_arr, P_arr):
    """Helper to invert tau_star back to aerodynamic diameter for an array of values."""
    da_sols = np.zeros_like(tau_star_arr)
    
    for i in range(len(tau_star_arr)):
        tau, mu, T, P = tau_star_arr[i], mu_arr[i], T_arr[i], P_arr[i]
        def equation(da):
            cc = Cc(np.array([da]), np.array([P]), np.array([T]))[0]
            return (cc * 1e3 * da**2 / (18 * mu)) - tau
            
        guess = np.sqrt(tau * 18 * mu / 1e3)
        da_sols[i] = fsolve(equation, guess)[0]
        
    return da_sols

def generate_aac_setpoints(aac_da=None, aac_spd=None, aac_qsh=None, aac_T=None, aac_P=None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate AAC setpoint diameters (d_star) and speeds (omega) across all data rows.
    """
    da = np.asarray(aac_da if aac_da is not None else [], dtype=float)
    spd = np.asarray(aac_spd if aac_spd is not None else [], dtype=float)

    da = np.atleast_1d(da) if da.size > 0 else da
    spd = np.atleast_1d(spd) if spd.size > 0 else spd

    n = max(len(da), len(spd))
    
    if n == 0:
        return np.array([]), np.array([])
        
    # 1. Standardize input arrays to ensure they are the correct length
    if len(da) < n: da = np.full(n, np.nan)
    if len(spd) < n: spd = np.full(n, np.nan)
    
    qsh = np.asarray(aac_qsh, dtype=float) if (aac_qsh is not None and np.size(aac_qsh) > 0) else np.full(n, 3.0)
    T = np.asarray(aac_T, dtype=float) if (aac_T is not None and np.size(aac_T) > 0) else np.full(n, 293.15)
    P = np.asarray(aac_P, dtype=float) if (aac_P is not None and np.size(aac_P) > 0) else np.full(n, 1.0)

    # Broadcast single values out to the full length of n if needed
    if qsh.size == 1 and n > 1: qsh = np.full(n, qsh.item())
    if T.size == 1 and n > 1: T = np.full(n, T.item())
    if P.size == 1 and n > 1: P = np.full(n, P.item())

    Qexh = qsh
    Qtot = qsh + Qexh

    # AAC properties
    r1 = 0.056  
    r2= 0.060  
    L = 0.206   
    
    # 2. Calculate constants and properties for the whole array
    mu_ref = 1.83245e-5   
    S = 110.4
    mu = mu_ref * (T / 296.15)**1.5 * (296.15 + S) / (T + S)
    area_term = np.pi * (r2 + r1)**2 * L

    # Final outputs
    tau_star_final = np.full(n, np.nan)
    omega_final = np.full(n, np.nan)
    da_final = np.full(n, np.nan)

    # 3. Vectorized Math Condition A: Diameter is provided
    mask_da = ~np.isnan(da)
    if np.any(mask_da):
        da_m = da[mask_da] * 1e-9
        cc_vals = Cc(da_m, P[mask_da], T[mask_da])
        tau_da = cc_vals * 1e3 * da_m**2 / (18 * mu[mask_da])
        
        tau_star_final[mask_da] = tau_da
        omega_final[mask_da] = np.sqrt(Qtot[mask_da] / (area_term * tau_da))
        da_final[mask_da] = da[mask_da]

    # 4. Vectorized Math Condition B: Speed is provided, but Diameter is missing
    mask_spd = (~np.isnan(spd)) & (~mask_da)
    if np.any(mask_spd):
        tau_spd = Qtot[mask_spd] / (area_term * spd[mask_spd]**2)
        
        # Invert tau to diameter using fsolve only for rows where diameter was missing
        da_m_spd = tau_to_da(tau_spd, mu[mask_spd], T[mask_spd], P[mask_spd])
        
        tau_star_final[mask_spd] = tau_spd
        omega_final[mask_spd] = spd[mask_spd]
        da_final[mask_spd] = da_m_spd * 1e9

    if np.all(np.isnan(da_final)):
        raise ValueError("For AAC you must specify valid arrays for either 'aac_da' (nm) or 'aac_spd' (rad/s).")

    # 5. Package results into the dictionary list format expected by the main script
    sp_list = []
    for i in range(n):
        sp = {
            'aac_da': da_final[i],  # Map expected key for main script
            'tau_star': tau_star_final[i],
            'omega': omega_final[i],
            'Qsh': qsh[i] if np.ndim(qsh) > 0 else qsh,
            'Qexh': Qexh[i] if np.ndim(Qexh) > 0 else Qexh,
            'T': T[i],
            'P': P[i]
        }
        sp_list.append(sp)

    return np.asarray(sp_list, dtype=object), da_final


# =============================================================
# Curve Fitting Functions
# =============================================================
def multi_lognormal_distribution(x: np.ndarray, *params) -> np.ndarray:
    """
    Sum of lognormal modes expressed as dN/dlog10(x).

    Parameterization:
        mu      = ln(CMD)
        sigma   = ln(GSD)
        N_total = integrated concentration of the mode
    """
    x = np.asarray(x, dtype=float)
    result = np.zeros_like(x)
    ln10 = np.log(10.0)

    for i in range(0, len(params), 3):
        mu, sigma, N_total = params[i:i + 3]
        result += (N_total * ln10 / (sigma* np.sqrt(2.0 * np.pi)) * np.exp(-0.5* ((np.log(x) - mu)/ sigma)**2))
    return result

def generate_initial_guesses(x: np.ndarray, y: np.ndarray, selector: 'InteractivePeakSelector') -> np.ndarray:
    """Generate initial parameter guesses based on user-selected peaks/valleys."""
    initial_guess = []
    
    if not selector.idx_max:
        area = max(0.0, np.trapezoid(y, x=np.log10(x)))

        logx = np.log(x)
        weights = np.maximum(y, 0.0)
        if np.sum(weights) > 0:
            mu0 = np.sum(weights * logx) / np.sum(weights)
            sigma0 = np.sqrt(np.sum(weights * (logx - mu0)**2) / np.sum(weights))

        else:
            mu0 = np.nanmean(logx)
            sigma0 = np.nanstd(logx)
        sigma0 = sigma0 if np.isfinite(sigma0) and sigma0 > 1e-3 else 0.3
        return np.array([mu0, sigma0, area])

    indices = sorted(list(set([0] + selector.idx_max + selector.idx_min + [len(y) - 1])))

    # Sub-segment guesses
    for j in range(0, len(indices) - 2, 2):
        idx_start, idx_end = indices[j], indices[j+2]
        x_seg = x[idx_start : idx_end + 1]
        y_seg = y[idx_start : idx_end + 1]

        if len(x_seg) < 3: 
            continue

        scale = max(1e-5, np.trapezoid(y_seg, x=np.log10(x_seg)))
        logx_seg = np.log(x_seg)
        weights_seg = np.maximum(y_seg, 0.0)

        if np.sum(weights_seg) > 0:
            mu_guess = np.sum(weights_seg * logx_seg) / np.sum(weights_seg)
            sigma_guess = np.sqrt(np.sum(weights_seg * (logx_seg - mu_guess)**2) / np.sum(weights_seg))

        else:
            mu_guess = np.mean(logx_seg)
            sigma_guess = 0.2

        sigma_guess = (sigma_guess if np.isfinite(sigma_guess) and sigma_guess > 1e-3 else 0.2)

        guess = (mu_guess, sigma_guess, scale)
        try:
            # Added physical bounds to sub-segment fits to prevent negative guesses
            popt_sub, _ = curve_fit(
                multi_lognormal_distribution, 
                x_seg, y_seg, 
                p0=guess, 
                bounds=([-np.inf, 1e-4, 0.0], [np.inf, np.inf, np.inf]),
                maxfev=10000
            )
            initial_guess.append(popt_sub)
        except (RuntimeError, ValueError):
            continue

    # Fallback if segment fits failed
    if not initial_guess:
        area = max(0.0, np.trapezoid(y, x=np.log10(x)))
        mu0 = float(np.nanmean(np.log(x)))
        sigma0 = float(np.nanstd(np.log(x)))
        sigma0 = sigma0 if np.isfinite(sigma0) and sigma0 > 1e-3 else 0.3
        return np.array([mu0, sigma0, area])
    
    return np.concatenate(initial_guess)

def make_constrained_lognorm(n_peaks: int) -> Callable:
    """Factory function to build a constrained lognormal wrapper for curve_fit."""
    def constrained_lognorm(x_data, *params):
        mu1, sigma1, scale1 = params[0:3]
        actual_params = [mu1, sigma1, scale1]
        
        for i in range(1, n_peaks):
            mu_i = params[i*3]
            delta_i = params[i*3 + 1]
            scale_i = params[i*3 + 2]
            
            # Enforce the relationship: sigma_i = sigma1 + delta_i
            sigma_i = sigma1 + delta_i  
            actual_params.extend([mu_i, sigma_i, scale_i])
            
        return multi_lognormal_distribution(x_data, *actual_params)
    return constrained_lognorm

def perform_constrained_fit(x: np.ndarray, y: np.ndarray, p0: np.ndarray) -> np.ndarray:
    """Executes a Poisson-weighted, multi-lognormal fit with relative GSD bounds."""
    n_peaks = len(p0) // 3
    fit_wrapper = make_constrained_lognorm(n_peaks)
    
    p0_constrained = []
    lower_bounds = []
    upper_bounds = []

    for i in range(n_peaks):
        mu_guess, sigma_guess, scale_guess = p0[i*3 : i*3+3]

        # Hard-clip incoming guesses to ensure they strictly respect bounds
        sigma_guess = max(1e-4, sigma_guess)
        scale_guess = max(0.0, scale_guess)
        
        if i == 0:
            p0_constrained.extend([mu_guess, sigma_guess, scale_guess])
            # Peak 1 Bounds: mu > -inf, sigma > 1e-5 (must be positive), scale > 0
            lower_bounds.extend([-np.inf, 1e-5, 0])
            upper_bounds.extend([np.inf, np.inf, np.inf])
        else:
            # For subsequent peaks, guess delta = 0 (assume identical GSD initially)
            delta_guess = 0.0 
            p0_constrained.extend([mu_guess, delta_guess, scale_guess])
            
            # Bounds for delta enforce the +/- 4% GSD constraint
            lower_bounds.extend([-np.inf, np.log(0.96), 0])
            upper_bounds.extend([np.inf, np.log(1.04), np.inf])

    popt_constrained, _ = curve_fit(
        fit_wrapper, 
        x, y, 
        p0=p0_constrained,
        bounds=(lower_bounds, upper_bounds),
        absolute_sigma=True,
        maxfev=10000
    )

    # Unpack the constrained parameters back into standard [mu, sigma, scale] format
    popt = []
    mu1_fit, sigma1_fit, scale1_fit = popt_constrained[0:3]
    popt.extend([mu1_fit, sigma1_fit, scale1_fit])
    
    for i in range(1, n_peaks):
        mu_i_fit = popt_constrained[i*3]
        delta_i_fit = popt_constrained[i*3 + 1]
        scale_i_fit = popt_constrained[i*3 + 2]
        popt.extend([mu_i_fit, sigma1_fit + delta_i_fit, scale_i_fit])
        
    return np.array(popt)

    
# =============================================================
# Interactive Plot Selection
# =============================================================
class InteractivePeakSelector:
    """Class to manage interactive peak selection using screen-pixel distance."""
    
    def __init__(self, fig, ax, x_data: np.ndarray, y_data: np.ndarray, title: str, x_label: str):
        self.fig = fig
        self.ax = ax
        self.x = x_data
        self.y = y_data
        
        self.sel_max = []
        self.sel_min = []
        self.idx_max = []
        self.idx_min = []
        
        self.title = title
        self.x_label = x_label
        self.cid = self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        self._update_plot()
        
    def _find_nearest_pixel(self, px_click: float, py_click: float) -> Tuple[int, float, float, float]:
        """Find the nearest data point based on display (pixel) distance."""
        # Convert data coordinates to pixel display coordinates
        pts = np.column_stack((self.x, self.y))
        pts_px = self.ax.transData.transform(pts)
        
        # Calculate Euclidean distance in pixel space
        distances = np.sqrt((pts_px[:, 0] - px_click)**2 + (pts_px[:, 1] - py_click)**2)
        idx = int(np.argmin(distances))
        return idx, self.x[idx], self.y[idx], distances[idx]

    def _update_plot(self):
        self.ax.cla()
        self.ax.set_xscale('log')
        self.ax.plot(self.x, self.y, 'o-', markersize=5, color='gray', label='Data')
        
        if self.sel_max:
            sx, sy = zip(*self.sel_max)
            self.ax.plot(sx, sy, 'ro', markersize=8, label='Max Points')
        if self.sel_min:
            sx2, sy2 = zip(*self.sel_min)
            self.ax.plot(sx2, sy2, 'go', markersize=8, label='Min Points')
            
        self.ax.set_xlabel(self.x_label)
        self.ax.set_ylabel('Concentration (#/cc)')
        self.ax.set_title(self.title + "\n(L-Click: Max | R-Click: Min | Click near point to remove)")
        self.ax.legend()
        self.ax.grid(True)
        self.fig.canvas.draw()

    def on_click(self, event):
        # Ensure click is inside the axes and has pixel coordinates
        if event.inaxes != self.ax: return
        if event.x is None or event.y is None: return

        idx, X, Y, dist_px = self._find_nearest_pixel(event.x, event.y)
        
        # Pixel threshold to remove an existing point
        rem_thresh = 15.0 
        removed = False

        # Attempt to remove existing Max point
        if self.sel_max and dist_px < rem_thresh:
            for i, (mx, my) in enumerate(self.sel_max):
                if mx == X and my == Y: # Data coordinates will match exactly because they snap to points
                    self.sel_max.pop(i)
                    self.idx_max.pop(i)
                    print(f"Removed MAX Point: (X={X:.2f}, Y={Y:.2f})")
                    removed = True
                    break
                    
        # Attempt to remove existing Min point
        if not removed and self.sel_min and dist_px < rem_thresh:
            for i, (mx, my) in enumerate(self.sel_min):
                if mx == X and my == Y:
                    self.sel_min.pop(i)
                    self.idx_min.pop(i)
                    print(f"Removed MIN Point: (X={X:.2f}, Y={Y:.2f})")
                    removed = True
                    break

        # If we didn't remove anything, add a new point based on the mouse button
        if not removed:
            if event.button == 1:  # Left click -> Max Point
                if idx not in self.idx_max and idx not in self.idx_min:
                    self.sel_max.append((X, Y))
                    self.idx_max.append(idx)
                    print(f"Selected MAX Point: (X={X:.2f}, Y={Y:.2f}, idx={idx})")
            elif event.button == 3:  # Right click -> Min Point
                if idx not in self.idx_min and idx not in self.idx_max:
                    self.sel_min.append((X, Y))
                    self.idx_min.append(idx)
                    print(f"Selected MIN Point: (X={X:.2f}, Y={Y:.2f}, idx={idx})")
                
        self._update_plot()

# =============================================================
# Main Analysis Routine
# =============================================================
def process_tandem_data(file_path: str):
    """Main script execution: read arrays natively, fit distributions, and save results."""
    
    # EXTRACT DATA FROM INPUTS DICTIONARY
    inputs = tandem_data_utils_2.main(file_path)
    df = inputs["df"]
    
    metadata = inputs["metadata"]
    qa_1 = metadata.get("qa_1", None)
    c1_name = metadata.get("classifier_1", "").lower()
    c2_name = metadata.get("classifier_2", "").lower()

    arrays = inputs["arrays"]
    cls1_data = arrays.get("cls1", {})
    cls2_data = arrays.get("cls2", {})
    bypass_cls2_data = arrays.get("bypass_cls2", {})

    is_c1_cpma = 'cpma' in c1_name
    is_c1_aac  = 'aac' in c1_name
    is_c1_dma = 'dma' in c1_name
    
    is_c2_cpma = 'cpma' in c2_name
    is_c2_aac  = 'aac' in c2_name
    is_c2_dma = 'dma' in c2_name

    cpma_mp = cpma_spd = cpma_v = cpma_T = cpma_p = np.array([], dtype=float)
    aac_da = aac_spd = aac_qsh = aac_T = aac_p = np.array([], dtype=float)
    dma_dm = dma_T = dma_p = np.array([], dtype=float)

    # Route Classifer 1
    if is_c1_cpma:
        cpma_mp = cls1_data.get("size", np.array([]))
        cpma_spd = cls1_data.get("speed", np.array([]))
        cpma_v = cls1_data.get("voltage", np.array([]))
        cpma_T = cls1_data.get("temperature", np.array([]))
        cpma_p = cls1_data.get("pressure", np.array([]))
    elif is_c1_aac:
        aac_da = cls1_data.get("size", np.array([]))
        aac_spd = cls1_data.get("speed", np.array([]))
        aac_qsh = cls1_data.get("sheath_flow", np.array([]))
        aac_T = cls1_data.get("temperature", np.array([]))
        aac_p = cls1_data.get("pressure", np.array([]))
    elif is_c1_dma:
        dma_dm = cls1_data.get("size", np.array([]))
        dma_T = cls1_data.get("temperature", np.array([]))
        dma_p = cls1_data.get("pressure", np.array([]))

    # Route Classifier 2
    if is_c2_cpma:
        cpma_mp = cls2_data.get("size", np.array([]))
        cpma_spd = cls2_data.get("speed", np.array([]))
        cpma_v = cls2_data.get("voltage", np.array([]))
        cpma_T = cls2_data.get("temperature", np.array([]))
        cpma_p = cls2_data.get("pressure", np.array([]))
    elif is_c2_aac:
        aac_da = cls2_data.get("size", np.array([]))
        aac_spd = cls2_data.get("speed", np.array([]))
        aac_qsh = cls2_data.get("sheath_flow", np.array([]))
        aac_T = cls2_data.get("temperature", np.array([]))
        aac_p = cls2_data.get("pressure", np.array([]))
    elif is_c2_dma:
        dma_dm = cls2_data.get("size", np.array([]))
        dma_T = cls2_data.get("temperature", np.array([]))
        dma_p = cls2_data.get("pressure", np.array([]))
    
    conc = arrays.get("conc", np.array([]))

    def _align_for_mask(arr: np.ndarray, target_len: int, fill_value: float = np.nan) -> np.ndarray:
        """Return array aligned to target_len or filled if missing."""
        arr = np.asarray(arr, dtype=float)
        if arr.size == target_len:
            return arr
        if arr.size == 0:
            return np.full(target_len, fill_value, dtype=float)
        raise ValueError(
            f"Array length mismatch: expected {target_len}, got {arr.size}. "
            "Check parsed columns in tandem_data_utils.main()."
        )

    _, cpma_m_star_arr = generate_cpma_setpoints(cpma_spd, cpma_v)

    if is_c1_aac:
        _, aac_da_star_arr_c1 = generate_aac_setpoints(aac_spd = aac_spd, aac_qsh = aac_qsh, aac_T = aac_T, aac_P = aac_p)

    if is_c2_aac:
        _, aac_da_star_arr_c2 = generate_aac_setpoints(aac_da=aac_da, aac_spd=aac_spd, aac_qsh=aac_qsh, aac_T=aac_T, aac_P=aac_p)

    if is_c1_cpma and cpma_mp.size > 0:
        setpoint_array, first_dim_label, dim_units = cpma_mp, 'm*', 'fg'
    elif is_c1_aac and aac_da.size > 0:
        setpoint_array, first_dim_label, dim_units = aac_da, 'da*', 'nm'
    elif is_c1_dma and dma_dm.size > 0:
        setpoint_array, first_dim_label, dim_units = dma_dm, 'dm*', 'nm'
    else:
        setpoint_array = cpma_mp if cpma_mp.size > 0 else (aac_da if aac_da.size > 0 else dma_dm)
        first_dim_label, dim_units = 'Setpoint', ''

    unique_setpoints = np.unique(setpoint_array)
    unique_setpoints = unique_setpoints[np.isfinite(unique_setpoints)]
    results = []

    target_len = setpoint_array.size
    if conc.size != target_len:
        raise ValueError(f"Concentration length ({conc.size}) does not match setpoint length ({target_len}).")
        
    if is_c2_cpma:
        raw_x = cpma_mp
        raw_p = cpma_p
        raw_T = cpma_T
    elif is_c2_aac:
        raw_x = aac_da
        raw_p = aac_p
        raw_T = aac_T
    elif is_c2_dma:
        raw_x = dma_dm
        raw_p = dma_p
        raw_T = dma_T
    else:
        raw_x = arrays.get("bypass_da", np.array([]))
        raw_p = np.array([])
        raw_T = np.array([])

    x_axis_arr = _align_for_mask(raw_x, target_len)
    pressure = _align_for_mask(raw_p, target_len, fill_value=1)
    Temp = _align_for_mask(raw_T, target_len, fill_value=293.15)

    # 2. Process each unique setpoint
    for sp in unique_setpoints:
        mask = (setpoint_array == sp)
        
        if is_c1_cpma:
            masked_m_star = cpma_m_star_arr[mask]
            c1_val, c1_label = (np.nanmean(masked_m_star) if masked_m_star.size > 0 else np.nan), 'm* (kg)'
        elif is_c1_aac:
            masked_da_star = aac_da_star_arr_c1[mask]
            c1_val, c1_label = (np.nanmean(masked_da_star) if masked_da_star.size > 0 else np.nan), 'da* (nm)'
        elif is_c1_dma:
            c1_val, c1_label = float(sp), 'dm* setpoint (nm)'
        else:
            c1_val, c1_label = sp, 'Setpoint_Val'

        x_raw, y_raw = x_axis_arr[mask], conc[mask]
        
        # Stricter NaN/Zero handling for curve fitting and plotting
        valid_mask = np.isfinite(x_raw) & np.isfinite(y_raw) & (x_raw > 0) & (y_raw >= 0)
        x, y = x_raw[valid_mask], y_raw[valid_mask]
        
        if len(x) == 0:
            print(f"Skipping setpoint {sp:.2f} {dim_units} due to lack of valid data.")
            continue
            
        p_mean = np.nanmean(pressure[mask][valid_mask]) if pressure[mask][valid_mask].size > 0 else 1.0
        t_mean = np.nanmean(Temp[mask][valid_mask]) if Temp[mask][valid_mask].size > 0 else 293.15
        
        title = f'{first_dim_label} = {sp:.2f} {dim_units}'

        if is_c2_dma:
            x_label_plot = 'Mobility Diameter setpoint (nm)'
        elif is_c2_aac:
            x_label_plot = 'Aerodynamic diameter setpoint (nm)'
        else:
            x_label_plot = 'Single particle mass setpoint (fg)'

        # ------------------------------------------------------------------
        # INNER LOOP: Allows retrying the fit for the current setpoint
        # ------------------------------------------------------------------
        popt = None  # Reset popt for each setpoint
        
        while True:
            # 3. Interactive Selection Phase
            fig, ax = plt.subplots()
            print('\n-----------------------------------------------------------')
            print(f"Selecting peaks for: {title}")
            print(" - Left-click to add a MAX peak.")
            print(" - Right-click to add a MIN valley.")
            print(" - Click near an existing point to remove it.")
            print(" - Close the plot window when finished to process fit.")
            print('-----------------------------------------------------------')
            
            selector = InteractivePeakSelector(fig, ax, x, y, title, x_label_plot)
            plt.show() # Automatically blocks execution until window is closed

            # 4. Fitting Phase
            p0 = generate_initial_guesses(x, y, selector)

            try:
                popt_temp = perform_constrained_fit(x, y, p0)
            except RuntimeError:
                print(f"\n[WARNING] Curve fit failed for setpoint {sp}.")
                user_choice = input("Press 'r' to Retry, or 's' to Skip this setpoint: ").strip().lower()
                if user_choice == 'r':
                    continue
                else:
                    break # Breaks the while loop, popt remains None
            
            # 5. Plot Results
            x_fit = np.logspace(np.log10(np.nanmin(x)), np.log10(np.nanmax(x)), 400)
            y_fit_total = multi_lognormal_distribution(x_fit, *popt_temp)
            
            plt.figure()
            plt.xscale('log')
            plt.plot(x, y, 'o', markersize=2,  color='gray', alpha=0.7, label='Data')
            plt.plot(x_fit, y_fit_total, '-', color='black', linewidth=1, label='Total Fitted Distribution')
            
            cmap = plt.get_cmap('tab10')
            for i in range(0, len(popt_temp), 3):
                peak_params = popt_temp[i:i+3]
                y_individual = multi_lognormal_distribution(x_fit, *peak_params)
                peak_num = (i // 3) + 1
                plt.plot(x_fit, y_individual, '--', color=cmap((peak_num - 1) % 10), label=f'Peak {peak_num}')
            
            plt.xlabel(x_label_plot)
            plt.ylabel('CPC Concentration (#/cc)')
            plt.title(f'Fit: {title}')
            plt.legend()
            plt.grid(True, linestyle=':', alpha=0.6, zorder=0)

            # Show the fit to the user, but do not save the figure yet
            plt.show()

            # PROMPT USER FOR APPROVAL
            print("\nDid the fit look acceptable?")
            user_choice = input("Press [Enter] to ACCEPT, 'r' to RETRY, or 's' to SKIP: ").strip().lower()
            
            if user_choice == 'r':
                print("Retrying fit...")
                continue
            elif user_choice == 's':
                print("Skipping setpoint...")
                break # popt remains None
            else:
                # Accepted! Save the popt and the figure
                popt = popt_temp
                plot_path = 'multi_lognormal_distribution'
                os.makedirs(plot_path, exist_ok=True)
                
                # Re-create and save the figure invisibly so we don't bother the user with it twice
                fig_save = plt.figure()
                plt.xscale('log')
                plt.plot(x, y, 'o', markersize=2, color='gray', alpha=0.7, label='Data')
                plt.plot(x_fit, y_fit_total, '-', color='black', linewidth=1, label='Total Fitted Distribution')
                for i in range(0, len(popt), 3):
                    peak_params = popt[i:i+3]
                    y_individual = multi_lognormal_distribution(x_fit, *peak_params)
                    peak_num = (i // 3) + 1
                    plt.plot(x_fit, y_individual, '--', color=cmap((peak_num - 1) % 10), label=f'Peak {peak_num}')
                plt.xlabel(x_label_plot)
                plt.ylabel('CPC Concentration (#/cc)')
                plt.title(f'Fit: {title}')
                plt.legend()
                plt.grid(True, linestyle=':', alpha=0.6, zorder=0)
                plt.savefig(os.path.join(plot_path, f"{title}.png"), dpi=300, bbox_inches='tight')
                plt.close(fig_save)
                
                break # Exit the while loop and proceed to aggregation

        # If the user chose to skip or the fit failed and was skipped, move to the next setpoint
        if popt is None:
            continue

        # 6. Aggregate Results
        params = popt.reshape(-1, 3)
        n_modes = params.shape[0]

        for i, (mu, sigma, N_total) in enumerate(params):
            median_m = np.exp(mu) * 1e-9
            GSD = np.exp(sigma)
            charge = (n_modes - i) if is_c2_dma else (i + 1)
            print('CMD:', median_m)

            if is_c2_aac:
                median_val, median_label = median_m, 'CMD (m)'
            elif is_c2_dma:
                _, median_val = compute_mobility_diameter(median_m, p_mean, t_mean, charge)
                median_label = 'CMD (m)'
            else:
                median_val, median_label = median_m, 'CMD (kg)'

            row = {
                'Setpoint_Orig': sp,
                c1_label: c1_val,
                'Charge': charge,
                median_label: median_val,
                'GSD': GSD,
                'N_total': N_total,
            }
            if is_c1_cpma and is_c2_aac:
                mobility_diam = float(mp_da2dm(c1_val * charge, median_val, t_mean, p_mean))
                if np.isfinite(median_val) and np.isfinite(c1_val) and c1_val != 0.0:
                    row['mobility diameter (m)'] = mobility_diam
                    row['effective_density (kg/m^3)'] = 6 * c1_val * charge / (np.pi * mobility_diam**3)
                else:
                    row['mobility diameter (m)'] = np.nan
                    row['effective_density (kg/m^3)'] = np.nan
            elif is_c1_cpma and is_c2_dma:
                mobility_diam = median_val
                if np.isfinite(median_val) and np.isfinite(c1_val) and c1_val != 0.0:
                    row['effective_density (kg/m^3)'] = 6 * c1_val * charge / (np.pi * mobility_diam**3)
                else:
                    row['effective_density (kg/m^3)'] = np.nan
            elif is_c1_aac and is_c2_dma:
                mobility_diam = median_val
                mass_p = float(dm_da2mp(mobility_diam, c1_val*1e-9, t_mean, p_mean))
                if np.isfinite(median_val) and np.isfinite(c1_val) and c1_val != 0.0:
                    row['particle_mass (kg)'] = mass_p
                    row['effective_density (kg/m^3)'] = 6 * mass_p / (np.pi * mobility_diam**3)
                else:
                    row['particle_mass (kg)'] = np.nan
                    row['effective_density (kg/m^3)'] = np.nan

            results.append(row)

    # 7. Save to CSV
    current_time = datetime.now().strftime('%Y-%m-%d_%H-%M')
    output_csv = f'results_{current_time}.csv'
    
    if results:
        with open(output_csv, mode='w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        print(f"\nResults successfully saved to: {output_csv}")
    else:
        print("\nNo results generated to save.")

    # 8. Plot Effective Density vs. Median Diameter
    d_m_c1, eff_rho_c1, d_m_c2, eff_rho_c2 = [], [], [], []

    for row in results:
        charge = row.get('Charge')
        density = row.get('effective_density (kg/m^3)', np.nan)
        if is_c2_dma:
            diam = row.get('CMD (m)', np.nan)
            x_label = 'Mobility diameter setpoint (nm)'
        elif is_c2_aac:
            diam = row.get('CMD (m)', np.nan)
            x_label = 'Aerodynamic diameter setpoint (nm)'
        else:
            print("\nNo effective density data available to plot.")

        if not np.isnan(density) and not np.isnan(diam):
            diam_nm = diam * 1e9 if diam < 1e-3 else diam 
            if charge == 1:
                d_m_c1.append(diam_nm)
                eff_rho_c1.append(density)
            elif charge == 2:
                d_m_c2.append(diam_nm)
                eff_rho_c2.append(density)

    if d_m_c1 or d_m_c2:     
        fig, ax1 = plt.subplots(figsize=(10, 6)) 
        if d_m_c1:
            ax1.plot(d_m_c1, eff_rho_c1, 'o', color='tab:blue', markersize=8, markeredgecolor='black', alpha=0.8, label='Singly_charged')
        if d_m_c2:
            ax1.plot(d_m_c2, eff_rho_c2, 's', color='tab:orange', markersize=8, markeredgecolor='black', alpha=0.8, label='Doubly_charged')
            
        ax1.set_xlabel(fr'{x_label}', fontweight='bold')
        ax1.set_ylabel(r'Effective Density (kg/m$^3$)', fontweight='bold')
        ax1.set_xscale('log')
        ax1.grid(True, linestyle=':', alpha=0.7)
        ax1.legend(loc='best', framealpha=0.9) 

        plt.tight_layout()
        os.makedirs(plot_path, exist_ok=True)
        plt.savefig(os.path.join(plot_path, f"Effective_Density_{current_time}.png"), dpi=300, bbox_inches='tight')
        plt.show()
    else:
        print("\nNo effective density data available to plot.")

# =============================================================
# Entry Point
# =============================================================
if __name__ == "__main__":
    FILE_PATH = r"Medresolution-1.csv"
    process_tandem_data(FILE_PATH)
# %%
