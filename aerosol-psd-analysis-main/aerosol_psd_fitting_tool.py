"""
==============================================================
CPMA–SMPS or AAC-SMPS Data Fitting and Multi-Lognormal Analysis Tool
==============================================================
Author: Morteza Kiasadegh
Date: 2025-10-16

Description:
------------
This script performs multi-lognormal fitting of particle size
distribution data obtained from CPMA–SMPS or AAC-SMPS experiments.
It allows manual selection of peaks using mouse clicks, performs
curve fitting, and calculates the equivalent mobility diameters
for multiply charged particles.

References:
-----------
- Kim, J. H., Mulholland, G. W., Kukuck, S. R., & Pui, D. Y. (2005)
  "Slip Correction Measurements of Certified PSL Nanoparticles..."
- Rader, D. J. (1990)
  "Momentum Slip Correction Factor for Small Particles in Nine Common Gases"
"""

# =============================================================
# Imports
# =============================================================
import numpy as np
import pandas as pd
import csv
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.optimize import curve_fit, fsolve
from datetime import datetime


# =============================================================
# Physical Helper Functions
# =============================================================

def cunningham_correction_factor(diameter_nm, pressure_atm, temperature_K):
    """
    Compute the Cunningham slip correction factor.

    Parameters
    ----------
    diameter_nm : float
        Particle diameter in nanometers.
    pressure_atm : float
        Ambient pressure in atmospheres.
    temperature_K : float
        Gas temperature in Kelvin.

    Returns
    -------
    float
        Cunningham slip correction factor.
    """
    # Constants from Kim et al. (2005)
    alpha = 1.165 * 2
    beta = 0.483 * 2
    gamma = 0.997 / 2

    # Mean free path at standard conditions
    lambda_air_std = 67.30e-9

    # Adjust mean free path for pressure and temperature
    lambda_air = lambda_air_std * (temperature_K / 296.15) ** 2 * (1 / pressure_atm) \
        * ((110.4 + 296.15) / (temperature_K + 110.4))

    d_m = diameter_nm * 1e-9  # convert nm → m

    return 1 + (lambda_air / d_m) * (alpha + beta * np.exp(-gamma * d_m / lambda_air))


def compute_mobility_diameter(d_m_initial, pressure_atm, temperature_K, charge_number):
    """
    Compute electrical mobility and the corresponding diameter for a given charge number.

    Parameters
    ----------
    d_m_initial : float
        Initial guess for mobility diameter (nm).
    pressure_atm : float
        Pressure in atm.
    temperature_K : float
        Temperature in Kelvin.
    charge_number : int
        Number of elementary charges.

    Returns
    -------
    z : float
        Electrical mobility (m²/V·s).
    d_m_solution : float
        Solved mobility diameter (nm).
    """
    e = 1.602e-19  # elementary charge (C)

    # Air viscosity (Rader, 1990)
    mu = 1.81809e-5 * (temperature_K / 293.15) ** 1.5 * \
        (293.15 + 110.4) / (temperature_K + 110.4)

    # Electrical mobility
    z = e * cunningham_correction_factor(d_m_initial, pressure_atm, temperature_K) / \
        (3 * np.pi * mu * d_m_initial)

    def equation(d_m_val):
        return d_m_val - (charge_number * e *
                          cunningham_correction_factor(d_m_val, pressure_atm, temperature_K)
                          ) / (3 * np.pi * mu * z)

    d_m_solution = fsolve(equation, charge_number * d_m_initial)
    return z, d_m_solution


# =============================================================
# Curve Fitting Functions
# =============================================================

def multi_lognormal_distribution(x, *params):
    """
    Sum of multiple lognormal distributions.

    Parameters
    ----------
    x : array-like
        Independent variable (mobility diameter).
    params : sequence
        Parameters in triplets (mu, sigma, scale).

    Returns
    -------
    np.ndarray
        Combined lognormal distribution.
    """
    result = np.zeros_like(x, dtype=float)
    for i in range(0, len(params), 3):
        mu, sigma, scale = params[i:i + 3]
        result += scale * np.exp(-0.5 * ((np.log(x) - mu) / sigma) ** 2) / \
            (x * sigma * np.sqrt(2 * np.pi))
    return result


# =============================================================
# Interactive Plot Selection Functions
# =============================================================

def find_nearest_point(x_data, y_data, x_click, y_click):
    """Return the nearest data point to a clicked location."""
    distances = np.sqrt((x_data - x_click) ** 2 + (y_data - y_click) ** 2)
    nearest_index = np.argmin(distances)
    return nearest_index, x_data[nearest_index], y_data[nearest_index]


def on_click(event, x_data, y_data, num_max_points, num_min_points, 
             selected_max_points, selected_min_points, indices_max, indices_min, 
             update_callback, cid):
    """
    Handle mouse click events for selecting maximum and minimum points with indices.
    
    Left-click: select a point
    Right-click: deselect last point
    """
    if event.xdata is None or event.ydata is None:
        return

    # -------------------------
    # Selecting maximum points
    # -------------------------
    if len(selected_max_points) < num_max_points:
        
        if event.button == 1:  # Left-click
            idx, X, Y = find_nearest_point(x_data, y_data, event.xdata, event.ydata)
            selected_max_points.append((X, Y))
            indices_max.append(idx)
            print("Selected MAX Point (X, Y, idx):", (X, Y, idx))
        elif event.button == 3:  # Right-click
            if selected_max_points:
                removed = selected_max_points.pop()
                removed_idx = indices_max.pop()
                print("Removed MAX Point (X, Y, idx):", (removed[0], removed[1], removed_idx))

    # -------------------------
    # Selecting minimum points
    # -------------------------
    else:
        print('Now select the minimum peaks with a left-click (right-click to undo).')
        if event.button == 1:  # Left-click
            idx, X, Y = find_nearest_point(x_data, y_data, event.xdata, event.ydata)
            selected_min_points.append((X, Y))
            indices_min.append(idx)
            print("Selected MIN Point (X, Y, idx):", (X, Y, idx))
            if len(selected_min_points) >= num_min_points:
                plt.gcf().canvas.mpl_disconnect(cid)
                print("Selection completed. Callback disconnected.")
        elif event.button == 3:  # Right-click
            if selected_min_points:
                removed = selected_min_points.pop()
                removed_idx = indices_min.pop()
                print("Removed MIN Point (X, Y, idx):", (removed[0], removed[1], removed_idx))

    # Update plot after each selection
    update_callback()



def update_plot(x, y, selected_points, title_label):
    """Update the current plot with selected points."""
    plt.cla()
    plt.xscale('log')
    plt.plot(x, y, 'o-', label='Data')
    if selected_points:
        sx, sy = zip(*selected_points)
        plt.plot(sx, sy, 'ro', label='Selected Points')
    plt.xlabel('Mobility Diameter (nm)')
    plt.ylabel('Concentration (#/cc)')
    plt.title(title_label)
    plt.legend()
    plt.grid(True)
    plt.draw()


# =============================================================
# Main Analysis Routine
# =============================================================

def main():
    """Main script execution: read data, fit distributions, and save results."""
    # -------------------------------------------------------------
    # File Paths
    # -------------------------------------------------------------
    input_csv = r'CPMA_SMPS_LF_LR_DOS_Rep1.csv'
    current_time = datetime.now().strftime('%Y-%m-%d_%H-%M')
    output_csv = rf'results_{current_time}.csv'

    # -------------------------------------------------------------
    # Data Import and Preparation
    # -------------------------------------------------------------
    df = pd.read_csv(input_csv, skiprows=7)
    df.columns = df.columns.str.strip()

    df['Conc'] = pd.to_numeric(df['Conc'], errors='coerce')

    # Create sorted setpoints
    setpoints = sorted(pd.to_numeric(df.iloc[2:, 0].unique()))

    results = []
    
    # -------------------------------------------------------------
    # Main Loop: Process Each Setpoint
    # -------------------------------------------------------------
    for i, sp in enumerate(setpoints):
        try:
            mask = pd.to_numeric(df['Mp (fg)']) == sp
            first_dim = 'm*'
        except KeyError:
            mask = pd.to_numeric(df['Da (nm)']) == sp
            first_dim = 'da*'

        Dm = pd.to_numeric(df['Dm (nm)2'][mask])
        P = pd.to_numeric(df['Pressure (kPa)2'][mask]) / 101.325
        T = pd.to_numeric(df['Temperature (C)2'][mask]) + 273.15
        conc = pd.to_numeric(df['Conc'][mask])

        x, y = Dm.values, conc.values
        # Drop NaNs before any plotting/fitting
        valid_mask = ~(np.isnan(x) | np.isnan(y))
        x = x[valid_mask]
        y = y[valid_mask]

        plt.figure()
        plt.xscale('log')
        plt.plot(x, y, 'o-', label='Data')
        plt.title(f'{first_dim} = {sp:.2f} fg')
        plt.xlabel('Mobility Diameter (nm)')
        plt.ylabel('Concentration (#/cc)')
        plt.legend()
        plt.grid(True)
        plt.show()

        selected_max_points = []
        selected_min_points = []
        indices_max = []
        indices_min = []

        # User interactive selection
        num_max_points = int(input("How many maximum peaks would you like to select? "))
        print('Now select the maximum peaks with a left-click (right-click to undo).')
        num_min_points = num_max_points - 1


        # Define global indices for on_click
        global indices
        indices = [0]

        fig = plt.figure()
        plt.xscale('log')
        plt.plot(x, y, marker='o', label='Data')
        plt.legend()
        plt.grid(True)

        def update():
            # Redraw with selected maxima (red) and minima (green)
            plt.cla()
            plt.xscale('log')
            plt.plot(x, y, 'o-', label='Data')
            if selected_max_points:
                sx, sy = zip(*selected_max_points)
                plt.plot(sx, sy, 'ro', label='Max Points')
            if selected_min_points:
                sx2, sy2 = zip(*selected_min_points)
                plt.plot(sx2, sy2, 'go', label='Min Points')
            plt.xlabel('Mobility Diameter (nm)')
            plt.ylabel('Concentration (#/cc)')
            plt.title(f'{first_dim} = {sp:.2f} fg')
            plt.legend()
            plt.grid(True)
            plt.draw()

        # Single event connection with late-bound cid via holder
        cid_holder = {'cid': None}
        cid_holder['cid'] = plt.gcf().canvas.mpl_connect(
            'button_press_event',
            lambda event: on_click(
                event,
                x, y,
                num_max_points, num_min_points,
                selected_max_points, selected_min_points,
                indices_max, indices_min,
                update,
                cid_holder['cid']
            )
        )

        plt.show()

        # Fit multi-lognormal to selected peaks
        initial_guess = []
        indices = indices + indices_max + indices_min
        indices.append(len(y))
        indices.sort()

        for j in range(0, len(indices) - 2, 2):
            x_segment = x[indices[j]:indices[j+2] + 1]
            y_segment = y[indices[j]:indices[j+2] + 1]

            scale = np.trapezoid(y_segment, x_segment)
            guess = (np.log(np.mean(x_segment)), 0.2, scale)
            try:
                popt_sub, _ = curve_fit(multi_lognormal_distribution, x_segment, y_segment, p0=guess, maxfev=10000)
                initial_guess.append(popt_sub)
            except RuntimeError:
                continue


        # Fit full curve (guard against empty initial_guess)
        if len(initial_guess) == 0:
            # Fallback: single-component guess from full data
            area = np.trapezoid(y, x)
            mu0 = float(np.nanmean(np.log(x)))
            sigma0 = float(np.nanstd(np.log(x))) if np.isfinite(np.nanstd(np.log(x))) and np.nanstd(np.log(x)) > 0 else 0.3
            p0 = np.array([mu0, sigma0, area])
        else:
            p0 = np.concatenate(initial_guess)

        popt, _ = curve_fit(multi_lognormal_distribution, x, y, p0=p0, maxfev=10000)

        # Plot original data and fitted distribution
        x_fit = np.logspace(np.log10(np.nanmin(x)), np.log10(np.nanmax(x)), 400)
        y_fit = multi_lognormal_distribution(x_fit, *popt)
        plt.figure()
        plt.xscale('log')
        plt.plot(x, y, 'o', label='Data')
        plt.plot(x_fit, y_fit, '-', label='Fitted Distribution')
        plt.xlabel('Mobility Diameter (nm)')
        plt.ylabel('Concentration (#/cc)')
        plt.title(f'Fit: {first_dim} = {sp:.2f} fg')
        plt.legend()
        plt.grid(True)
        plt.show()

        # Store results
        for j in range(0, len(popt), 3):
            mu, sigma, N_total = popt[j:j + 3]
            median = np.exp(mu)
            charge = len(popt) // 3 - j // 3
            _, d_m_n = compute_mobility_diameter(median, P.mean(), T.mean(), charge)

            results.append({
                'Setpoint': sp,
                'Charge': charge,
                'Median_d_m': d_m_n[0],
                'N_total': N_total
            })

    # -------------------------------------------------------------
    # Save Results
    # -------------------------------------------------------------
    results_df = pd.DataFrame(results)
    results_df.to_csv(output_csv, index=False)
    print(f"\nResults saved to: {output_csv}")
    print(results_df)


# =============================================================
# Entry Point
# =============================================================
if __name__ == "__main__":
    main()
