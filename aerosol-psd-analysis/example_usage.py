#!/usr/bin/env python3
"""
Example usage of the CPMA-SMPS Multi-Lognormal Fitting Tool

This script demonstrates how to use the main fitting routine with a sample dataset.
"""

import os
import sys
from aerosol_psd_fitting_tool import main

def run_example():
    """
    Run the fitting routine with example data.
    
    Make sure you have a CSV file with the required columns:
    - Mp (fg) or Da (nm): Setpoint values
    - Dm (nm)2: Mobility diameter measurements  
    - Pressure (kPa)2: Pressure measurements
    - Temperature (C)2: Temperature measurements
    - Conc: Particle concentration
    """
    
    print("CPMA-SMPS Multi-Lognormal Fitting Tool")
    print("=" * 50)
    print()
    print("This tool will:")
    print("1. Load your particle size distribution data")
    print("2. Display each setpoint for interactive peak selection")
    print("3. Fit multi-lognormal distributions to selected peaks")
    print("4. Calculate mobility diameters for different charge states")
    print("5. Save results to a timestamped CSV file")
    print()
    print("Instructions:")
    print("- Left-click to select maximum peaks (red markers)")
    print("- Left-click to select minimum points between peaks (green markers)")
    print("- Right-click to undo the last selection")
    print("- The fitting will proceed automatically after selection")
    print()
    
    # Check if the data file exists
    data_file = "CPMA_SMPS_LF_LR_DOS_Rep1.csv"
    if not os.path.exists(data_file):
        print(f"Warning: Data file '{data_file}' not found.")
        print("Please ensure your CSV file is in the same directory as this script.")
        print("The file should contain the required columns as described in the README.")
        return
    
    try:
        # Run the main fitting routine
        main()
        print("\nAnalysis completed successfully!")
        print("Check the results_*.csv file for your fitted parameters.")
        
    except FileNotFoundError as e:
        print(f"Error: Could not find required file - {e}")
        print("Please check that your data file exists and has the correct format.")
        
    except Exception as e:
        print(f"An error occurred during analysis: {e}")
        print("Please check your data format and try again.")

if __name__ == "__main__":
    run_example()
