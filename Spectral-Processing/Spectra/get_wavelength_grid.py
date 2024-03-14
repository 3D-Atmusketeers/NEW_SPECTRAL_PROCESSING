import pandas as pd
import numpy as np

def find_closest_wavelength_indices(opacity_set_id, full_wavelength_range, WAVELENGTH_START_APPROX=None, WAVELENGTH_END_APPROX=None):
    file_path = f'DATA/{opacity_set_id}/Wavelengths.txt'

    # Assuming wavelengths are in the first column
    df = pd.read_csv(file_path, header=None, usecols=[0])
    wavelengths = df.iloc[:, 0]

    # Error handling for input wavelengths outside the possible range
    min_wavelength, max_wavelength = wavelengths.min(), wavelengths.max()
    if not full_wavelength_range:
        if WAVELENGTH_START_APPROX < min_wavelength or WAVELENGTH_END_APPROX > max_wavelength:
            print(f"Error: Input wavelengths must be within the range {min_wavelength} - {max_wavelength}.")
            exit(0)

    if full_wavelength_range:
        LAMBDA_START, LAMBDA_END = 0, len(wavelengths) - 1
        START_WAVELENGTH, END_WAVELENGTH = wavelengths.iloc[LAMBDA_START], wavelengths.iloc[LAMBDA_END]
    else:
        # Find the closest indices to the specified wavelengths
        abs_diff_start = np.abs(wavelengths - WAVELENGTH_START_APPROX)
        abs_diff_end = np.abs(wavelengths - WAVELENGTH_END_APPROX)
        
        LAMBDA_START = abs_diff_start.idxmin()
        LAMBDA_END = abs_diff_end.idxmin()
        
        START_WAVELENGTH = wavelengths.iloc[LAMBDA_START]
        END_WAVELENGTH = wavelengths.iloc[LAMBDA_END]

    return LAMBDA_START, LAMBDA_END + 1, START_WAVELENGTH, END_WAVELENGTH