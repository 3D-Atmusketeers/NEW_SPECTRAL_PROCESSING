import pandas as pd
import numpy as np

def find_closest_wavelength_indices(opacity_set_id, full_wavelength_range, WAVELENGTH_START_APPROX=None, WAVELENGTH_END_APPROX=None):
    """
    Finds the closest wavelength indices in a dataset for a specified range.

    Parameters:
    - opacity_set_id: Identifier for the dataset.
    - full_wavelength_range: Boolean indicating whether to use the full range of the dataset.
    - WAVELENGTH_START_APPROX: The approximate starting wavelength (used if full_wavelength_range is False).
    - WAVELENGTH_END_APPROX: The approximate ending wavelength (used if full_wavelength_range is False).

    Returns:
    - LAMBDA_START: Index of the start wavelength.
    - LAMBDA_END: Index of the end wavelength (inclusive).
    - START_WAVELENGTH: Value of the start wavelength.
    - END_WAVELENGTH: Value of the end wavelength.
    """
    # Construct the file path from the given opacity set ID
    file_path = f'DATA/{opacity_set_id}/Wavelengths.txt'

    # Read the first column as wavelengths from the file
    df = pd.read_csv(file_path, header=None, usecols=[0])

    # Get the number of lambda points in the wavelength files
    NLAMBDA = len(df)

    wavelengths = df.iloc[:, 0]

    # Check if the specified wavelength range is within the dataset's bounds
    min_wavelength, max_wavelength = wavelengths.min(), wavelengths.max()
    if not full_wavelength_range:
        if WAVELENGTH_START_APPROX < min_wavelength or WAVELENGTH_END_APPROX > max_wavelength:
            print(f"Error: Input wavelengths must be within the range {min_wavelength} - {max_wavelength}.")
            exit(0)

    # Determine the indices and values for the wavelength range
    if full_wavelength_range:
        # Use the entire range of wavelengths in the dataset
        LAMBDA_START, LAMBDA_END = 0, len(wavelengths) - 1
        START_WAVELENGTH, END_WAVELENGTH = wavelengths.iloc[LAMBDA_START], wavelengths.iloc[LAMBDA_END]
    else:
        # Find the closest wavelengths to the specified approximations
        abs_diff_start = np.abs(wavelengths - WAVELENGTH_START_APPROX)
        abs_diff_end = np.abs(wavelengths - WAVELENGTH_END_APPROX)
        
        LAMBDA_START = abs_diff_start.idxmin()
        LAMBDA_END = abs_diff_end.idxmin()
        
        START_WAVELENGTH = wavelengths.iloc[LAMBDA_START]
        END_WAVELENGTH = wavelengths.iloc[LAMBDA_END]

    # Return the calculated indices and wavelength values
    return LAMBDA_START, LAMBDA_END + 1, START_WAVELENGTH, END_WAVELENGTH, NLAMBDA
