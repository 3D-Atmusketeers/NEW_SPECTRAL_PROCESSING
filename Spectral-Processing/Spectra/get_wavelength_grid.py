import sys

def find_closest_wavelength_indices(opacity_set_id, wavelength_start, wavelength_end):
    """
    Finds the indices of the closest wavelengths to the specified start and end wavelengths
    in a data file. Exits the program with an error if the start is not >= the smallest wavelength
    or the end is not <= the largest wavelength in the file.

    Parameters:
    - opacity_set_id: The identifier for the set of opacity data.
    - wavelength_start: The starting wavelength to find the closest match for.
    - wavelength_end: The ending wavelength to find the closest match for.
    
    Returns:
    - A tuple containing the indices (NLAMBDA_START, NLAMBDA_END) of the closest matches
      to the wavelength_start and wavelength_end, respectively.
    """
    file_path = f'DATA/{opacity_set_id}/opacCIA.dat'
    
    # Initialize variables to track the first and last wavelength seen
    first_wavelength_seen, last_wavelength_seen = None, None
    closest_start, closest_end = float('inf'), float('inf')
    NLAMBDA_START, NLAMBDA_END = -1, -1
    index = -1  # Start with -1 to account for skipped header rows

    with open(file_path, 'r') as file:
        # Skip the first two header lines
        next(file), next(file)
        
        for line in file:
            index += 1
            wavelength = float(line.split()[0])
            
            if first_wavelength_seen is None:
                first_wavelength_seen = wavelength
            
            last_wavelength_seen = wavelength
            
            # Update if this wavelength is closer to the start value
            if abs(wavelength - wavelength_start) < abs(closest_start - wavelength_start):
                closest_start, NLAMBDA_START = wavelength, index
                
            # Update if this wavelength is closer to the end value
            if abs(wavelength - wavelength_end) < abs(closest_end - wavelength_end):
                closest_end, NLAMBDA_END = wavelength, index

    # After reading all lines, check the first and last seen wavelengths
    if first_wavelength_seen is None or last_wavelength_seen is None:
        print("Error: File appears to be empty or improperly formatted.")
        sys.exit(1)
    
    if wavelength_start < first_wavelength_seen or wavelength_end > last_wavelength_seen:
        print("Error: Specified wavelength range is out of bounds.")
        sys.exit(1)
    
    if NLAMBDA_START == -1 or NLAMBDA_END == -1:
        print("Error: Specified wavelength range did not match any data points.")
        sys.exit(1)

    return NLAMBDA_START, NLAMBDA_END
