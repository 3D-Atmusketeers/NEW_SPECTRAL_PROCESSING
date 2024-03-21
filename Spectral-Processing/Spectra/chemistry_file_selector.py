import glob
import math

def find_closest_chemistry_file(MET_X_SOLAR):
    # Calculate the log10 of MET_X_SOLAR
    log_MET_X_SOLAR = math.log10(MET_X_SOLAR)

    # Use glob to find all files that match the pattern in the DATA/chemistry_grid folder
    file_paths = glob.glob('DATA/chemistry_grid/fastchem_grid_allspecies_ions_lotemp_Z_*.dat')

    # Function to extract the Z value from a file path
    def extract_z_value(file_path):
        try:
            start = file_path.find("Z_") + 2
            end = file_path.find("_C_to_O", start)
            return float(file_path[start:end])
        except ValueError:
            # In case there's an issue converting to float, return None
            return None

    # Calculate the difference between log_MET_X_SOLAR and each file's Z value, then find the file with the minimum difference
    closest_file = None
    min_difference = float('inf')
    for file_path in file_paths:
        z_value = extract_z_value(file_path)
        if z_value is not None:
            difference = abs(log_MET_X_SOLAR - z_value)
            if difference < min_difference:
                min_difference = difference
                closest_file = file_path

    return closest_file
