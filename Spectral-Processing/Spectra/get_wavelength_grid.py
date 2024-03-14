import sys

def find_closest_wavelength_indices(opacity_set_id, WAVELENGTH_START=5.0e-6, WAVELENGTH_END=5.2e-6):
    print(f"Initializing search for closest wavelength indices in '{opacity_set_id}' dataset...")
    file_path = f'DATA/{opacity_set_id}/opacCIA.dat'

    closest_start_diff = float('inf')
    closest_end_diff = float('inf')
    closest_start_wavelength = None
    closest_end_wavelength = None
    NLAMBDA_START = None
    NLAMBDA_END = None
    index = -1

    with open(file_path, 'r') as file:
        print("File opened successfully. Processing data...")
        # Skip the first two header lines
        next(file), next(file)

        for line in file:
            index += 1
            wavelength = float(line.split()[0])

            # Calculate the difference for start wavelength
            if abs(wavelength - WAVELENGTH_START) < closest_start_diff:
                closest_start_diff = abs(wavelength - WAVELENGTH_START)
                NLAMBDA_START = index
                closest_start_wavelength = wavelength

            # Calculate the difference for end wavelength
            if abs(wavelength - WAVELENGTH_END) < closest_end_diff:
                closest_end_diff = abs(wavelength - WAVELENGTH_END)
                NLAMBDA_END = index
                closest_end_wavelength = wavelength

            # Stop condition based on the wavelength range
            if wavelength > WAVELENGTH_END:
                break

    if NLAMBDA_START is None or NLAMBDA_END is None:
        print("Error: Specified wavelength range did not match any data points.")
        sys.exit(1)

    print("\n" + "="*60)
    print(f"Start Index: {NLAMBDA_START}, End Index: {NLAMBDA_END}")
    print(f"Actual Start Wavelength: {closest_start_wavelength}, Actual End Wavelength: {closest_end_wavelength}")
    print("\n" + "="*60)
    return NLAMBDA_START, NLAMBDA_END, closest_start_wavelength, closest_end_wavelength
