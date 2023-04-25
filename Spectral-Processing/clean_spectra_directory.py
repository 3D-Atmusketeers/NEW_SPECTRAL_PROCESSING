import os
import glob

def delete_files_with_pattern(pattern):
    file_list = glob.glob(pattern)
    for file in file_list:
        try:
            os.remove(file)
            print(f"Deleted {file}")
        except OSError as e:
            print(f"Error: {e}")

directory = "Spectra/"
patterns = ["rt_emission_aerosols_*", "Run_sbatch_*", "run_spectra_*"]

for pattern in patterns:
    full_pattern = os.path.join(directory, pattern)
    delete_files_with_pattern(full_pattern)