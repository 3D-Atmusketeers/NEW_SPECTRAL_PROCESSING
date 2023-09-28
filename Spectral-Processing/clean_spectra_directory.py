import os
import glob
import shutil

def delete_files_with_pattern(pattern):
    file_list = glob.glob(pattern)
    for path in file_list:
        try:
            if os.path.isfile(path):  # Check if it's a file
                os.remove(path)
                print(f"Deleted file {path}")
            elif os.path.isdir(path):  # Check if it's a directory
                shutil.rmtree(path)  # Use rmtree to remove the entire directory tree
                print(f"Deleted directory {path}")
        except OSError as e:
            print(f"Error: {e}")

def clean_spectra_directory(directory="Spectra/"):
    patterns = [
        "rt_emission_aerosols_*",
        "Run_sbatch_*",
        "run_spectra_*",
        "*fortfiles*",  # Match any path with the string "fortfiles" within it
        "*lock*"  # Match any path with the string "lock" within it
    ]

    for pattern in patterns:
        full_pattern = os.path.join(directory, pattern)
        delete_files_with_pattern(full_pattern)

if __name__ == "__main__":
    clean_spectra_directory()
