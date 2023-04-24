import os
import sys
import shutil
import fileinput
import contextlib
import time

@contextlib.contextmanager
def change_directory(path):
    """Context manager to temporarily change the working directory."""
    prev_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_dir)


print("MAKE SURE TO SET THE PHASES AND THE OPACITY SET IN RUN_SPECTRA.PY!!!")
print("")
print("This is an important file, and has the capability of generating a lot of data in the wrong places")
print("")
print("Put this folder in the same folder as the directory of GCMS that you want to be run.")
print("This will copy Spectral-Processing folder several times, with one GCM in each folder and then run each of them")

print ("")
print ("")


# Get the names of all the folders that you want to run
#gcm_folder = input("Enter the name of the directory with all the GCMS you want to run: ")
gcm_folder = 'GCM-OUTPUT'
finished_gcms = [name for name in os.listdir(gcm_folder) if os.path.isdir(os.path.join(gcm_folder, name))]

phases = [0.0]
RUN_REGRIDDING = False

source_file_name = "Run_sbatch"
for i in range(len(finished_gcms)):
    for j, phase in enumerate(phases):
        # copy the file run_spectra.py to a new file with the extension finished_gcms[i]
        shutil.copyfile("Spectra/run_spectra.py", "Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + ".py")

        # Replace the planet name in the spectra folder with the correct one
        for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + ".py"], inplace=True):
            if line.strip().startswith('planet_names'):
                line = 'planet_names = ["' + str(finished_gcms[i]) + '"]\n'
            sys.stdout.write(line)

        # Replace the phase in the spectra folder with the correct one
        for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + ".py"], inplace=True):
            if line.strip().startswith('phases'):
                line = 'phases = [' + str(phase) + ']\n'
            sys.stdout.write(line)
        
        if RUN_REGRIDDING == True:
            # Replace the phase in the run_spectra file with the correct one
            for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + ".py"], inplace=True):
                if line.strip().startswith('RUN_REGRIDDING'):
                    line = '    RUN_REGRIDDING = True\n'
                sys.stdout.write(line)
        else:
                # Replace the phase in the run_spectra file with the correct one
            for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + ".py"], inplace=True):
                if line.strip().startswith('RUN_REGRIDDING'):
                    line = '    RUN_REGRIDDING = False\n'
                sys.stdout.write(line)               
            

        # Run the spectra for that folder
        with change_directory("../Spectral-Processing/Spectra/"):
            new_file_name = "Run_sbatch_" + finished_gcms[i] + "_" + str(phase)
            shutil.copy(source_file_name, new_file_name)
            temp_file_name = f"{new_file_name}_temp"
            with open(new_file_name, "r") as file, open(temp_file_name, "w") as temp_file:
                for line in file:
                    if line.startswith("#SBATCH --job-name=spectra"):
                        line = f"#SBATCH --job-name={finished_gcms[i]}_{str(phase)}\n"
                    elif line.startswith("python run_spectra.py"):
                        line = f"python run_spectra_{finished_gcms[i]}_{str(phase)}.py\n"
                    temp_file.write(line)

            shutil.move(temp_file_name, new_file_name)
            sbatch_command = f"sbatch {new_file_name}"

            print(sbatch_command)
            #os.system(sbatch_command)

            python_command = f"python run_spectra_{finished_gcms[i]}_{str(phase)}.py"
            os.system(python_command)

        # Wait until the file is created before proceeding to the next iteration
        executable = f"Spectra/rt_emission_aerosols_{finished_gcms[i]}_phase_{str(phase)}.exe"

        while not os.path.exists(executable):
            print("You've gotta wait for the jobs to be executed in series unfortunately")
            time.sleep(1)