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

phases = [0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0, 195.0, 210.0, 225.0, 240.0, 255.0, 270.0, 285.0, 300.0, 315.0, 330.0, 345.0]

STEP_ONE = False
STEP_TWO = True
STEP_THREE = False

# Count the number of true statements
true_count = sum([STEP_ONE, STEP_TWO, STEP_THREE])

# Check if more than one statement is true
if true_count > 1:
    # Code to be executed if more than one statement is true
    print("You can't have more than one step be true")
    exit(0)
else:
    pass

# For the first step you're only generating the non-rotated stuff
# So you don't need to have all the jobs send off for the different phases
# If you do, they'll save over each other
if (STEP_ONE == True):
    phases = [0.0]

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
        

        if STEP_ONE == True:
            # Replace the phase in the run_spectra file with the correct one
            for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + ".py"], inplace=True):
                if line.strip().startswith('STEP_ONE'):
                    line = '    STEP_ONE = True\n'
                sys.stdout.write(line)

            for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + ".py"], inplace=True):
                if line.strip().startswith('STEP_TWO'):
                    line = '    STEP_TWO = False\n'
                sys.stdout.write(line)


            for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + ".py"], inplace=True):
                if line.strip().startswith('STEP_THREE'):
                    line = '    STEP_THREE = False\n'
                sys.stdout.write(line)


        
        if STEP_TWO == True:
            # Replace the phase in the run_spectra file with the correct one
            for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + ".py"], inplace=True):
                if line.strip().startswith('STEP_ONE'):
                    line = '    STEP_ONE = False\n'
                sys.stdout.write(line)

            for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + ".py"], inplace=True):
                if line.strip().startswith('STEP_TWO'):
                    line = '    STEP_TWO = True\n'
                sys.stdout.write(line)

            for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + ".py"], inplace=True):
                if line.strip().startswith('STEP_THREE'):
                    line = '    STEP_THREE = False\n'
                sys.stdout.write(line)
        

        if STEP_THREE == True:
            for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + ".py"], inplace=True):
                if line.strip().startswith('STEP_ONE'):
                    line = '    STEP_ONE = False\n'
                sys.stdout.write(line)

            for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + ".py"], inplace=True):
                if line.strip().startswith('STEP_TWO'):
                    line = '    STEP_TWO = False\n'
                sys.stdout.write(line)


            for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + ".py"], inplace=True):
                if line.strip().startswith('STEP_THREE'):
                    line = '    STEP_THREE = True\n'
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
            os.system(sbatch_command)

            print("Runnng", new_file_name)


            #python_command = f"python run_spectra_{finished_gcms[i]}_{str(phase)}.py"
            #os.system(python_command)
        
        if STEP_THREE == True:
            # Wait until the file is created before proceeding to the next iteration
            executable = f"Spectra/rt_emission_aerosols_{finished_gcms[i]}_phase_{str(phase)}.exe"

            while not os.path.exists(executable):
                print("You've gotta wait for the jobs to be executed in series unfortunately")
                time.sleep(1)