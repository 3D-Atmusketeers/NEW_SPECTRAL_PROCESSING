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


def query_yes_no(question, default="no"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
            It must be "yes" (the default), "no" or None (meaning
            an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True, "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == "":
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' " "(or 'y' or 'n').\n")

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

# Ask the user for some input
#question = "Do you want to run " + str(finished_gcms) 
#answer = query_yes_no(question)
answer = True

phases = [0.0]
source_file_name = "Run_sbatch"


if answer == True:
    for i in range(len(finished_gcms)):
        time.sleep(60)  
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
                time.sleep(60)  

else:
    print("Exit, you chose not to execute this program")