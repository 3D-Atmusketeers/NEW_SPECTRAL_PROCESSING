import os
import sys
import shutil
from shutil import copytree, ignore_patterns
import fileinput


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


print("This is an important file, and has the capability of generating a lot of data in the wrong places")
print("")
print("Put this folder in the same folder as the directory of GCMS that you want to be run.")
print("This will copy Spectral-Processing folder several times, with one GCM in each folder and then run each of them")

print ("")
print ("")

# Get the names of all the folders that you want to run
gcm_folder = input("Enter the name of the directory with all the GCMS you want to run: ")
finished_gcms = [name for name in os.listdir(gcm_folder) if os.path.isdir(os.path.join(gcm_folder, name))]

# Ask the user for some input
question = "Do you want to run " + str(finished_gcms) 
answer = query_yes_no(question)

if answer == True:
    for i in range(len(finished_gcms)):
        # Make the top folder
        os.mkdir("../Spectral-Processing-" + finished_gcms[i])
        
        # Within each one make the GCM-OUTPUT and put in one single GCM
        shutil.copytree(gcm_folder + "/" + finished_gcms[i], "../Spectral-Processing-" + finished_gcms[i] + "/GCM-OUTPUT/" + finished_gcms[i])
        
        # Create empty folders
        os.mkdir("../Spectral-Processing-" + finished_gcms[i] + "/FINISHED_SPECTRA")
        os.mkdir("../Spectral-Processing-" + finished_gcms[i] + "/PLANET_MODELS")


        # Copy the spectra folder but ignore the DATA folder and the OUT folder
        copytree("Spectra", "../Spectral-Processing-" + finished_gcms[i] + "/Spectra/", ignore=ignore_patterns('init_*', 'Spec_*'))
        
        # Replace the planet name in the spectra folder with the correct one
        for line in fileinput.input(["../Spectral-Processing-" + finished_gcms[i] + "/Spectra/run_spectra.py"], inplace=True):
            if line.strip().startswith('planet_names'):
                line = 'planet_names = ["' + str(finished_gcms[i]) + '"]\n'
            sys.stdout.write(line)

        # Run the spectra for that folder
        print("Running " + "../Spectral-Processing-" + str(finished_gcms[i]) + " spectra")
        os.chdir("../Spectral-Processing-" + str(finished_gcms[i]) + "/Spectra/")        
        os.system("sbatch ./Run_sbatch")
        os.chdir("../../Spectral-Processing")

else:
    print("Exit, you chose not to execute this program")


