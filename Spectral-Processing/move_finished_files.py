import os
import sys
import shutil
from shutil import copytree, ignore_patterns
import fileinput
import re

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

print("THIS HAS TO BE RUN IN PYTHON3, SORRY")
print ("All the files are run, and now you want to put everything is one folder?")
print("")
print("")

# Get the names of all the folders that were run
# Get the names of all the folders that you want to run
gcm_folder = input("Enter the name of the directory with all the GCMS that have been run: ")
finished_gcms = [name for name in os.listdir(gcm_folder) if os.path.isdir(os.path.join(gcm_folder, name))]

# Ask the user for some input
question = "Do you want to move " + str(finished_gcms) 
answer = query_yes_no(question)

if answer == True:
    for finished_gcm in finished_gcms:
        for filename in os.listdir('../Spectral-Processing-' + str(finished_gcm) +'/PLANET_MODELS'):
            if re.match(str(finished_gcm) + r'*', filename):
                print("Moving " + str(filename) + " out of " + 'Spectral-Processing-' + str(finished_gcm) +'/PLANET_MODELS/')
                shutil.move(os.path.join('../Spectral-Processing-' + str(finished_gcm) +'/PLANET_MODELS', filename), 'PLANET_MODELS')
            elif re.match(r'init_*', filename):
                print("Moving " + str(filename) + " out of " + 'Spectral-Processing-' + str(finished_gcm) +'/PLANET_MODELS/')
                shutil.move(os.path.join('../Spectral-Processing-' + str(finished_gcm) +'/PLANET_MODELS', filename), 'PLANET_MODELS')
            else:
                print("No models in this dir to move")

        for filename in os.listdir('../Spectral-Processing-' + str(finished_gcm) +'/FINISHED_SPECTRA'):
            if re.match(r'Spec_*', filename):
                print("Moving " + str(filename) + " out of " + 'Spectral-Processing-' + str(finished_gcm) +'/FINISHED_SPECTRA/')
                shutil.move(os.path.join('../Spectral-Processing-' + str(finished_gcm) +'/FINISHED_SPECTRA', filename), 'FINISHED_SPECTRA')
            else:
                print("No models in this dir to move")
else:
    print("Exit, you chose not to execute this program")


