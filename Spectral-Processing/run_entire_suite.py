import os
import sys
import shutil
import fileinput
import contextlib
import time
import subprocess
import glob
import numpy as np

phases = [0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0, 120.0, 135.0, 150.0, 165.0, 180.0, 195.0, 210.0, 225.0, 240.0, 255.0, 270.0, 285.0, 300.0, 315.0, 330.0, 345.0]
gcm_folder = 'GCM-OUTPUT'
source_file_name = "Run_sbatch"

lowest_phase = np.min(phases)

@contextlib.contextmanager
def change_directory(path):
    """Context manager to temporarily change the working directory."""
    prev_dir = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_dir)

def runsbatch(phases, source_file_name, finished_gcm, step, dependency = 'none'):
    """ The function that submits jobs """
    jobnum = []
    for phase in enumerate(phases):
        with change_directory("../Spectral-Processing/Spectra/"):
            new_file_name = "Run_sbatch_" + finished_gcm + "_" + str(phase[1]) + '_' + step
            shutil.copy(source_file_name, new_file_name)
            temp_file_name = f"{new_file_name}_temp"
            with open(new_file_name, "r") as file, open(temp_file_name, "w") as temp_file:
                for line in file:
                    if line.startswith("#SBATCH --job-name=spectra"):
                        stepname = step.replace('step','')
                        print(stepname, 'stepname')
                        line = f"#SBATCH --job-name={str(phase[1])}_{stepname}_{finished_gcm}\n"
                    elif line.startswith("python run_spectra.py"):
                        line = f"python run_spectra_{finished_gcm}_{str(phase[1])}_{step}.py\n"
                    temp_file.write(line)

            shutil.move(temp_file_name, new_file_name)
            
            if dependency == 'none':
                sbatch = 'sbatch'
                name = new_file_name
                run = subprocess.run([sbatch, name], stdout = subprocess.PIPE)
            else:
                dependstring = ''
                for job in dependency:
                    dependstring+= ':' + job
                sbatch = 'sbatch'
                d = '-d'
                afterok = f'afterok{dependstring}'
                name = new_file_name
                run = subprocess.run([sbatch, d, afterok, name], stdout = subprocess.PIPE)
            job = run.stdout.decode('utf-8')
            print(job)
            jobnum.append(str(job[20:-1]).strip())
            print("Running", new_file_name)
    return jobnum


def copyfiles(phases, step, finished_gcms):
    print('phases are ', phases)   
    for phase in enumerate(phases):
        shutil.copyfile("Spectra/run_spectra.py", "Spectra/run_spectra_" + finished_gcms + "_" + str(phase[1]) + '_' + step + ".py")
        for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms + "_" + str(phase[1]) + '_' + step + ".py"], inplace=True):
            if line.strip().startswith('phases'):
                line = 'phases = [' + str(phase[1]) + ']\n'
            if line.strip().startswith('planet_names'):
                line = 'planet_names = ["' + str(finished_gcms) + '"]\n'
            if line.strip().startswith('Planet_name'):
                line = 'Planet_name = "'+ str(finished_gcms) + '"\n'
            sys.stdout.write(line)
            
    print('copied phase files successfully')
    return None


print("MAKE SURE TO SET THE OPACITY SET IN RUN_SPECTRA.PY!!!")
print("")
print("This will run the entire post processing suite")
print("")
print("Put this folder in the same folder as the directory of GCMS that you want to be run.")
print("This will copy Spectral-Processing folder several times, with one GCM in each folder and then run each of them")

print ("")
print ("")


# Get the names of all the folders that you want to run
#gcm_folder = input("Enter the name of the directory with all the GCMS you want to run: ")

finished_gcms = [name for name in os.listdir(gcm_folder) if os.path.isdir(os.path.join(gcm_folder, name))]
reassign = phases.copy()

#please do not touch these, the code will automatically determine which steps are necessary
STEP_ONE = False
STEP_TWO = False
STEP_THREE = True

os.chdir('PLANET_MODELS')
for j in range (len(finished_gcms)):
    nonrotatedcount = len(glob.glob(finished_gcms[j] + '*'))
    if nonrotatedcount !=4:
        STEP_ONE = True

os.chdir('..')
os.chdir('Spectra/DATA')

for x in range (len(finished_gcms)):
    initfilecount = len(glob.glob('init_' + str(finished_gcms[x]) + '*'))
    if initfilecount != len(phases):
        STEP_TWO = True

os.chdir('..')
os.chdir('..')

for i in range(len(finished_gcms)):
    if STEP_ONE:
        phases = [0.0]
        step = 'stepone'
        copyfiles(phases, step, finished_gcms[i])

        for j, phase in enumerate(phases):    
            for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + '_' + step + ".py"], inplace=True):
                if line.strip().startswith('STEP_ONE'):
                    line = '    STEP_ONE = True\n'
                sys.stdout.write(line)
            for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + '_' + step + ".py"], inplace=True):
                if line.strip().startswith('STEP_TWO'):
                    line = '    STEP_TWO = False\n'
                sys.stdout.write(line)
            for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + '_' + step + ".py"], inplace=True):
                if line.strip().startswith('STEP_THREE'):
                    line = '    STEP_THREE = False\n'
                sys.stdout.write(line)
            for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + '_' + step + ".py"], inplace=True):
                if line.strip().startswith('LOWEST_PHASE'):
                    line = 'LOWEST_PHASE = True\n'
                sys.stdout.write(line)
                
        step1jobnum = runsbatch(phases, source_file_name, finished_gcms[i], step)

    if STEP_TWO:
        phases = reassign
        step = 'steptwo'
        copyfiles(phases, step, finished_gcms[i])
        
        if STEP_ONE:
            dependency = step1jobnum
        else:
            dependency = 'none'
            
        for j, phase in enumerate(phases):
            for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + '_' + step + ".py"], inplace=True):
                if line.strip().startswith('STEP_ONE'):
                    line = '    STEP_ONE = False\n'
                sys.stdout.write(line)
            for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + '_' + step + ".py"], inplace=True):
                if line.strip().startswith('STEP_TWO'):
                    line = '    STEP_TWO = True\n'
                sys.stdout.write(line)
            for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + '_' + step + ".py"], inplace=True):
                if line.strip().startswith('STEP_THREE'):
                    line = '    STEP_THREE = False\n'
                sys.stdout.write(line)
            if phase == lowest_phase:
                for line in fileinput.input(
                        ["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + '_' + step + ".py"],
                        inplace=True):
                    if line.strip().startswith('LOWEST_PHASE'):
                        line = 'LOWEST_PHASE = True\n'
                    sys.stdout.write(line)
                    
        step2jobnum = runsbatch(phases, source_file_name, finished_gcms[i], step,dependency)

    if STEP_THREE:
        phases = reassign
        step = 'stepthree'
        copyfiles(phases, step, finished_gcms[i])
        
        if STEP_ONE:
            dependency = step1jobnum
        if STEP_TWO:
            dependency = step2jobnum
        else:
            dependency = 'none'
            
        for j, phase in enumerate(phases):
            for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + '_' + step + ".py"], inplace=True):
                if line.strip().startswith('STEP_ONE'):
                    line = '    STEP_ONE = False\n'
                sys.stdout.write(line)
            for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + '_' + step + ".py"], inplace=True):
                if line.strip().startswith('STEP_TWO'):
                    line = '    STEP_TWO = False\n'
                sys.stdout.write(line)
            for line in fileinput.input(["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + '_' + step + ".py"], inplace=True):
                if line.strip().startswith('STEP_THREE'):
                    line = '    STEP_THREE = True\n'
                sys.stdout.write(line)
                
            if phase == lowest_phase:
                for line in fileinput.input(
                        ["Spectra/run_spectra_" + finished_gcms[i] + "_" + str(phase) + '_' + step + ".py"],
                        inplace=True):
                    if line.strip().startswith('LOWEST_PHASE'):
                        line = 'LOWEST_PHASE = True\n'
                    sys.stdout.write(line)
            
            phaseslist = [phase]
            step3jobnum = runsbatch(phaseslist, source_file_name, finished_gcms[i], step, dependency)
            
            executable = f"Spectra/rt_emission_aerosols_{finished_gcms[i]}_phase_{str(phase)}.exe"
            
            while not os.path.exists(executable):
                print('Waiting for previous phase to create the .exe file')
                time.sleep(30)

