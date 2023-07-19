import glob
import os
#removes the successfully ran python, bash, and exe scripts that were run through
#the run_entire_suite.py

#you can run this by typing 'python3 Clean_suite.py' in the command line

bash = glob.glob('Run_sbatch_*')

y = False
for x in bash:
    os.remove(x)
    print('removed', x)
    y = True
if y == True:
    print('removed all redundant bash files')
else:
    print('no redundant bash files found')

pythonfiles = glob.glob('run_spectra_*')

z = False

for a in pythonfiles:
    os.remove(a)
    print('removed', a)
    z = True
if z == True:
    print('removed all redundant python files')
else:
    print('no redundant python files found')

exefiles = glob.glob('*.exe')

b = False

for c in exefiles:
    os.remove(c)
    print('removed', c)
    b = True
if b == True:
    print('removed all exe files')
else:
    print('no exe files found')

#below is the function that automatically deletes the files as they are made 
def automaticclean(x):
    namestring = x

    if namestring == 'run_spectra.py':
        print('its a good thing I added this toggle, or you wouldve removed run_spectra.py and Run_sbatch')
    else:
        batchmid = namestring.replace('run_spectra', 'Run_sbatch')
        batch = batchmid.replace('.py','')
        os.remove(batch)

        if STEP_THREE:
            pos = -1
            #print('namestring = ', namestring)
            #print('to string = ', str(namestring))
            while pos < 0:
                if namestring[pos] == '_':
                    namedstring = str(namestring)
                    named = namedstring[:pos]
                    pos = 1
                else:
                    pos -= 1
            exename = namestring.replace(named, named + '_phase')
            exemid = exename.replace('run_spectra','rt_emission_aerosols')
            exe = exemid.replace('.py', '.exe')
            os.remove(exe)

        os.remove(__file__)
