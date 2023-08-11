import glob
import os

files = glob.glob('slurm*')
#print(files)
y = False
for x in files:
    os.remove(x)
    print('removed', x)
    y = True

if y == True:
    print('removed all slurm files')
else:
    print('no slurm files found')
