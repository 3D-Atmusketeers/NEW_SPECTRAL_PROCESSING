import os
os.system('-m numpy.f2py -c testing1.f90 -m fortranpythonfile')

import fortranpythonfile


testarrays = fortranpythonfile.testing1(1)
print(testarrays)
