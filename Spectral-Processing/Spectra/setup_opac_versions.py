import re
import os
import shutil
import time

def replace_files(opacity_files, MET_X_SOLAR):
    if os.path.exists("opac.h"):
        os.remove("opac.h")
    shutil.copy('OPAC_CODE_VERSIONS/' + opacity_files + '/opac.h', "opac.h")

    if os.path.exists("input.h"):
        os.remove("input.h")
    shutil.copy('OPAC_CODE_VERSIONS/' + opacity_files + '/input.h', "input.h")

    if os.path.exists("readchemtable.c"):
        os.remove("readchemtable.c")
    shutil.copy('OPAC_CODE_VERSIONS/' + opacity_files + '/readchemtable.c', "readchemtable.c")

    if os.path.exists("readopactable.c"):
        os.remove("readopactable.c")
    shutil.copy('OPAC_CODE_VERSIONS/' + opacity_files + '/readopactable.c', "readopactable.c")

    if os.path.exists("template_inputs.h"):
        os.remove("template_inputs.h")
    shutil.copy('OPAC_CODE_VERSIONS/' + opacity_files + '/template_inputs.h', "template_inputs.h")

    if os.path.exists("totalopac.c"):
        os.remove("totalopac.c")
    shutil.copy('OPAC_CODE_VERSIONS/' + opacity_files + '/totalopac.c', "totalopac.c")

    return None
