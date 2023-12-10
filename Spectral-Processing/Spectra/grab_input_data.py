#below it imports the compiled fortran subroutine as a python module
#import fortrantopythonfile as fp
import re
import pandas as pd
import numpy as np


def create_dict():
    # Call the fortran

    print('here!')
    fp.readinvals()
    print('done!')

    # Create dictionary
    fort7dict = {}

    #assigns the dictionary values depending on the location from the compiled fortran
    fort7dict['THECOMMENT'] = str(fp.commen.thecomment.tolist())[2:-2].strip()

    fort7dict['GA'] = fp.ppl.ga.tolist()
    fort7dict['GASCON'] = fp.ppl.gascon.tolist()
    fort7dict['RADEA'] = fp.ppl.radea.tolist()
    fort7dict['AKAP'] = fp.ppl.akap.tolist()
    fort7dict['WW'] = fp.ppl.ww.tolist()
    fort7dict['P0'] = fp.ppl.p0.tolist()
    fort7dict['RV'] = fp.ppl.rv.tolist()
    fort7dict['CLATNT'] = fp.ppl.clatnt.tolist()

    fort7dict['OOM_IN'] = fp.varparam.oom_in.tolist()
    fort7dict['LPLOTMAP'] = bool(fp.varparam.lplotmap.tolist())
    fort7dict['NLPLOTMAP_IN'] = fp.varparam.nlplotmap_in.tolist()
    fort7dict['RFCOEFF_IN'] = fp.varparam.rfcoeff_in.tolist()
    fort7dict['NTSTEP_IN'] = fp.varparam.ntstep_in.tolist()
    fort7dict['BOTRELAXTIME'] = fp.varparam.botrelaxtime.tolist()
    fort7dict['FBASEFLUX'] = fp.varparam.fbaseflux.tolist()
    fort7dict['FORCE1DDAYS'] = fp.varparam.force1ddays.tolist()
    fort7dict['OPACIR_POWERLAW'] = fp.varparam.opacir_powerlaw.tolist()
    fort7dict['OPACIR_REFPRES'] = fp.varparam.opacir_refpres.tolist()
    fort7dict['SOLC_IN'] = fp.varparam.solc_in.tolist()
    fort7dict['TOAALB'] = fp.varparam.toaalb.tolist()
    fort7dict['PORB'] = fp.varparam.porb.tolist()
    fort7dict['OBLIQ'] = fp.varparam.obliq.tolist()
    fort7dict['ECCEN'] = fp.varparam.eccen.tolist()

    fort7dict['LBDRAG'] = bool(fp.mag.lbdrag.tolist())
    fort7dict['BFIELD'] = fp.mag.bfield.tolist()
    fort7dict['TDRAG_MIN'] = fp.mag.tdrag_min.tolist()
    fort7dict['RAMPUP'] = fp.mag.rampup.tolist()

    fort7dict['LBIN'] = bool(fp.binval.lbin.tolist())
    fort7dict['PORBST'] = fp.binval.porbst.tolist()
    fort7dict['ECCPL'] = fp.binval.eccpl.tolist()
    fort7dict['ECCST'] = fp.binval.eccst.tolist()
    fort7dict['SMAPL'] = fp.binval.smapl.tolist()
    fort7dict['SMAST'] = fp.binval.smast.tolist()
    fort7dict['STMASS1'] = fp.binval.stmass1.tolist()
    fort7dict['STMASS2'] = fp.binval.stmass2.tolist()
    fort7dict['STRAD1'] = fp.binval.strad1.tolist()
    fort7dict['STRAD2'] = fp.binval.strad2.tolist()
    fort7dict['STTEMP1'] = fp.binval.sttemp1.tolist()
    fort7dict['STTEMP2'] = fp.binval.sttemp2.tolist()

    fort7dict['LLOGPLEV'] = bool(fp.simprad.llogplev.tolist())
    fort7dict['LFLUXDIAG'] = bool(fp.simprad.lfluxdiag.tolist())
    fort7dict['L1DZENITH'] = bool(fp.simprad.l1dzenith.tolist())
    fort7dict['LDIUR'] = bool(fp.simprad.ldiur.tolist())
    fort7dict['JSKIPLON'] = int(fp.simprad.jskiplon.tolist())
    fort7dict['JSKIPLAT'] = int(fp.simprad.jskiplat.tolist())
    fort7dict['DOSWRAD'] = bool(fp.simprad.doswrad.tolist())
    fort7dict['DOLWRAD'] = bool(fp.simprad.dolwrad.tolist())
    fort7dict['LWSCAT'] = bool(fp.simprad.lwscat.tolist())
    fort7dict['FLXLIMDIF'] = bool(fp.simprad.flxlimdif.tolist())
    fort7dict['SURFEMIS'] = fp.simprad.surfemis.tolist()
    fort7dict['RAYSCAT'] = bool(fp.simprad.rayscat.tolist())
    fort7dict['RAYSCATLAM'] = fp.simprad.rayscatlam.tolist()
    fort7dict['AEROSOLS'] = bool(fp.simprad.aerosols.tolist())
    fort7dict['ABSSW'] = fp.simprad.abssw.tolist()
    fort7dict['ABSLW'] = fp.simprad.abslw.tolist()
    fort7dict['ALBSW'] = fp.simprad.albsw.tolist()
    fort7dict['NEWTB'] = fp.simprad.newtb.tolist()
    fort7dict['NEWTE'] = fp.simprad.newte.tolist()
    fort7dict['with_TiO_and_VO'] = fp.simprad.with_tio_and_vo.tolist()
    fort7dict['picket_fence_optical_depths'] = bool(fp.simprad.picket_fence_optical_depths.tolist())

    fort7dict['AEROSOLMODEL'] = str(fp.cloudy.aerosolmodel.tolist())[2:-2].strip()
    fort7dict['AEROSOLCOMP'] = str(fp.cloudy.aerosolcomp.tolist())[2:-2].strip()
    fort7dict['HAZETYPE'] = str(fp.cloudy.hazetype.tolist())[2:-2].strip()
    fort7dict['MTLX'] = fp.cloudy.mtlx.tolist()
    fort7dict['METALLICITY'] = fp.cloudy.metallicity.tolist()
    fort7dict['HAZES'] = bool(fp.cloudy.hazes.tolist())
    fort7dict['PICKET_FENCE_CLOUDS'] = bool(fp.cloudy.picket_fence_clouds.tolist())
    fort7dict['MOLEF'] = fp.cloudy.molef.tolist()
    fort7dict['AERLAYERS'] = fp.cloudy.aerlayers.tolist()
    fort7dict['AERTOTTAU'] = fp.cloudy.aertottau.tolist()
    fort7dict['CLOUDBASE'] = fp.cloudy.cloudbase.tolist()
    fort7dict['CLOUDTOP'] = fp.cloudy.cloudtop.tolist()
    fort7dict['AERHFRAC'] = fp.cloudy.aerhfrac.tolist()
    fort7dict['PI0AERSW'] = fp.cloudy.pi0aersw.tolist()
    fort7dict['ASYMSW'] = fp.cloudy.asymsw.tolist()
    fort7dict['EXTFACTLW'] = fp.cloudy.extfactlw.tolist()
    fort7dict['PI0AERLW'] = fp.cloudy.pi0aerlw.tolist()
    fort7dict['ASYMLW'] = fp.cloudy.asymlw.tolist()
    fort7dict['DELTASCALE'] = bool(fp.cloudy.deltascale.tolist())
    fort7dict['SIG_AREA'] = fp.cloudy.sig_area.tolist()
    fort7dict['PHI_LON'] = fp.cloudy.phi_lon.tolist()
    return fort7dict

def get_input_data(path, runname, input_file, input_param):
    # define the input_param and the regex pattern for numbers and booleans
    pattern = r"\b" + input_param + r"\s*=\s*(.*?)\s*(?:$|\||&)"

    # define the regex pattern for strings
    string_pattern = r"\b" + input_param + r"\s*=\s*'([^']*)'"

    # compile the regex patterns
    regex = re.compile(pattern)
    string_regex = re.compile(string_pattern)

    # open the file and read its contents
    with open(path + runname + "/" + input_file, "r") as f:
        text = f.readlines()

    # find all the matches in the text
    for line in text:
        # skip lines that contain an exclamation point
        if "!" in line:
            continue

        # find matches for strings first
        string_matches = string_regex.findall(line)
        if string_matches:
            return string_matches[0]

        # find matches for numbers and booleans
        matches = regex.findall(line)

        for match in matches:
            # extract the values from the match, including scientific notation
            if (input_param == 'RADEA'):
                values = re.findall(r"[+\-]?[^A-Za-z]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)", match)
            else:
                values = re.findall(r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?", match)

            if len(values) == 0:
                # define the regex pattern
                pattern = r"\b(T|F|True|False)\b"

                # compile the regex pattern
                bool_regex = re.compile(pattern)

                # find all the matches in the string
                bool_match = bool_regex.findall(line)
                return bool_match

            elif len(values) == 1:
                values = [float(i) for i in values]

                # if there is only one number, return it
                return values[0]
            else:
                values = [float(i) for i in values]

                # if there are multiple values, return the list
                return values

    return None



def read_planet_and_star_params(planet_name, column_name_str):
    # read the Excel sheet into a DataFrame
    df = pd.read_excel("eplanetpars.xlsx")

    # Get the planet name from the input string
    planet_name_base = re.split(r"[_|-]", planet_name)[0]

    # Calculate Levenshtein distance between input string and each string in the column
    distances = df['Name'].apply(lambda x: levenshtein_distance(planet_name_base, x))

    # get index of row with smallest distance
    best_match_index = distances.idxmin()

    value = df.iloc[best_match_index][column_name_str]

    return value


# define a function to calculate the Levenshtein distance between two strings
def levenshtein_distance(s, t):
    # create a matrix to store the distance values
    m, n = len(s), len(t)
    d = np.zeros((m + 1, n + 1))

    # initialize the first row and column
    for i in range(m + 1):
        d[i, 0] = i
    for j in range(n + 1):
        d[0, j] = j

    # fill in the matrix
    for j in range(1, n + 1):
        for i in range(1, m + 1):
            if s[i - 1] == t[j - 1]:
                cost = 0
            else:
                cost = 1
            d[i, j] = min(d[i - 1, j] + 1, d[i, j - 1] + 1, d[i - 1, j - 1] + cost)

    # return the distance value
    return d[m, n]
