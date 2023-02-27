import re
import pandas as pd
import numpy as np

def get_input_data(path, runname, input_file, input_param):
    # define the input_param and the regex pattern
    pattern = r"\b" + input_param + r"\s+(.*)"

    # compile the regex pattern
    regex = re.compile(pattern)

    # open the file and read its contents
    with open(path + runname + "/" + input_file, "r") as f:
        text = f.readlines()

    # find all the matches in the text
    for line in text:
        # skip lines that contain an exclamation point
        if "!" in line:
            continue

        # find matches on the current line
        matches = regex.findall(line)

        # print the matches
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

                # if there is only one number, print it
                return values[0]
            else:
                values = [float(i) for i in values]

                # if there are multiple values, print the list
                return values



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
    d = np.zeros((m+1, n+1))

    # initialize the first row and column
    for i in range(m+1):
        d[i, 0] = i
    for j in range(n+1):
        d[0, j] = j

    # fill in the matrix
    for j in range(1, n+1):
        for i in range(1, m+1):
            if s[i-1] == t[j-1]:
                cost = 0
            else:
                cost = 1
            d[i, j] = min(d[i-1, j]+1, d[i, j-1]+1, d[i-1, j-1]+cost)

    # return the distance value
    return d[m, n]