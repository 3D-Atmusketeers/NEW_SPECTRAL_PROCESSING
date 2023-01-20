import re
import pandas as pd


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

    # define the input string
    input_string = planet_name

    # get the row names from the DataFrame
    row_names = df.Name.values

    # Get the planet name from the input string
    planet_name_base = re.split(r"[_|-]", planet_name)[0]

    # check if any of the row names are a substring of the input string
    found = False
    for row_name in row_names:
        if row_name in planet_name_base or planet_name_base in row_name:
            found = True
            break

    if not found:
        print("The planet name isn't in the dataframe")
    else:
        # return the value for the row with the planet name
        # and the column name set as the string calling name
        value = list(df.loc[df['Name'] == row_name][column_name_str])[0]

        return value



