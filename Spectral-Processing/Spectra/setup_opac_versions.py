import os
import fcntl
from shutil import copyfile
import shutil
import sys

def replace_files(opacity_files):
    lock_file_path = f'file_lock_{opacity_files}.lock'
    with open(lock_file_path, 'w') as lock_file:
        # Acquire the lock
        fcntl.flock(lock_file, fcntl.LOCK_EX)

        # Define the source directory path
        src_directory_path = os.path.join('OPAC_CODE_VERSIONS', opacity_files)

        # List and print files in the source directory for debugging
        src_files = os.listdir(src_directory_path)

        # Define template files expected to be replaced
        template_files_to_replace = [
            "template_opac.h", "template_input.h", "template_readchemtable.c",
            "template_readopactable.c", "template_totalopac.c"
        ]

        # Loop through each template file to replace
        for template_file in template_files_to_replace:
            # Construct source path
            src_path = os.path.join(src_directory_path, template_file)

            # Remove 'template_' prefix for destination file name
            dst_file_name = template_file.replace('template_', '')
            dst_path = os.path.join('.', dst_file_name)

            # Perform the copy if source file exists
            try:
                if os.path.exists(dst_path):
                    os.remove(dst_path)  # Remove existing file if present

                if not os.path.exists(src_path):
                    print(f"Error: Source file {src_path} does not exist!")
                    continue

                #print(f"Copying {src_path} to {dst_path}")
                shutil.copy(src_path, dst_path)
            except Exception as e:
                print(f"Unable to copy file from {src_path} to {dst_path}. Error: {e}")

        # Release the lock
        fcntl.flock(lock_file, fcntl.LOCK_UN)

def modify_input_h(opacity_files, modifications):
    inputs_file = 'input.h'
    lock_file_path = 'file_lock.lock'

    # Ensure atomic modifications by acquiring a file lock
    with open(lock_file_path, 'a') as lock_file:  # 'a' mode to append or create if it doesn't exist
        fcntl.flock(lock_file, fcntl.LOCK_EX)  # Acquire an exclusive lock

        # Copy the template for inputs as the basis for modifications
        try:
            template_input_path = os.path.join('OPAC_CODE_VERSIONS', opacity_files, 'template_input.h')
            copyfile(template_input_path, inputs_file)
            #print(f"Copied template_input.h to {inputs_file}")
        except IOError as e:
            print(f"Unable to copy file: {e}")
            sys.exit(1)
        except Exception as e:
            print(f"Unexpected error: {e}", sys.exc_info())
            sys.exit(1)

        # Read in the file
        try:
            with open(inputs_file, 'r') as file:
                filedata = file.read()

            # Perform replacements based on modifications dictionary
            for placeholder, value in modifications.items():
                filedata = filedata.replace(placeholder, value)

            # Write the modified content back to the input file
            with open(inputs_file, 'w') as file:
                file.write(filedata)
            #print("input.h has been successfully modified.")
        except IOError as e:
            print(f"Error while modifying {inputs_file}: {e}")
            sys.exit(1)
        except Exception as e:
            print(f"Unexpected error while modifying {inputs_file}: {e}", sys.exc_info())
            sys.exit(1)

        finally:
            # Ensure the lock is always released
            fcntl.flock(lock_file, fcntl.LOCK_UN)

def insert_opacity_definitions(filepath, directory, opacity_species):
    lock_file_path = 'global_file_lock.lock'  # Use a global lock file for shared resource management

    # Open or create the lock file and acquire an exclusive lock
    with open(lock_file_path, 'w') as lock_file:
        fcntl.flock(lock_file, fcntl.LOCK_EX)

        try:
            # Open the file for reading and writing
            with open(filepath, 'r+') as file:
                filedata = file.read()

                endif_position = filedata.rfind("#endif")
                if endif_position == -1:
                    print("Could not find #endif marker in the file.")
                    return

                # Iterate over the provided list of opacity species
                for species in opacity_species:
                    # Construct the line to insert for each species
                    new_line = f'#define {species}_FILE   "{directory}/opac{species}.dat"\n'
                    # Insert the new line before the #endif marker
                    filedata = filedata[:endif_position] + new_line + filedata[endif_position:]
                    # Update the position of #endif to include the newly added line
                    endif_position += len(new_line)

                # Move the pointer to the start of the file and write the modified content
                file.seek(0)
                file.write(filedata)
                # Truncate the file in case the new content is shorter than the old
                file.truncate()

            print(f"File '{filepath}' updated successfully.")

        finally:
            # Release the lock
            fcntl.flock(lock_file, fcntl.LOCK_UN)


import fcntl

def modify_totalopac(opacity_files, species_to_include):
    totalopac_path = 'totalopac.c'
    lock_file_path = 'file_lock_totalopac.lock'

    # Attempt to open or create the lock file, then acquire an exclusive lock
    with open(lock_file_path, 'w') as lock_file:
        fcntl.flock(lock_file, fcntl.LOCK_EX)

        try:
            # Read in the totalopac.c file
            with open(totalopac_path, 'r') as file:
                lines = file.readlines()

            # Find the index where new lines should be added for opac structures
            index_to_add_opac_structs = lines.index("// Add in the new lines for each opacity species\n")

            # Add new lines for each opacity species
            for species in species_to_include:
                species_lines = [
                    f"struct Opac opac{species};\n"
                ]
                lines.insert(index_to_add_opac_structs + 1, species_lines[0])

            # Find the index where new lines should be added for filling opacities
            index_to_add_fill_opacities = lines.index("  // Fill in the opacities for each species\n")

            # Add new lines for filling opacities for each species
            for species in species_to_include:
                species_lines = [
                    f"  /* Fill in {species} opacities */\n",
                    f"  opac{species}.T = dvector(0, NTEMP-1);\n",
                    f"  opac{species}.P = dvector(0, NPRESSURE-1);\n",
                    f"  opac{species}.Plog10 = dvector(0, NPRESSURE-1);\n",
                    f"  opac{species}.kappa = malloc(NLAMBDA*sizeof(double));\n",
                    f"  for(i=0; i<NLAMBDA; i++){{\n",
                    f"    opac{species}.kappa[i] = malloc(NPRESSURE*sizeof(double));\n",
                    f"    for(j=0; j<NPRESSURE; j++){{\n",
                    f"      opac{species}.kappa[i][j] = malloc(NTEMP*sizeof(double));\n",
                    f"    }}\n",
                    f"  }}\n",
                    f"  opac{species}.abundance = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);\n",
                    f"  for(j=0; j<NPRESSURE; j++){{\n",
                    f"    for(k=0; k<NTEMP; k++){{\n",
                    f"      opac{species}.abundance[j][k] = chem.{species}[j][k];\n",
                    f"    }}\n",
                    f"  }}\n",
                    f"  strcpy(filename, {species}_FILE);\n",
                    f"  ReadOpacTable(opac{species}, filename);\n",
                    f"  printf(\"Read {species} Opacity done\\n\");\n\n"
                ]
                lines.insert(index_to_add_fill_opacities + 1, "".join(species_lines))

            # Find the index where new lines should be added for final summation
            index_to_add_summation = lines.index("  // Do a final summation over the species\n")

            # Add new lines for final summation over species
            summation_lines = [
                "  // Do a final summation over the species\n",
                "  for (i=0; i<NLAMBDA; i++) {\n",
                "    for (j=0; j<NPRESSURE; j++) {\n",
                "      for (k=0; k<NTEMP; k++) {\n",
                "          opac.kappa[i][j][k] = "
            ]
            for species in species_to_include:
                summation_lines.append(f"                              + opac{species}.kappa[i][j][k]\n")
            summation_lines.append("                              ;\n")
            summation_lines.extend([
                "      }\n",
                "    }\n",
                "  }\n"
            ])
            lines.insert(index_to_add_summation + 1, "".join(summation_lines))

            # Find the index where new lines should be added for freeing uneeded opacity structures
            index_to_add_free_opacities = lines.index("  // Free not needed opacity structures and chemistry table\n")

            # Add new lines for freeing uneeded opacity structures
            for species in species_to_include:
                free_lines = [
                    f"  FreeOpacTable(opac{species});\n"
                ]
                lines.insert(index_to_add_free_opacities + 1, free_lines[0])

            # Write the modified content back
            with open(totalopac_path, 'w') as file:
                file.writelines(lines)

            print(f"'{totalopac_path}' has been successfully modified.")

        except Exception as e:
            print(f"Error modifying '{totalopac_path}': {e}")

        finally:
            # Ensure the lock is always released
            fcntl.flock(lock_file, fcntl.LOCK_UN)
