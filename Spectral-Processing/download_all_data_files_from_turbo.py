import os
import shutil
import subprocess

def file_exists_on_remote(remote_user, remote_host, remote_file_path, control_path):
    """Check if a file exists on the remote server."""
    check_cmd = [
        "ssh", "-o", f"ControlPath={control_path}",
        f"{remote_user}@{remote_host}",
        f"test -e {remote_file_path} && echo True || echo False"
    ]
    result = subprocess.run(check_cmd, capture_output=True, text=True)
    return result.stdout.strip() == "True"

def download_file(remote_user, remote_host, remote_file_path, local_file_path, control_path):
    """Download file from remote server."""
    cmd = [
        "scp", "-o", f"ControlPath={control_path}",
        f"{remote_user}@{remote_host}:{remote_file_path}",
        local_file_path
    ]
    subprocess.run(cmd)

def unpack_downloaded_file(downloaded_file):
    """Unpack the downloaded file if it's zipped."""
    unpacked_name = downloaded_file.replace(".zip", "")

    if os.path.exists(unpacked_name):
        response = input(f"The directory {unpacked_name} already exists. If you are sure you want to overwrite it, type 'yes, really': ")
        if response == 'yes, really':
            shutil.rmtree(unpacked_name)
            print(f"Unpacking {downloaded_file}...")
            shutil.unpack_archive(downloaded_file)
        else:
            print(f"Skipping unpacking of {downloaded_file} as the directory exists.")
            return
    else:
        print(f"Unpacking {downloaded_file}...")
        shutil.unpack_archive(downloaded_file)

def move_files(source_dir, target_dir):
    """Move files from source directory to target directory."""
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)

    for file_name in os.listdir(source_dir):
        source = os.path.join(source_dir, file_name)
        target = os.path.join(target_dir, file_name)

        if os.path.exists(target):
            response = input(f"The file {target} already exists. Do you want to overwrite it? (yes or no): ")
            if response.lower() == 'yes':
                shutil.move(source, target)
                print(f"Moved file: {file_name}")
            else:
                print(f"Skipped moving file: {file_name}")
        else:
            shutil.move(source, target)
            print(f"Moved file: {file_name}")



def prompt_to_delete_files(file_paths):
    """Prompt user to delete specified files."""
    response = input("Do you want to delete the downloaded zipped files? (yes or no): ")
    if response.lower() == 'yes':
        for file_path in file_paths:
            if os.path.exists(file_path):
                os.remove(file_path)
                print(f"Deleted: {file_path}")
            else:
                print(f"File not found: {file_path}")
    else:
        print("Skipped deleting files.")

def download_and_process_files(remote_user, remote_host, files):
    """Download and unpack multiple files if they're zipped."""
    control_path_dir = "/tmp/{}@{}:{}"
    control_path = control_path_dir.format(remote_user, remote_host, 22)

    # Start the master connection
    master_cmd = [
        "ssh", "-Nf", "-o", f"ControlPath={control_path}",
        "-o", "ControlMaster=yes",
        f"{remote_user}@{remote_host}"
    ]
    subprocess.run(master_cmd)

    for remote_file in files:
        if not file_exists_on_remote(remote_user, remote_host, remote_file, control_path):
            print(f"Error: {remote_file} does not exist on the remote server.")
            continue

        downloaded_file = os.path.basename(remote_file)
        download_file(remote_user, remote_host, remote_file, downloaded_file, control_path)

        if downloaded_file.endswith(".zip"):
            unpack_downloaded_file(downloaded_file)

            unpacked_dir = os.path.splitext(downloaded_file)[0]

            # Move the files to respective folders
            if "SET_" in downloaded_file:
                move_files(unpacked_dir, "Spectra/DATA")
            elif "SCATTERING_DATA" in downloaded_file:
                move_files(unpacked_dir, "Spectra/SCATTERING_DATA")

    # Close the master connection
    close_cmd = ["ssh", "-O", "exit", "-o", f"ControlPath={control_path}", f"{remote_user}@{remote_host}"]
    subprocess.run(close_cmd)

    # Prompt to delete the downloaded zipped files
    prompt_to_delete_files([os.path.basename(f) for f in files if f.endswith(".zip")])

# List of files to download
files_to_download = [
    "/nfs/turbo/lsa-erausche/Post-Processing-Data-Files/SET_1.zip",
    "/nfs/turbo/lsa-erausche/Post-Processing-Data-Files/SET_2.zip",
    "/nfs/turbo/lsa-erausche/Post-Processing-Data-Files/SET_3.zip",
    "/nfs/turbo/lsa-erausche/Post-Processing-Data-Files/SCATTERING_DATA.zip"
]

# Prompt for username
remote_user = input("Please enter your username: ")

# Run the main function
download_and_process_files(
    remote_user=remote_user,
    remote_host="greatlakes.arc-ts.umich.edu",
    files=files_to_download
)
