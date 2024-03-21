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

    # Close the master connection
    close_cmd = ["ssh", "-O", "exit", "-o", f"ControlPath={control_path}", f"{remote_user}@{remote_host}"]
    subprocess.run(close_cmd)

# List of files to download
files_to_download = [
#    "/nfs/turbo/lsa-erausche/Post-Processing-Data-Files/SET_1.zip",
#    "/nfs/turbo/lsa-erausche/Post-Processing-Data-Files/SCATTERING_DATA.zip"
]

# Prompt for username
remote_user = input("Please enter your username: ")

# Run the main function
download_and_process_files(
    remote_user=remote_user,
    remote_host="greatlakes.arc-ts.umich.edu",
    files=files_to_download
)
