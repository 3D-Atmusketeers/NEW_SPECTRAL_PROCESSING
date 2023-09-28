import os
import shutil
import fcntl

def replace_files(opacity_files, MET_X_SOLAR):
    # Create a lock file. This can be any file, but it's just used to handle the locking mechanism.
    lock_file_path = f'file_lock_{opacity_files}.lock'
    with open(lock_file_path, 'w') as lock_file:
        # Acquire the lock
        fcntl.flock(lock_file, fcntl.LOCK_EX)

        # Perform the copy operations
        files_to_replace = [
            "opac.h", "input.h", "readchemtable.c", 
            "readopactable.c", "template_inputs.h", "totalopac.c"
        ]
        for file_name in files_to_replace:
            dst_path = os.path.join('.', file_name)
            if os.path.exists(dst_path):
                os.remove(dst_path)
                
            print("replacing", file_name, "in", dst_path)
            src_path = os.path.join('OPAC_CODE_VERSIONS', opacity_files, file_name)
            
            if not os.path.exists(src_path):
                print(f"Error: {src_path} does not exist!")
                continue

            shutil.copy(src_path, dst_path)

        # Release the lock
        fcntl.flock(lock_file, fcntl.LOCK_UN)

    return None