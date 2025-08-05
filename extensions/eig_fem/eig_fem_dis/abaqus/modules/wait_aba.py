import os
import time

def wait_aba(foldername_ABAQUS, ABAQUS_jobname, num_cpus, flag_check_interval=1.0):
    """
    Wait until ABAQUS fortran file changes ABAQUSworking to Matlabworking
    """


    log_path = os.path.join(foldername_ABAQUS, f"{ABAQUS_jobname}.log")
    last_checked_mtime = None


    while True:
        # 1. if log file exists, check error in abaqus
        if os.path.isfile(log_path):
            mtime = os.path.getmtime(log_path)

            if last_checked_mtime is None or mtime != last_checked_mtime:
                last_checked_mtime = mtime

                with open(log_path, "r") as f:
                    contents = f.read()
                    if "Error" in contents:
                        if "end-of-file" in contents:
                            print("  +++ ABAQUS didn't work. end-of-file detected. File may be missing or not copied properly.")
                            raise RuntimeError()
                        else:
                            print("  +++ ABAQUS didn't work")
                            raise RuntimeError()


        # 2. check ABAQUSworking files
        flags_exist = any(os.path.isfile(f"{foldername_ABAQUS}/ABAQUSworking{j:04d}.flag") for j in range(1, num_cpus + 1))
        if not flags_exist:
            break

        time.sleep(flag_check_interval)