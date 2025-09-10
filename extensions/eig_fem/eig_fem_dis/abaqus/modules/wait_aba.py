import os
import time

def wait_aba(foldername_ABAQUS, ABAQUS_jobname, num_cpus, flag_type = "stress_ready", flag_check_interval=1.0):
    """
    Wait until ABAQUS fortran calculation
    (08/14/2025, kyeongmi) [Completed] Modified to work with OpenDiS.eig_fem
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
        # flags_exist = any(os.path.isfile(f"{foldername_ABAQUS}/ABAQUSworking{j:04d}.flag") for j in range(1, num_cpus + 1))
        # if not flags_exist:
        #     break

        # time.sleep(flag_check_interval)
        
        # (08/14/2025, kyeongmi) check wheter ABAQUS_stress_ready.flag exists
        ABAQUS_stress_ready = f"{foldername_ABAQUS}/ABAQUS_stress_ready.flag"
        ABAQUS_pause = f"{foldername_ABAQUS}/ABAQUS_pause.flag"
        
        if flag_type == "stress_ready":
            if os.path.exists(ABAQUS_stress_ready):
                stress_ready = True
            else:
                stress_ready = False
                print(f"wait_aba : {ABAQUS_stress_ready} does not exist")
            
            if stress_ready:
                print(f"wait_aba : {ABAQUS_stress_ready} exists, break the loop")
                break
        
        if flag_type == "pause":
            if os.path.exists(ABAQUS_pause):
                pause = True
            else:
                pause = False
                print(f"wait_aba : {ABAQUS_pause} does not exist")
            
            if pause:
                print(f"wait_aba : {ABAQUS_pause} exists, break the loop")
                break


        time.sleep(flag_check_interval)