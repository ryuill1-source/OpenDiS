import subprocess

def run_abaqus(config):
    ''' Input : config (dictionary) 
                with the following keys: jobname_head, ABAQUS_input_filename, umatname, num_cpus
    '''
    # Set local variables
    jobname_head = config.get("jobname_head", None)
    ABAQUS_input_filename = config.get("ABAQUS_input_filename", None)
    umatname = config.get("umatname", None)
    num_cpus = config.get("num_cpus", None)
    foldername_ABAQUS = config.get("foldername_ABAQUS", None)
    
    ABAQUS_jobname = jobname_head + ABAQUS_input_filename
    
    # string to run abaqus
    cmd_str = f"abaqus job={ABAQUS_jobname} user={umatname} cpus={num_cpus} double=both"
    
    # run ABAQUS
    result = subprocess.run(cmd_str, shell=True, cwd=foldername_ABAQUS,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)    
    
    print(f'    + Abaqus Job-{ABAQUS_input_filename} Execution Command: {cmd_str}')