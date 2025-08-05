# //------------------------------------------------------------------------------------+
# //         COPYRIGHT (c) 2019, Taejoon Park and Ill Ryu.                              |
# //                          All rights reserved.                                      |
# //                                                                                    |
# // This work was produced at UT Dallas and The Ohio State University                  |
# // through collaborative research project on multiscale model for plasticity          |
# // 20190613                                                                           |
# //------------------------------------------------------------------------------------+
# Developement of Python version for DDM code
# colaborators: Gyu-Jang Sim, Kyeongmi Yeon, Jun-Ho Seong, Seong-Jin Yoon

import shutil
import subprocess
import time
import os
import glob

from modules.import_elements_nodes import import_elements_nodes
from modules.generate_input_ABAQUS_timeSync import generate_input_ABAQUS_timeSync
from modules.check_plastic_strain import check_plastic_strain
from modules.generate_input_ABAQUS_stressControl import generate_input_ABAQUS_stressControl
from modules.generate_input_ABAQUS_erateControl import generate_input_ABAQUS_erateControl
from modules.generate_segment_files_v1 import generate_segment_files_v1
from modules.wait_aba import wait_aba
from modules.generate_postprocessing_data12 import generate_postprocessing_data12
from modules.generate_postprocessing_data14 import generate_postprocessing_data14

def DDM_Linux_v0(config):

    # confige from input file (ex. D1000_FRsource_qstatic.m)
    # config from param.m

    ### Unpack variables TODO: remove unnecessary, and add necessary variables
    ## inside input m file (ex. D1000_FRsource_qstatic.m)
    burger = config.get("burger", None)
    L_elem = config.get("L_elem", None)
    density = config.get("density", None)
    E = config.get("E", None)
    rmax_DD = config.get("rmax_DD", None)
    sync_dd_fem_dt = config.get("sync_dd_fem_dt", None)
    LoadType = config.get("LoadType", None)
    erate = config.get("erate", None)
    InelasticStep = config.get("InelasticStep",1) # TODO: Do something...
    num_cpus = config.get("num_cpus", None)
    num_max_cycle = config.get("num_max_cycle", None)
    ABAQUS_input_filename = config.get("ABAQUS_input_filename", None)
    DD_Coupling_All_Element = config.get("DD_Coupling_All_Element", None)
    enabled_temp_disp = config.get("enabled_temp_disp", None)
    plot_freq = config.get("plot_freq", None)
    critical_eps_rate = config.get("critical_eps_rate", None)
    FEM_Volume = config.get("FEM_Volume", None)
    dd_size = config.get("dd_size", None)
    DD_Volume = config.get("DD_Volume", None)
    Vol_dd_fem_ratio = config.get("Vol_dd_fem_ratio", None)
    def_ParaDiS_DEBUG = config.get("def_ParaDiS_DEBUG", None)
    def_ABAQUS_DEBUG = config.get("def_ABAQUS_DEBUG", None)
    # systemtype_Linux = config.get("systemtype_Linux", None)
    # foldername_root = config.get("foldername_root", None)
    foldername_DDM = config.get("foldername_DDM", None)
    # foldername_wsl = config.get("foldername_wsl", None)
    foldername_temp_win = config.get("foldername_temp_win", None)
    # foldername_temp_wsl = config.get("foldername_temp_wsl", None)
    foldername_tests = config.get("foldername_tests", None)
    foldername_Job = config.get("foldername_Job", None)
    foldername_segs = config.get("foldername_segs", None)
    foldername_arm = config.get("foldername_arm", None)
    foldername_restart = config.get("foldername_restart", None)
    foldername_properties = config.get("foldername_properties", None)
    foldername_FEMstress = config.get("foldername_FEMstress", None)
    Element_type = config.get("Element_type", None)
    cmd_env = config.get("cmd_env", None)
    FEM_Stress_Scale = config.get("FEM_Stress_Scale", None)
    xyzlimit = config.get("xyzlimit", None)
    ParaDiS_Input = config.get("ParaDiS_Input", None)
    ParaDiS_Print = config.get("ParaDiS_Print", None)
    Jobname = config.get("Jobname", None)
    Jobname_fullpath = config.get("Jobname_fullpath", None)
    foldername_tests = config.get("foldername_tests", None)

    ## inside param.m
    # Time step control
    dt_FEM_each_step = config.get("dt_FEM_each_step", None)
    DD_maxstep_each_cycle = config.get("DD_maxstep_each_cycle", None)
    # ABAQUS setup
    jobname_head = config.get("jobname_head", None)
    envname = config.get("envname", None)
    umatname = config.get("umatname", None)
    cmd_env = config.get("cmd_env", None)
    foldername_ABAQUS = config.get("foldername_ABAQUS", None)
    filename_merge_seg_shell = config.get("filename_merge_seg_shell", None)
    # Effective slip zone
    EffSlip = config.get("EffSlip", None)
    # Surface element/node sets
    surface_element_set_name = config.get("surface_element_set_name", None)
    surface_node_set_name = config.get("surface_node_set_name", None)
    formatSpec_Elset = config.get("formatSpec_Elset", None)
    Elset_col_num = config.get("Elset_col_num", None)
    # Dislocation segment info
    foldername_segs = config.get("foldername_segs", None)
    foldername_arm = config.get("foldername_arm", None)
    foldername_restart = config.get("foldername_restart", None)
    foldername_properties = config.get("foldername_properties", None)
    foldername_cmdout = config.get("foldername_cmdout", None)
    filename_seg_head = config.get("filename_seg_head", None)
    filename_seg_num = config.get("filename_seg_num", None)
    filename_seg_temp = config.get("filename_seg_temp", None)
    filename_seg_lastCycle = config.get("filename_seg_lastCycle", None)
    # ABAQUS FEM output
    filename_FEMoutput_stress = config.get("filename_FEMoutput_stress", None)
    filename_FEMoutput_nodes = config.get("filename_FEMoutput_nodes", None)
    filename_link_element = config.get("filename_link_element", None)
    filename_center_element = config.get("filename_center_element", None)
    filename_surface_stress = config.get("filename_surface_stress", None)
    filename_rn_check = config.get("filename_rn_check", None)
    foldername_FEMstress = config.get("foldername_FEMstress", None)
    filename_FEMstress_head = config.get("filename_FEMstress_head", None)
    filename_NUCsite_head = config.get("filename_NUCsite_head", None)
    filename_surface_int = config.get("filename_surface_int", None)
    # ParaDiS executables
    paradis_bulk_exe = config.get("paradis_bulk_exe", None)
    paradis_nuc_exe = config.get("paradis_nuc_exe", None)








    # Initialize time and loading variables
    N_DD = 1 # For the first DDM cycle, just run 1 time step in DD (for both with and without DD-FEM synchronization)
    Total_N_DD = 0 # The total number of cycles run in DD

    paradis_timeNow = 0. # Current time in ParaDiS
    Dt_paradis_perCycle = 0. # the time step (delta T) in ParaDiS within this cycle
    FEM_Du = 0.
    FEM_Dsigma = 0.
    current_disp = 0.
    current_stress = 0.
    Load_Stage = 1 # 1: load 2: hold 3: unload



    

    
    # Import Elements and Nodes, then calculate min, max, and centroid of elements
    Nodes, Elements, SurfaceElementSet, surface_element_flag, Elements_Min, Elements_Max, Elements_Center, Num_Nodes, Num_Elements, Num_SurfaceElement = \
        import_elements_nodes( 
            foldername_tests, foldername_Job, 
            foldername_ABAQUS,
            jobname_head, ABAQUS_input_filename,
            config.get("Element_type"), config.get("surface_element_set_name"),
            config.get("formatSpec_Elset"), config.get("Elset_col_num")
        )

    
    # Copy UMAT TODO: check final directory
    # shutil.copy(foldername_wsl + foldername_root + '/' + foldername_DDM + '/' + umatname, 
    #             foldername_ABAQUS)
    shutil.copy(os.path.join('..','..',foldername_DDM,umatname), foldername_ABAQUS)
    # shutil.copy(
    #     os.path.join(foldername_wsl, foldername_root, foldername_tests, foldername_ABAQUS, umatname),
    #     os.path.join(foldername_ABAQUS, umatname)
    # )

    
    with open(f"./{foldername_ABAQUS}/surface_elements.txt", "w") as file_elem_surface:
        for j in range(Num_Elements):
            file_elem_surface.write(f"{j + 1}  {surface_element_flag[j]}\n")


    original_length = config.get('original_length', 0)


    #
    ## ABAQUS Execution
    for i in range(1, num_cpus+1):
        fileID = f"{foldername_ABAQUS}/ABAQUSworking{i:04d}.flag"
        open(fileID, "w").close()
    
    print("  +++ Abaqus Job-" + ABAQUS_input_filename + " EXECUTION +++")
    
    ABAQUS_jobname = jobname_head + ABAQUS_input_filename



    # Terminate previous abaqus job, if exist
    cmd_str = f"abaqus terminate job={ABAQUS_jobname}"
    result = subprocess.run(cmd_str, shell=True, cwd=foldername_ABAQUS,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if "Sent Termination message to ABAQUS job" in result.stdout:
        print(f"    + Previous ABAQUS execution was terminated: {result.stdout.strip()} +")



    #!
    cmd_str = f"abaqus job={ABAQUS_jobname} user={umatname} cpus={num_cpus} double=both"
    result = subprocess.run(cmd_str, shell=True, cwd=foldername_ABAQUS,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)    
    print(f'    + Abaqus Job-{ABAQUS_input_filename} Execution Command: {cmd_str}')
    if result.returncode != 0:
        print(f"  +++ ABAQUS execution Failed: {result.stderr.strip()} +++")
        raise RuntimeError()

    # Wait until ABAQUS fortran file changes ABAQUSworking to Matlabworking
    wait_aba(foldername_ABAQUS, ABAQUS_jobname, num_cpus)

    # Now, continue Matlab script.
    ######################################################################################################################################################################################################################################
    ######################################################################################################################################################################################################################################
    ########################################################################################################## MAIN LOOP STARTS ##########################################################################################################
    ######################################################################################################################################################################################################################################
    ######################################################################################################################################################################################################################################

    for i in range(1, num_max_cycle+1):
        ## Display step number
        print(f"===== DDM STEP {i:04d} START =====")
        ## Execute ParaDiS to calculate the dislcation movements
        # if i == 1:
        #     paradisrunstr = f"{foldername_root}/bin/{paradis_bulk_exe} {foldername_root}/{foldername_tests}/{foldername_Job}/{Jobname}.ctrl"
        # else:
        #     paradisrunstr = f"{foldername_root}/bin/{paradis_nuc_exe} -d {foldername_root}/{foldername_tests}/{foldername_Job}/{foldername_restart}/restart.data {foldername_root}/{foldername_tests}/{foldername_Job}/{Jobname}.ctrl"
        # job_dir = os.path.join(foldername_root,foldername_tests, foldername_Job)
        job_dir = os.path.join(foldername_tests, foldername_Job)
        ctrl_file_path = os.path.join(job_dir, f'{Jobname}.ctrl')
        restart_data_path = os.path.join(job_dir, 'restart', 'restart.data')
        restart_cn_path = os.path.join(job_dir, 'restart', 'restart.cn')
        step_txt_path = os.path.join(job_dir, foldername_cmdout, f'step{i:04d}.txt')
        # exe_FEM_DD = os.path.join(foldername_root,'bin','paradis_FEM_DD')
        exe_FEM_DD = os.path.abspath(os.path.join('..','..','bin','paradis_FEM_DD'))
        # exe_FEM_DD_Nuc = os.path.join(foldername_root,'bin','paradis_FEM_DD_Force_Nuc')
        exe_FEM_DD_Nuc = os.path.abspath(os.path.join('..','..','bin','paradis_FEM_DD_Force_Nuc'))

        if i == 1:

            paradisrunstr = f'{exe_FEM_DD} {ctrl_file_path}'

        else:

            paradisrunstr = (
                f'{exe_FEM_DD_Nuc} -d {restart_data_path} {restart_cn_path} '
                f'|& tee {step_txt_path}'
            )

        # paradis execution
        print(f"+++ ParaDis {i:04d} Cycle EXECUTION +++")
        print(f"+++ ParaDis Execution Command: {paradisrunstr} +++")

        # breakpoint()
        # result = subprocess.run(paradisrunstr, shell=True, cwd=foldername_root,
        #                         stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        command = f"cd {os.path.abspath(os.path.join('..','..'))} && {paradisrunstr}" # This doesn't actually changes main direcotry of code.
        return_code = os.system(command)
          
        Total_N_DD=Total_N_DD+N_DD

        if result.returncode != 0:

            print(f"  +++ ParaDiS execution Failed: {result.stderr.strip()} +++")
            print(f" {result.stdout}")
            raise RuntimeError()
        
        print(f"+++ ParaDis {i:04d} Cycle Finished +++")
        
        ## Format change (ParaDis to ABAQUS)

        # ArmFolderFullName = os.path.join(foldername_wsl, foldername_root, foldername_tests, foldername_Job, foldername_arm)
        ArmFolderFullName = os.path.join('..', '..', foldername_tests, foldername_Job, foldername_arm)
        # SegFolderFullName = os.path.join(foldername_wsl, foldername_root, foldername_tests, foldername_Job, foldername_segs)
        SegFolderFullName = os.path.join('..', '..', foldername_tests, foldername_Job, foldername_segs)

        # Get the deltaT in paradis in this cycle
        if sync_dd_fem_dt: # Synchronize timestep in FEM and DD simulation
            Dt_paradis_perCycle, FEM_Du = generate_input_ABAQUS_timeSync(foldername_restart, foldername_ABAQUS, erate, Load_Stage, original_length, paradis_timeNow)
            paradis_timeNow = paradis_timeNow + Dt_paradis_perCycle # Update the current time in ParaDiS
            current_disp = current_disp + FEM_Du # The total displacement applied in FEM
            print(f'\t\t+++ Dt in ParaDiS: {Dt_paradis_perCycle:.6e}      +++ Current Disp {current_disp:.6e}')
        else: # Not synchronize timestep 
            # Get the deltaT in paradis in this cycle and update the applied stress in FEM
            if LoadType == 1:
                # Read the plastic strain rate from ParaDiS
                # epsdot_path = os.path.join(foldername_wsl, foldername_root, foldername_tests, foldername_Job,foldername_properties,'epsdot')
                epsdot_path = os.path.join('..', '..', foldername_tests, foldername_Job,foldername_properties,'epsdot')
                rate_plastic_strain = check_plastic_strain(epsdot_path, Total_N_DD, N_DD)
                rate_plastic_strain = rate_plastic_strain*Vol_dd_fem_ratio
                print(f'+++ Plastic Strain Rate: {rate_plastic_strain:.6e}')
                Dt_paradis_perCycle, FEM_Dsigma = generate_input_ABAQUS_stressControl(
                    foldername_restart,
                    foldername_ABAQUS,
                    rate_plastic_strain,
                    paradis_timeNow,
                    N_DD,
                    critical_eps_rate
                )
            elif LoadType == 2:
                Dt_paradis_perCycle, FEM_Dsigma = generate_input_ABAQUS_erateControl(
                    foldername_restart,
                    foldername_ABAQUS,
                    paradis_timeNow,
                    N_DD,
                    erate,
                    original_length
                )
            else:
                print("  +++ LoadType should be 1 or 2.")
                raise RuntimeError()           


            paradis_timeNow = paradis_timeNow + Dt_paradis_perCycle # Update the current time in ParaDiS
            current_stress = current_stress + FEM_Dsigma # The total displacement applied in FEM

            print(f'\t\t+++ Dt in ParaDiS: {Dt_paradis_perCycle:.6e}      +++ Current Stress (Or Disp) {current_stress:.6e}')



        ## Continue ABAQUS
        
        # ???
        # dt update based on ParaDis deltaTT in restart.cn file
        # if(i>1)
        #     if(systemtype_Linux==1)% For WSL system
        #         [status,cmdout]=system(['bash -c "grep deltaTT ' foldername_root '/' foldername_tests '/' foldername_Job '/' foldername_restart '/restart.cn"']);
        #     else% For Linux system or MacOS
        #         [status,cmdout]=system(['grep deltaTT ' foldername_root '/' foldername_tests '/' foldername_Job '/' foldername_restart '/restart.cn']);
        #     end
        # end


        ## Merge the data  segment files in armdata and delete the segment files(Dis_Segs*) in armdata folder
        # arm_dir = os.path.join(foldername_root, foldername_tests, foldername_Job, foldername_arm)
        arm_dir = os.path.join('..', '..', foldername_tests, foldername_Job, foldername_arm)
        Dis_Segs_all = os.path.join(arm_dir, 'Dis_Segs_all')
        with open(Dis_Segs_all, 'wb') as outfile:
            for fname in sorted(glob.glob(os.path.join(arm_dir, 'Dis_Segs0*'))):
                with open(fname, 'rb') as infile:
                    outfile.write(infile.read())
        
        for fname in glob.glob(os.path.join(arm_dir, 'Dis_Segs0*')):
            os.remove(fname)

    
        ## Copy segment files from armdata folder to ABAQUS directory
        N_segments=generate_segment_files_v1(foldername_ABAQUS,foldername_segs,filename_seg_head,filename_seg_temp,filename_seg_lastCycle,filename_seg_num,dt_FEM_each_step,EffSlip,InelasticStep,i)

        # Empty previous cycle data
        for j in range(1,num_cpus+1):
            filename = f"{filename_surface_int}{j:04d}.txt"
            with open(filename, 'w') as f:
                f.write('')  # Deleting contents
            filename = f"{filename_surface_stress}{j:04d}.txt"
            with open(filename, 'w') as f:
                f.write('')  # Deleting contents
                
        # Move MATLAB flags to ABAQUS flag 
        for j in range(1, num_cpus + 1):
            src = os.path.join(foldername_ABAQUS, f"Matlabworking{j:04d}.flag")
            dst = os.path.join(foldername_ABAQUS, f"ABAQUSworking{j:04d}.flag")
            shutil.move(src, dst)

        print('\t\t+++ Abaqus CONTINUE +++')

        # Wait until ABAQUS fortran file changes ABAQUSworking to Matlabworking
        wait_aba(foldername_ABAQUS, ABAQUS_jobname, num_cpus)


        ## Read the stable time increment from ABAQUS
        with open(os.path.join(foldername_ABAQUS,'FEM_Dt0001.txt'), 'r') as f:
            lines = f.readlines()
            Abaqus_dt = [float(line.strip()) for line in lines if line.strip()]

        with open(os.path.join(foldername_ABAQUS,'FEM_Dt0001.txt'), 'w') as f:
            f.write('')  # Deleting contents
        
        if Abaqus_dt:
            FEM_time_inc = Abaqus_dt[-1]

        print(f'\t\t+++ Dt in ABAQUS : {dt_FEM_each_step:.6e}')


        # Perform post-processing before ParaDis execution.
        print('\t\t+ Post-Processing ABAQUS results... +')

        if DD_Coupling_All_Element == 0:
            generate_postprocessing_data14(Num_Elements,Num_DDM_Elements,Element_type,Elements,Elements_Center,Elements_Min,Elements_Max,
                Jobname,filename_FEMoutput_stress,filename_FEMoutput_nodes,foldername_FEMstress,filename_FEMstress_head,filename_NUCsite_head,foldername_restart,filename_center_element,
                foldername_tests,foldername_temp_win,foldername_ABAQUS,EffSlip,dt_FEM_each_step,filename_rn_check,filename_link_element,
                num_cpus,xyzlimit,SurfaceElementSet,Num_SurfaceElement,FEM_Stress_Scale,i,enabled_temp_disp)
        else:
            generate_postprocessing_data12(Num_Elements,Element_type,Elements,Elements_Center,Elements_Min,Elements_Max,
                Jobname,filename_FEMoutput_stress,filename_FEMoutput_nodes,foldername_FEMstress,filename_FEMstress_head,filename_NUCsite_head,foldername_restart,filename_center_element,
                foldername_tests,foldername_temp_win,foldername_ABAQUS,EffSlip,dt_FEM_each_step,filename_rn_check,filename_link_element,
                num_cpus,xyzlimit,SurfaceElementSet,Num_SurfaceElement,FEM_Stress_Scale,i,enabled_temp_disp)


        ## Change the number of DD time steps for the next ParaDiS run
        print(f'\t\t+ Updating "restart.cn": setting maxstep = {N_DD} +')
        if sync_dd_fem_dt: # Dynamic simulation. Time synchronization between DD and FEM
            dt_FEM_each_step = FEM_time_inc # Update the time step in FEM
            deltaT_DD = Dt_paradis_perCycle/N_DD # the average time increment in DD within this DDM cycle
            N_DD = round(min(FEM_time_inc/deltaT_DD, L_elem/(2*rmax_DD))) # Update the number of time increments for the next DD run
        else: # Quasi-static simulation. Number of DD timesteps per cycle is defined in param.m
             N_DD = round(L_elem/(2*rmax_DD)) # Update the number of time increments for the next DD run

        restart_dir = os.path.join('..', '..', foldername_tests, foldername_Job, foldername_restart)
        restart_file = os.path.join(restart_dir, 'restart.cn')
        command_str = f'sed -i "s/^maxstep.*$/maxstep={N_DD}/" "{restart_file}"'
        subprocess.run(command_str, shell=True, check=True)