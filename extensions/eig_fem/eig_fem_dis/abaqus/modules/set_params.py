def set_params(config):
    # (25/07/07 kyeongmi) import generate_ParaDiS_folders function
    from generate_ParaDiS_folders import generate_ParaDiS_folders
    
    config["foldername_tests"] = 'tests'

    # Time step control setup
    if config.get("sync_dd_fem_dt", False):  # Dynamic simulation
        config["dt_FEM_each_step"] = 1e-10
    else:  # Quasi-static simulation
        config["dt_FEM_each_step"] = 1.0e-4
        config["DD_maxstep_each_cycle"] = 1

    # ABAQUS setup
    config["jobname_head"] = 'Job-'
    config["envname"] = 'abaqus_v6.env'

    umatname = 'VUMAT_DDM_multiple_DDcycle.for'
    if "VUMAT_dualPhase" in config:
        umatname = 'VUMAT_DDM_multiple_DDcycle_dualPhase.for'
    if "VUMAT_indent" in config:
        umatname = 'VUMAT_DDM_multiple_DDcycle_indent.for'
    if "VUMAT_damage" in config:
        umatname = 'VUMAT_DDM_multiple_DDcycle_damage.for'
    config["umatname"] = umatname

    config["cmd_env"] = 'setenv.bat' #???
    config["foldername_ABAQUS"] = 'ABAQUS'

    config["filename_merge_seg_shell"] = 'merge_DisSegs.sh' #???

    # For dislocation slip, the radius of plastic zone.
    burger = config.get("burger")
    L_elem = config.get("L_elem")
    if burger is not None and L_elem is not None:
        config["EffSlip"] = (L_elem / burger) * 2.
    else:
        raise ValueError(f"burger and L_elem should be defined in config")

    # Surface elements and format: For dislocation nucleation, need to set surface elements 
    config["surface_element_set_name"] = 'elset=Set-SurfaceEl'
    config["surface_node_set_name"] = 'nset=Set-SurfaceNode'
    config["formatSpec_Elset"] = '%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f'
    config["Elset_col_num"] = 16

    # Dislocation segment information files
    config["foldername_segs"] = 'armdata'
    config["foldername_arm"] = 'armdata'
    config["foldername_restart"] = 'restart'
    config["foldername_properties"] = 'properties'
    config["foldername_cmdout"] = 'paradis_output'
    config["filename_seg_head"] = 'Dis_Segs' # For nodal points on dislocation segments

    # Temporary input files for ABAQUS
    config["filename_seg_num"] = 'Data_segs_num.txt'
    config["filename_seg_temp"] = 'Data_segs.txt'
    config["filename_seg_lastCycle"] = 'Data_segs_final.txt'

    # ABAQUS FEM results output setup
    config["filename_FEMoutput_stress"] = 'Output_stress'
    config["filename_FEMoutput_nodes"] = 'Output_nodes'
    config["filename_link_element"] = 'Output_link'
    config["filename_center_element"] = 'Output_element_center'
    config["filename_surface_stress"] = 'Surface_stress'
    config["filename_rn_check"] = 'rn_check'
    config["foldername_FEMstress"] = 'FEMstress'
    config["filename_FEMstress_head"] = 'FEMSTRESS' # For segment stress output
    config["filename_NUCsite_head"] = 'NUCSITE' # for nucleation site output
    config["filename_surface_int"] = 'Element_center'


    # ParaDiS executables
    config["paradis_bulk_exe"] = 'paradis_FEM_DD'
    config["paradis_nuc_exe"] = 'paradis_FEM_DD_Force_Nuc'

    ## (25/07/07 kyeongmi) move parameters from python run file to set_params function
    # set folders
    config["foldername_Job"] = config["Jobname"] + '_results'

    # ABAQUS Input
    config["FEM_Stress_Scale"] = 1000000.

    config["xyzlimit"] = config["dd_size"]/2.
    
    # ParaDiS Input filenames
    config["ParaDiS_Input"] = config["Jobname"] + '.ctrl'
    config["ParaDiS_Print"] = config["Jobname"] + '.out'

    # generate ParaDiS folders
    generate_ParaDiS_folders(
        config["foldername_tests"],
        config["foldername_Job"],
        config["foldername_segs"],
        config["foldername_arm"],
        config["foldername_restart"],
        config["foldername_properties"],
        config["foldername_FEMstress"]
    )
    ##
    return config
