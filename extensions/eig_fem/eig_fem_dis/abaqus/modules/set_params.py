def set_params(state):
    '''
    (09/01/2025, kyeongmi) [..] Modified to work with OpenDiS.eig_fem
                                More specifically, I changed config to state dictionary and removed unneccessary values
    '''
    # (25/07/07 kyeongmi) import generate_ParaDiS_folders function
    from .generate_ParaDiS_folders import generate_ParaDiS_folders
    
    state["foldername_tests"] = 'tests'

    # Time step control setup
    if state.get("sync_dd_fem_dt", False):  # Dynamic simulation
        state["dt_FEM_each_step"] = 1e-10
    else:  # Quasi-static simulation
        state["dt_FEM_each_step"] = 1.0e-4
        state["DD_maxstep_each_cycle"] = 1

    # ABAQUS setup
    # jobname_head is '' for now
    #state["jobname_head"] = 'Job-'
    state["envname"] = 'abaqus_v6.env'

    umatname = 'VUMAT_DDM_multiple_DDcycle.for'
    if "VUMAT_dualPhase" in state:
        umatname = 'VUMAT_DDM_multiple_DDcycle_dualPhase.for'
    if "VUMAT_indent" in state:
        umatname = 'VUMAT_DDM_multiple_DDcycle_indent.for'
    if "VUMAT_damage" in state:
        umatname = 'VUMAT_DDM_multiple_DDcycle_damage.for'
    
    # umatname is test_eig_fem_elastic_VUMAT_control.for for now
    umatname = 'test_eig_fem_elastic_VUMAT_control.for'
    state["umatname"] = umatname

    state["cmd_env"] = 'setenv.bat' #???
    state["foldername_ABAQUS"] = 'ABAQUS'

    # For dislocation slip, the radius of plastic zone.
    burger = state.get("burger")
    L_elem = state.get("L_elem")
    if burger is not None and L_elem is not None:
        state["EffSlip"] = (L_elem / burger) * 2.
    else:
        raise ValueError(f"burger and L_elem should be defined in state")

    # Surface elements and format: For dislocation nucleation, need to set surface elements 
    state["surface_element_set_name"] = 'elset=Set-SurfaceEl'
    state["surface_node_set_name"] = 'nset=Set-SurfaceNode'
    state["formatSpec_Elset"] = '%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f'
    state["Elset_col_num"] = 16

    # Dislocation segment information files
    state["foldername_segs"] = 'armdata'
    state["foldername_arm"] = 'armdata'
    state["foldername_restart"] = 'restart'
    state["foldername_properties"] = 'properties'
    state["foldername_cmdout"] = 'paradis_output'
    state["filename_seg_head"] = 'Dis_Segs' # For nodal points on dislocation segments

    # Temporary input files for ABAQUS
    state["filename_seg_num"] = 'Data_segs_num.txt'
    state["filename_seg_temp"] = 'Data_segs.txt'
    state["filename_seg_lastCycle"] = 'Data_segs_final.txt'

    # ABAQUS FEM results output setup
    state["filename_FEMoutput_stress"] = 'Output_stress'
    state["filename_FEMoutput_nodes"] = 'Output_nodes'
    state["filename_link_element"] = 'Output_link'
    state["filename_center_element"] = 'Output_element_center'
    state["filename_surface_stress"] = 'Surface_stress'
    state["filename_rn_check"] = 'rn_check'
    state["foldername_FEMstress"] = 'FEMstress'
    state["filename_FEMstress_head"] = 'FEMSTRESS' # For segment stress output
    state["filename_NUCsite_head"] = 'NUCSITE' # for nucleation site output
    state["filename_surface_int"] = 'Element_center'


    # ParaDiS executables
    state["paradis_bulk_exe"] = 'paradis_FEM_DD'
    state["paradis_nuc_exe"] = 'paradis_FEM_DD_Force_Nuc'

    ## (25/07/07 kyeongmi) move parameters from python run file to set_params function
    # set folders
    state["foldername_Job"] = state["Jobname"] + '_results'

    # ABAQUS Input
    state["FEM_Stress_Scale"] = 1000000.

    state["xyzlimit"] = state["dd_size"]/2.
    
    # ParaDiS Input filenames
    state["ParaDiS_Input"] = state["Jobname"] + '.ctrl'
    state["ParaDiS_Print"] = state["Jobname"] + '.out'

    # generate ParaDiS folders
    generate_ParaDiS_folders(
        state["foldername_tests"],
        state["foldername_Job"],
        state["foldername_segs"],
        state["foldername_arm"],
        state["foldername_restart"],
        state["foldername_properties"],
        state["foldername_FEMstress"]
    )
    ##
    return state
