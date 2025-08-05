import os

def generate_ParaDiS_folders(foldername_tests,foldername_Job,foldername_segs,foldername_arm,foldername_restart,foldername_properties,foldername_FEMstress):
    folder_paths = [foldername_segs, foldername_arm, foldername_restart, foldername_properties, foldername_FEMstress]
    for folder_path in folder_paths:
        os.makedirs(folder_path, exist_ok=True)