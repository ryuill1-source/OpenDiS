"""@package docstring
CalForce: class for calculating forces on dislocation network

Provide force calculation functions given a DisNet object
"""
import sys, os
pydis_paths = ['../../../../python', '../../../../lib', '../../../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]

import numpy as np
from typing import Tuple
from pydis.disnet import DisNet, Tag
from framework.disnet_manager import DisNetManager
from framework.calforce_base import CalForce_Base
from eig_fem_dis.abaqus.modules.wait_aba import wait_aba

class CalForce(CalForce_Base):
    """CalForce_DisNet: class for calculating forces on dislocation network
    """
    def __init__(self, state: dict={}, calforce_bulk=None) -> None:
        self.calforce_bulk = calforce_bulk
        self.segment_stress: dict[tuple[int,int], np.ndarray] = {} # Save data from FEMSTRESS
                                                                   # Format: Dictionary (key: nodes at both ends, value: stress as voigt order)

    def AddRemoteForce(self, DM: DisNetManager, state: dict) -> dict:
        """AddRemoteForce: add image force from image stress on nodes
           We could use 'userstress' module by Kyeongmi

        """
        print("CalForce: AddRemoteForce")
        state = self.ReadRemoteStress(DM, state)
        # Add remote force to dislocation nodes
        return state

    def CalEigstrainField(self, DM: DisNetManager, state: dict) -> dict:
        # calculate eigstrain field from dislocation motion from previous step
        print("CalForce: CalEigstrainField")
        return state

    def FEMRemoteStress(self, DM: DisNetManager, state: dict) -> dict:
        """AddRemoteForce: add remote force from remote stress on nodes

        """
        # export eig strain increment (from previous step)
        print("CalForce: FEMRemoteStress")
        state = self.CalEigstrainField(DM, state)
       
        foldername_ABAQUS = state.get("foldername_ABAQUS", None)
        jobname_head = state.get("jobname_head", None)
        ABAQUS_input_filename = state.get("ABAQUS_input_filename", None)
        num_cpus = state.get("num_cpus")
        ABAQUS_jobname = jobname_head + ABAQUS_input_filename

        # (08/19/2025, kyeongmi) Use wait_aba function which is located in abaqus/modules/wait_aba.py
        wait_aba(foldername_ABAQUS, ABAQUS_jobname, num_cpus, flag_type = "stress_ready", flag_check_interval = 1)
        wait_aba(foldername_ABAQUS, ABAQUS_jobname, num_cpus, flag_type = "pause", flag_check_interval = 1)
        
        # (08/19/2025, kyeongmi) Remove the ABAQUS_stress_ready.flag for the next step
        ABAQUS_stress_ready = f"{foldername_ABAQUS}/ABAQUS_stress_ready.flag"
        if os.path.exists(ABAQUS_stress_ready):
            os.remove(ABAQUS_stress_ready)
            print(f"remove {ABAQUS_stress_ready} file")
        else:
            print(f"{ABAQUS_stress_ready} does not exist after wait_aba (??)")
        
        # (08/14/2025, kyeongmi) remove ABAQUS_pause.flag
        ABAQUS_pause = f"{foldername_ABAQUS}/ABAQUS_pause.flag"
        if os.path.exists(ABAQUS_pause):
            os.remove(ABAQUS_pause)
            print(f"remove {ABAQUS_pause} file")
        else:
            print(f"{ABAQUS_pause} does not exist after wait_aba (??)")    

        # (08/19/2025, kyeongmi) ABAQUS starts the calculation of stress field.
        ABAQUS_running = f"{foldername_ABAQUS}/ABAQUS_running.flag"
        with open(ABAQUS_running, "w") as f:
            f.write("")    
        print(f"created {ABAQUS_running} file, ABAQUS starts the stress calculation ..")
        
        return state

    def ReadRemoteStress(self, DM: DisNetManager, state: dict) -> dict:
        """AddRemoteForce: add remote force from remote stress on nodes

        """
        print("FEMRemoteStress: ReadRemoteStress")
        print("Check whether ABAQUS_stress_ready.flag file exists")
        
        # (08/18/2025, kyeongmi) Check whether ABAQUS_stress_ready.flag exists
        # if not wait for 250 ms and check again (To test, 3 s for now)
        foldername_ABAQUS = state.get("foldername_ABAQUS", None)
        jobname_head = state.get("jobname_head", None)
        ABAQUS_input_filename = state.get("ABAQUS_input_filename", None)
        num_cpus = state.get("num_cpus")
        ABAQUS_jobname = jobname_head + ABAQUS_input_filename

        # (08/18/2025, kyeongmi) Use wait_aba function which is located in abaqus/modules/wait_aba.py
        wait_aba(foldername_ABAQUS, ABAQUS_jobname, num_cpus, flag_check_interval = 1)
        
        # TODO: If stress ready, read stress from file and evaluate remote force
        #       (or read remote force directly from file?)
        # (08/18/2025, kyeongmi) Load FEMSTRESS file and save in self.segment_stress
        seg_stress = {}
        try:
            with open("FEMSTRESS","r") as f:
                for line in f:
                    if not line.strip(): # skip empty lines
                        continue
                    s = line.strip().split()
                    if len(s) < 9:
                        continue # skip malformed line
                    
                    n1, n2 = int(s[0]), int(s[1])
                    # Use (min, max) tuple as the key -> easier for serialization/debugging
                    key = (n1, n2) if n1 <= n2 else (n2, n1)
                    
                    stress_voigt = np.array([float(s[2]), float(s[3]), float(s[4]),
                                            float(s[5]), float(s[6]), float(s[7])])
                    seg_stress[key] = stress_voigt
                print("FEMSTRESS file has been read")
        except FileNotFoundError:
            print("FEMSTRESS file not found. User stress will default to zero.")
            seg_stress = {}

        self.segment_stress = seg_stress

        # (08/18/2025, kyeongmi) To check whether FEMSTRESS has been loaded properly or not 
        for key, stress in self.segment_stress.items():
            print(f"Segment {key}: {stress}")
        
        return state

    def NodeForce(self, DM: DisNetManager, state: dict, pre_compute: bool=True) -> dict:
        """NodeForce: compute all nodal forces and store them in the state dictionary
        """
        istep = state['istep']
        state = self.calforce_bulk.NodeForce(DM, state, pre_compute)

        if istep == 0: # skip the first step
            pass
        else: # starting from second step
            state = self.AddRemoteForce(DM, state)
    
        return state

    def PreCompute(self, DM: DisNetManager, state: dict) -> dict:
        """PreCompute: pre-compute some data for force calculation
        """
        return state

    def OneNodeForce(self, DM: DisNetManager, state: dict, tag: Tag, update_state: bool=True) -> np.array:
        """OneNodeForce: compute and return the force on one node specified by its tag
        """
        state = self.calforce_bulk.OneNodeForce(DM, state, tag, update_state)

        #ToDo: Add Remote force contribution to OneNodeForce
        return state
    

    def NodeForce_Elasticity_SBA_Cutoff(self, G: DisNet, applied_stress: np.ndarray) -> Tuple[dict, dict]:
        """ It is copied from NodeForce_Elasticity_SBA(calforce_disnet.py)
            We want to compute the elastic interation within cutoff.
            Neighbour list needs to be useful (pydis/nbrlist)
        """
        return ({},{})

