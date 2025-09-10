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
        # Call ReadRemoteStress function to read FEMSTRESS and get stress tensor,
        # which is saved in self.segment_stress
        state = self.ReadRemoteStress(DM, state)
        
        # (08/25/2025, kyeongmi) Remove the ABAQUS_stress_ready.flag for the next step
        foldername_ABAQUS = state.get("foldername_ABAQUS", None)
        jobname_head = state.get("jobname_head", None)
        ABAQUS_input_filename = state.get("ABAQUS_input_filename", None)
        num_cpus = state.get("num_cpus")
        ABAQUS_jobname = jobname_head + ABAQUS_input_filename

        # (08/25/2025, kyeongmi) Use wait_aba function which is located in abaqus/modules/wait_aba.py
        wait_aba(foldername_ABAQUS, ABAQUS_jobname, num_cpus, flag_type = "stress_ready", flag_check_interval = 1)
        
        ABAQUS_stress_ready = f"{foldername_ABAQUS}/ABAQUS_stress_ready.flag"
        if os.path.exists(ABAQUS_stress_ready):
            os.remove(ABAQUS_stress_ready)
            print(f"remove {ABAQUS_stress_ready} file")
        else:
            print(f"{ABAQUS_stress_ready} does not exist after wait_aba (??)")
        

        # (09/03/2025, kyeongmi) Add remote force to dislocation nodes (To be updated: UserStress function)
        # Access and get information from DisNet class
       
        # (09/08/2025, kyeongmi) The following only works for ExaDisNet Class
        #G = DM.get_disnet(DisNet) 
        #cell = G.cell
        #rn = G.get_nodes_data()["positions"]
        #segs = G.get_segs_data()
        #segsnid = segs["nodeids"]
        #b_ij = segs["burgers"]

        # (09/08/2025, kyeongmi) For DisNet Class
        G = DM.get_disnet(DisNet) 
        data = DM.export_data()
        rn = np.array(data["nodes"]["positions"])         # shape: (N_nodes, 3)
        segsnid = np.array(data["segs"]["nodeids"])       # shape: (N_segs, 2)
        b_ij = np.array(data["segs"]["burgers"])          # shape: (N_segs, 3)
        cell = DM.cell
            
        # Calculate positions of segment middle points
        r1 = np.array(cell.closest_image(Rref=np.array(cell.center()), R=rn[segsnid[:,0]]))
        r2 = np.array(cell.closest_image(Rref=r1, R=rn[segsnid[:,1]]))
        r_ij = r2 - r1
        Rseg = 0.5 * (r1 + r2)        
        
        # Get stress in the midpoint of segment (segment-wise stress)
        segment_midpoint_stress = self.SearchSegmentStress(Rseg, cell, state, seg_node_ids=segsnid)

        # Calculate segment-wise current PK force
        sig_remote = segment_midpoint_stress[:,[0,5,4,5,1,3,4,3,2]].reshape(-1,3,3) # Transform user_stress to 3by3 matrix
        sigb = np.einsum('kij,kj->ki', sig_remote, b_ij) # sigma*burgers for each segments
        fpkremote = 0.5 * np.cross(sigb, r_ij) # equally split calculated PKforce
        
        print("[AddUserStress] PK forces from user stress:")
        print(fpkremote)

        # Add to nodal force
        #f = G.get_forces() # current nodal forces
        f = state.get("nodeforces", np.zeros_like(rn))  # for DisNet
        print("[AddUserStress] Original node forces:")
        print(f)
        
        np.add.at(f, segsnid[:,0], fpkremote)
        np.add.at(f, segsnid[:,1], fpkremote)
        
        print("[AddUserStress] Updated node forces:")
        print(f)

        # store new forces in state dictionnary
        state["nodeforces"] = f

        return state

    def CalEigstrainField(self, DM: DisNetManager, state: dict) -> dict:
        # calculate eigstrain field from dislocation motion from previous step
        print("CalForce: CalEigstrainField")
        return state

    def FEMRemoteStress(self, DM: DisNetManager, state: dict) -> dict:
        """FEMRemoteForce: removes ABAQUS_pause.flag and generates ABAQUS_running.flag

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
        wait_aba(foldername_ABAQUS, ABAQUS_jobname, num_cpus, flag_type = "pause", flag_check_interval = 1)
        
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
    
    # (09/04/2025, kyeongmi) Search segment id and return appropriate stress in that segment 
    def SearchSegmentStress(self, R, cell, state, seg_node_ids=None):
        seg_stress_dict = self.segment_stress
        R = np.atleast_2d(R)
        seg_midpoint_stress = np.zeros((R.shape[0],6)) # xx,yy,zz,yz,xz,xy in Pa
        if seg_node_ids is not None:
            for i, (n1, n2) in enumerate(seg_node_ids):
                key = (n1, n2) if n1 <= n2 else (n2, n1)
                stress = seg_stress_dict.get(key, np.zeros(6)) # fallback to zero stress
                print(f"[UserStress] Segment ({n1},{n2}) → Stress: {stress}")
                seg_midpoint_stress[i, :] = stress
        return seg_midpoint_stress
    
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

