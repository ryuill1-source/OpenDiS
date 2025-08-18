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

class CalForce(CalForce_Base):
    """CalForce_DisNet: class for calculating forces on dislocation network
    """
    def __init__(self, state: dict={}, calforce_bulk=None) -> None:
        self.calforce_bulk = calforce_bulk

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
       
        # (08/14/2025, kyeongmi) remove ABAQUS_pause.flag
        foldername_ABAQUS = state.get("foldername_ABAQUS")
        
        ABAQUS_pause = f"{foldername_ABAQUS}/ABAQUS_pause.flag"

        if os.path.exists(ABAQUS_pause):
            os.remove(ABAQUS_pause)
            print(f"remove {ABAQUS_pause} file")
        else:
            print(f"{ABAQUS_pause} does not exist")    

        # (08/14/2025, kyeongmi) since ABAQUS_pause.flag is removed, ABAQUS will start the calculation of stress field (should be modified)
        print("ABAQUS starts the stress calculation ..")
        
        return state

    def ReadRemoteStress(self, DM: DisNetManager, state: dict) -> dict:
        """AddRemoteForce: add remote force from remote stress on nodes

        """
        print("FEMRemoteStress: ReadRemoteStress")
        print("  check whether ABAQUS_stress_ready.flag file exists")
        # Check whether ABAQUS_stress_ready.flag exists
        # if not wait for 250 ms and check again
        # If stress ready, ready stress from file and evaluate remote force
        #  (or read remote force directly from file?)
        return state

    def NodeForce(self, DM: DisNetManager, state: dict, pre_compute: bool=True) -> dict:
        """NodeForce: compute all nodal forces and store them in the state dictionary
        """
        state = self.calforce_bulk.NodeForce(DM, state, pre_compute)

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

