"""@package docstring
CalRemoteForce: class for calculating forces on dislocation network

Provide force calculation functions given a DisNet object
"""
import sys, os
pydis_paths = ['../../../../python', '../../../../lib', '../../../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]

import numpy as np
from pydis.disnet import DisNet, Tag
from framework.disnet_manager import DisNetManager

class CalRemoteStress():
    """CalForce_DisNet: class for calculating forces on dislocation network
    """
    def __init__(self, state: dict={}) -> None:
        pass

    def CalEigstrainField(self, DM: DisNetManager, state: dict) -> dict:
        return state

    def CalRemoteStress(self, DM: DisNetManager, state: dict) -> dict:
        """AddRemoteForce: add remote force from remote stress on nodes

        """
        print("CalRemoteStress: CalRemoteStress")
        state = self.CalEigstrainField(DM, state)
        # export eig strain increment (from previous step)
        # remove ABAQUS_pause.flag
        # call function in the remote_stress module
        return state

    def ReadRemoteStress(self, DM: DisNetManager, state: dict) -> dict:
        """AddRemoteForce: add remote force from remote stress on nodes

        """
        print("CalRemoteStress: ReadRemoteStress")
        # Check whether ABAQUS_stress_ready.flag exists
        # if not wait for 250 ms and check again
        # If stress ready, ready stress from file and evaluate remote force
        #  (or read remote force directly from file?)
        return state
