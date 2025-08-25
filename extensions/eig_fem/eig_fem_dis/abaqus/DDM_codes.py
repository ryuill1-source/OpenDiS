# (08/14/2025, kyeongmi) Define new module, DDM_codes (should be modified)
import sys, os, time
pydis_paths = ['../../../../python', '../../../../lib', '../../../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]

import numpy as np
from typing import Tuple
from pydis.disnet import DisNet, Tag
from framework.disnet_manager import DisNetManager

class DDM_codes:
    """
    (08/18/2025, kyeongmi) Dummy RemoteStress provider for simulation testing.
    To be replaced with real ABAQUS coupling and flag file I/O implementation.
    """
    def __init__(self, state: dict | None = None) -> None:
        # store config/state so that CalForce can access keys like foldername_ABAQUS
        self.state = dict(state or {})

    def prepare(self, DM=None, state=None):
        """Dummy prepare function (does nothing)."""
        return state

    def run(self, DM=None, state=None):
        """Dummy run function (does nothing)."""
        return state

    def wait(self, flag_check_interval: float = 1.0):
        """Dummy wait function (does nothing)."""
        return
