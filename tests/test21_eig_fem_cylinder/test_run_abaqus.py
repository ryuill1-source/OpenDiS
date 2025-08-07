''' This python script can run abaqus inp file with user subroutine which are included in subfolder ABAQUS
'''
import numpy as np
import sys, os

pydis_paths = ['../../python', '../../lib', '../../core/pydis/python', '../../extensions/eig_fem']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]
np.set_printoptions(threshold=20, edgeitems=5)

from framework.disnet_manager import DisNetManager
from pydis import DisNode, DisNet, Cell, CellList
from pydis import CalForce as CalForce_Bulk, MobilityLaw as MobilityLaw_Bulk, TimeIntegration, Topology
from pydis import Collision, Remesh, VisualizeNetwork

from eig_fem_dis import SimulationDriver, CalForce, Surface_Topology
from eig_fem_dis import CalForce as CalForce_withSurface, MobilityLaw as MobilityLaw_withSurface

from eig_fem_dis import run_abaqus

# ABAQUS input file and user subroutine file should be placed in foldername_ABAQUS
config = {
        "jobname_head": '',
        "ABAQUS_input_filename": 'test_elastic_VUMAT_control',
        "num_cpus": 1,
        "umatname": 'test_elastic_VUMAT_control.for',
        "foldername_ABAQUS": 'ABAQUS',
        # "cmd_env": 'setenv.bat',
    }

run_abaqus(config)
