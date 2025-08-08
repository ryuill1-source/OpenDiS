"""@package docstring
SimulationDriver: class for simulating dislocation network with boundary value problem (BVP)

Provide simulation functions based on other utlitity classes
"""
import sys, os
pydis_paths = ['../../../../python', '../../../../lib', '../../../../core/pydis/python']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]

import numpy as np
import pickle
from pydis.disnet import DisNet
from pydis.calforce.calforce_disnet import CalForce
from pydis.mobility.mobility_disnet import MobilityLaw
from pydis.simulate.sim_disnet import SimulateNetwork
from pydis.timeint.timeint_disnet import TimeIntegration
from pydis.visualize.vis_disnet import VisualizeNetwork
from framework.disnet_manager import DisNetManager

try:
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Line3DCollection
except ImportError:
    print('-----------------------------------------')
    print(' cannot import matplotlib or mpl_toolkits')
    print('-----------------------------------------')

class SimulationDriver(SimulateNetwork):
    """SimulationDriver: class for simulating dislocation network with boundary value problem (BVP)
           extends SimulateNetwork from PyDiS
    """
    def __init__(self, state: dict,
                 #image_stress=None, image_force=None,
                 surface_mobility=None, surface_topology=None,
                 **kwargs) -> None:
        super().__init__(state, **kwargs)

        #self.image_stress = image_stress
        #self.image_force = image_force
        self.surface_mobility = surface_mobility
        self.surface_topology = surface_topology

        # Launch abaqus
        # or put at the beginning of test case

    def step_begin(self, DM: DisNetManager, state: dict):
        """step_begin: invoked at the begining of each time step
        """
        #if self.image_stress is not None:
            # skip first step
            # starting from second step
            # should call FEMRemotStress(DM, state)
            # export eig strain increment (from previous step)
            # remove ABAQUS_pause.flag
            # call function in the remote_stress module
        self.calforce.CalRemoteStress(DM, state)
        pass

    def step_integrate(self, DM: DisNetManager, state: dict):
        """step_integrate: invoked for time-integration at each time step
        """
        #self.save_old_nodes(DM, state)

        state = self.calforce.NodeForce(DM, state)
        state = self.mobility.Mobility(DM, state)
        state = self.timeint.Update(DM, state)
        #self.plastic_strain(DM, state)
        # calculate eigstrain field
        # call function in the remote_stress module

    def step_topological_operations(self, DM: DisNetManager, state: dict):
        """step_topological_operations: invoked for handling topological events at each time step
        """
        if self.cross_slip is not None:
            self.cross_slip.Handle(DM, state)

        # The order of topology vs collision is opposite to ExaDiS
        if self.topology is not None:
            self.topology.Handle(DM, state)

        #ToDo: may need to write a combined Topology class for both bulk and surface nodes
        if self.surface_topology is not None:
            self.surface_topology.Handle(DM, state)

        if self.collision is not None:
            self.collision.HandleCol(DM, state)

        if self.remesh is not None:
            self.remesh.Remesh(DM, state)

    # (2025/08/07, Ill Ryu)
    # For data communication, we need to keep the same dislocation network in plastic strain and rn_links
    # It requires to write a files(Dis_Segs) after topology operations.
    # Output format:
    def write_Dis_Segs(self, DM: DisNetManager, filename_Dis_Segs):
        """Write DisNetManager data to 'Dis_Segs' file
           [Format] NodeID, NBBID, Node_old_position, NBR_old_position, Node_current_position, NBR_current_position, B, N
        """
        data = DM.export_data()
        data['version'] = '1.0'
        data['nodes_attr'] = ['domain', 'index', 'x', 'y', 'z', 'constraint']
        data['segs_attr'] = ['node1', 'node2', 'bx', 'by', 'bz', 'nx', 'ny', 'nz']


    def step_write_files(self, DM: DisNetManager, state: dict):
        if self.write_freq != None:
            istep = state['istep']
            if istep % self.write_freq == 0:
                DM.write_json(os.path.join(self.write_dir, f'disnet_{istep}.json'))
                if self.save_state:
                    with open(os.path.join(self.write_dir, f'state_{istep}.pickle'), 'wb') as file:
                        pickle.dump(state, file)
            # Here we call a new function to print out the data to Dis_Segs
            self.write_Dis_Segs(os.path.join(self.write_dir, f'Dis_Segs_{istep}'))  # create with the step number

    #def my_function(self, DM: DisNetManager, state: dict):
    #    print("my_new_function is called")
    #    pass

