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
# (08/14/2025, kyeongmi) Add arguments remote_stress, remote_force
    def __init__(self, state: dict,
                 remote_stress=None, remote_force=None,
                 surface_mobility=None, surface_topology=None,
                 **kwargs) -> None:
        super().__init__(state, **kwargs)

        self.remote_stress = remote_stress
        self.remote_force = remote_force
        self.surface_mobility = surface_mobility
        self.surface_topology = surface_topology

        # Launch abaqus
        # or put at the beginning of test case

    def step_begin(self, DM: DisNetManager, state: dict):
        """step_begin: invoked at the begining of each time step
           (08/14/2025, kyeongmi) skip the first step and call FEMRemoteStress,
                                  which gets eigenstrain field from OpenDiS and
                                  runs ABAQUS to calculate stress field. 
        """
        print("------step_begin-------")
        if self.remote_stress != None:
            istep = state['istep']
            if istep == 0: # skip the first step
                pass
            else: # starting from second step
            # should call FEMRemotStress(DM, state)
            # export eig strain increment (from previous step)
            # remove ABAQUS_pause.flag
            # call function in the remote_stress module
                self.calforce.FEMRemoteStress(DM, state)
        pass

    def step_integrate(self, DM: DisNetManager, state: dict):
        """step_integrate: invoked for time-integration at each time step
        """
        #self.save_old_nodes(DM, state)

        print("------step_integrate-------")
        state = self.calforce.NodeForce(DM, state)
        state = self.mobility.Mobility(DM, state)
        state = self.timeint.Update(DM, state)
        #self.plastic_strain(DM, state)
        # calculate eigstrain field
        # call function in the remote_stress module

    def step_topological_operations(self, DM: DisNetManager, state: dict):
        """step_topological_operations: invoked for handling topological events at each time step
        """
        print("------step_topological_operations-------")
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
    # We are waiting for the update which could export old_position of nodes.
    def write_Dis_Segs(self, DM: DisNetManager, filename_Dis_Segs: str):
        """Write DisNetManager data to 'Dis_Segs' file
           [Format] NodeID, NBBID, Node_old_position, NBR_old_position, Node_current_position, NBR_current_position, B, N
        """
        # This is copied from disnet.py (get_segs_data_with_positions)
        G = DM.get_disnet(DisNet)
        _, ntags = G.get_nodes_data()
        Nseg = G.num_segments()
        nodeids = np.zeros((Nseg, 2), dtype=int)
        tag1 = np.zeros((Nseg, 2), dtype=int)
        tag2 = np.zeros((Nseg, 2), dtype=int)
        burgers = np.zeros((Nseg, 3))
        planes = np.zeros((Nseg, 3))
        R1 = np.zeros((Nseg, 3))
        R2 = np.zeros((Nseg, 3))

        R1_old = np.zeros((Nseg, 3))
        R2_old = np.zeros((Nseg, 3))
        dyad_flat = np.zeros((Nseg, 9))

        i = 0
        for (source, target), edge_attr in G.all_segments_dict().items():
            nodeids[i,:] = ntags[source], ntags[target]
            tag1[i,:] = source
            tag2[i,:] = target
            burgers[i,:] = edge_attr.burg_vec_from(source).copy()
            planes[i,:] = getattr(edge_attr, "plane_normal", np.zeros(3)).copy()
            r1_local = G.nodes(source).R
            r2_local = G.nodes(target).R
            # apply PBC
            r2_local = G.cell.closest_image(Rref=r1_local, R=r2_local)
            R1[i,:] = r1_local
            R2[i,:] = r2_local
            
            # (08/19/2025, kyeongmi) For testing, I assumed the old position as the current position shifted by -b
            R1_old[i,:] = R1[i,:] - burgers[i,:]
            R2_old[i,:] = R2[i,:] - burgers[i,:]

            dyad = np.outer(planes[i,:], burgers[i,:])
            dyad_flat[i,:] = dyad.ravel()

            i += 1
        '''
        with open(filename_Dis_Segs, "w") as file:
            for i in range(Nseg):
                file.write("%d %d %e %e %e %e %e %e\n"%(nodeids[i,0], nodeids[i,1],
                            R1[i,0], R1[i,1], R1[i,2],
                            R2[i,0], R2[i,1], R2[i,2]))
        '''
        '''
        with open(filename_Dis_Segs, "w") as file:
            for i in range(Nseg):
                file.write("%d %d %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n"
                            %(nodeids[i,0], nodeids[i,1],
                            R1[i,0], R1[i,1], R1[i,2], R2[i,0], R2[i,1], R2[i,2],
                            R1_old[i,0], R1_old[i,1], R1_old[i,2], R2_old[i,0], R2_old[i,1], R2_old[i,2],
                            dyad_flat[i,0], dyad_flat[i,1], dyad_flat[i,2],
                            dyad_flat[i,3], dyad_flat[i,4], dyad_flat[i,5],
                            dyad_flat[i,6], dyad_flat[i,7], dyad_flat[i,8]))
        '''
        # (08/19/2025, kyeongmi) Write Dis_Segs file with desired format
        with open(filename_Dis_Segs, "w") as file:
            for i in range(Nseg):
                values = [
                    nodeids[i,0], nodeids[i,1],
                    R1[i,0], R1[i,1], R1[i,2],
                    R2[i,0], R2[i,1], R2[i,2],
                    R1_old[i,0], R1_old[i,1], R1_old[i,2],
                    R2_old[i,0], R2_old[i,1], R2_old[i,2],
                    *dyad_flat[i,:]   # 9 elements
                ]
                file.write(
                    f"{values[0]:d} {values[1]:d} " + 
                    " ".join(f"{v:.6e}" for v in values[2:]) + "\n"
                )
        print("Dis_Segs file has been written")

    def read_Dis_Segs(self, filename_Dis_Segs: str):
        """read from 'Dis_Segs' file
        """
        # This is copied from disnet.py (get_segs_data_with_positions)
        data = np.loadtxt(filename_Dis_Segs, delimiter=' ', dtype=float)

        return data

    def change_surface_node_flag(self, DM: DisNetManager, state: dict):
        """ Change the flag of node outside the FEM domain
        """
        G = DM.get_disnet(DisNet)
        node1 = G.nodes((0,0))
        node1.constraint = 6 
        return state

    def step_write_files(self, DM: DisNetManager, state: dict):
        print("------step_write_files-------")
        if self.write_freq != None:
            istep = state['istep']
            if istep % self.write_freq == 0:
                DM.write_json(os.path.join(self.write_dir, f'disnet_{istep}.json'))
                if self.save_state:
                    with open(os.path.join(self.write_dir, f'state_{istep}.pickle'), 'wb') as file:
                        pickle.dump(state, file)
            # Here we call a new function to print out the data to Dis_Segs
            self.write_Dis_Segs(DM, os.path.join(self.write_dir, f'Dis_Segs_{istep}'))  # create with the step number

    #def my_function(self, DM: DisNetManager, state: dict):
    #    print("my_new_function is called")
    #    pass

