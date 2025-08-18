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

from eig_fem_dis import run_abaqus, DDM_codes

def init_frank_read_src_loop(arm_length=1.0, box_length=8.0, burg_vec=np.array([1.0,0.0,0.0]), pbc=False):
    '''Generate an initial Frank-Read source configuration
    '''
    print("init_frank_read_src_loop: length = %f" % (arm_length))
    cell = Cell(h=box_length*np.eye(3), is_periodic=[pbc,pbc,pbc])

    rn    = np.array([[0.0, -arm_length/2.0, 0.0,         DisNode.Constraints.PINNED_NODE],
                      [0.0,  0.0,            0.0,         DisNode.Constraints.UNCONSTRAINED],
                      [0.0,  arm_length/2.0, 0.0,         DisNode.Constraints.PINNED_NODE],
                      [0.0,  arm_length/2.0, -arm_length, DisNode.Constraints.PINNED_NODE],
                      [0.0, -arm_length/2.0, -arm_length, DisNode.Constraints.PINNED_NODE]])
    rn[:,0:3] += cell.center()

    N = rn.shape[0]
    links = np.zeros((N, 8))
    for i in range(N):
        pn = np.cross(burg_vec, rn[(i+1)%N,:3]-rn[i,:3])
        pn = pn / np.linalg.norm(pn)
        links[i,:] = np.concatenate(([i, (i+1)%N], burg_vec, pn))

    return DisNetManager(DisNet(cell=cell, rn=rn, links=links))

def main():
    global net, sim, state

    Lbox = 1000.0
    net = init_frank_read_src_loop(box_length=Lbox, arm_length=0.125*Lbox, pbc=True)
    nbrlist = CellList(cell=net.cell, n_div=[8,8,8])

    vis = VisualizeNetwork()

    state.update({"burgmag": 3e-10, "mu": 50e9, "nu": 0.3, "a": 1.0, "maxseg": 0.04*Lbox, "minseg": 0.01*Lbox, "rann": 3.0, "mob": 1.0})

    calforce_bulk  = CalForce_Bulk(force_mode='LineTension', state=state)
    mobility_bulk  = MobilityLaw_Bulk(mobility_law='SimpleGlide', state=state)
    timeint   = TimeIntegration(integrator='EulerForward', dt=1.0e-8, state=state)
    collision = Collision(collision_mode='Proximity', state=state, nbrlist=nbrlist)
    remesh    = Remesh(remesh_rule='LengthBased', state=state)

    # (08/14/2025, kyeongmi) import DDM_codes as remote_stress
    # Put keys and values for DDM simulation in the state dictionary (originally in the config dictionary)
    state.update({"foldername_ABAQUS": 'ABAQUS'})

    remote_stress = DDM_codes(state)
    calforce = CalForce_withSurface(state=state, calforce_bulk=calforce_bulk)
    mobility = MobilityLaw_withSurface(state=state, mobility_bulk=mobility_bulk)

    # may need to define new Topology class to use compatible force module and mobility module
    topology  = Topology(split_mode='MaxDiss', state=state, force=calforce_bulk, mobility=mobility_bulk)
    surface_topology  = Surface_Topology(state=state, force=calforce, mobility=mobility)

    sim = SimulationDriver(calforce=calforce, mobility=mobility, timeint=timeint,
                          topology=topology, collision=collision, remesh=remesh, vis=vis,
                          remote_stress=remote_stress,
                          surface_topology=surface_topology,
                          state=state, max_step=2, loading_mode="stress",
                          applied_stress=np.array([0.0, 0.0, 0.0, 0.0, -4.0e8, 0.0]),
                          print_freq=1, plot_freq=1, plot_pause_seconds=0.01,
                          write_freq=1, write_dir='output', save_state=False)
    sim.run(net, state)

    return net.is_sane()


if __name__ == "__main__":
    
    # (08/14/2025, kyeongmi) Put keys and values for DDM simulation in the state dictionary (originally in the config dictionary)
    # ABAQUS input file and user subroutine file should be placed in foldername_ABAQUS
    state = {
            "jobname_head": '',
            "ABAQUS_input_filename": 'test_eig_fem_elastic_VUMAT_control',
            "num_cpus": 1,
            "umatname": 'test_eig_fem_elastic_VUMAT_control.for',
            "foldername_ABAQUS": 'ABAQUS',
            # "cmd_env": 'setenv.bat',
        }

    # Initiate ABAQUS run, ABAQUS starts and waits until "ABAQUS_running" flag is shown.
    run_abaqus(state)

    main()

    # explore the network after simulation
    G  = net.get_disnet()

    os.makedirs('output', exist_ok=True)
    net.write_json('output/frank_read_src_pydis_final.json')
