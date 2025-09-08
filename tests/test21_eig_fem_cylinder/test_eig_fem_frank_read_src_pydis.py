import numpy as np
import sys, os

# (09/01/2025, kyeongmi) For time measure in DDM
import time
time_elapsed0 = time.time()

pydis_paths = ['../../python', '../../lib', '../../core/pydis/python', '../../extensions/eig_fem']
[sys.path.append(os.path.abspath(path)) for path in pydis_paths if not path in sys.path]
np.set_printoptions(threshold=20, edgeitems=5)

from framework.disnet_manager import DisNetManager
from pydis import DisNode, DisNet, Cell, CellList
from pydis import CalForce as CalForce_Bulk, MobilityLaw as MobilityLaw_Bulk, TimeIntegration, Topology
from pydis import Collision, Remesh, VisualizeNetwork

from eig_fem_dis import SimulationDriver, CalForce, Surface_Topology
from eig_fem_dis.simulate.sim_eig_fem_dis import Nucleation, MAT_TYPE_BCC, MAT_TYPE_FCC #nucleation
from eig_fem_dis import CalForce as CalForce_withSurface, MobilityLaw as MobilityLaw_withSurface

from eig_fem_dis import run_abaqus, DDM_codes
from eig_fem_dis.abaqus.modules.set_params import set_params

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
    net = init_frank_read_src_loop(box_length=Lbox, arm_length=0.4*Lbox, pbc=True)
    nbrlist = CellList(cell=net.cell, n_div=[8,8,8])

    vis = VisualizeNetwork()

    # (08/14/2025, kyeongmi) Since state dictionary has already been defined for ABAQUS initiation,
    #                        one should "update" the state dict. (DO NOT overwrite it) 
    state.update({"burgmag": 3e-10, "mu": 50e9, "nu": 0.3, "a": 1.0, "maxseg": 0.4*Lbox, "minseg": 0.1*Lbox, "rann": 3.0, "mob": 1.0})

    calforce_bulk  = CalForce_Bulk(force_mode='LineTension', state=state)
    mobility_bulk  = MobilityLaw_Bulk(mobility_law='SimpleGlide', state=state)
    timeint   = TimeIntegration(integrator='EulerForward', dt=1.0e-8, state=state)
    collision = Collision(collision_mode='Proximity', state=state, nbrlist=nbrlist)
    remesh    = Remesh(remesh_rule='LengthBased', state=state)

    # (08/14/2025, kyeongmi) import DDM_codes as remote_stress
    remote_stress = DDM_codes(state=state)
    
    calforce = CalForce_withSurface(state=state, calforce_bulk=calforce_bulk)
    mobility = MobilityLaw_withSurface(state=state, mobility_bulk=mobility_bulk)
    # (08//2025, junho) inserted the nucleation module
    nucleation = Nucleation(workdir=os.path.dirname(os.path.abspath(__file__)), dir_femstress = "", stress_filename="SurfaceStress",  material_type=MAT_TYPE_FCC)

    # may need to define new Topology class to use compatible force module and mobility module
    topology  = Topology(split_mode='MaxDiss', state=state, force=calforce_bulk, mobility=mobility_bulk)
    surface_topology  = Surface_Topology(state=state, force=calforce, mobility=mobility)

    sim = SimulationDriver(calforce=calforce, mobility=mobility, timeint=timeint,
                          topology=topology, collision=collision, remesh=remesh, vis=vis,
                          remote_stress=remote_stress,
                          surface_topology=surface_topology, nucleation = nucleation, 
                          state=state, max_step=5, loading_mode="stress",
                          applied_stress=np.array([0.0, 0.0, 0.0, 0.0, 4.0e10, 0.0]),
                          print_freq=1, plot_freq=1, plot_pause_seconds=0.01,
                          write_freq=1, write_dir='output', save_state=False)
    sim.run(net, state)

    return net.is_sane()


if __name__ == "__main__":
    
    # (08/14/2025, kyeongmi) Put keys and values for DDM simulation in the state dictionary (originally in the config dictionary)
    # ABAQUS input file and user subroutine file should be placed in foldername_ABAQUS
    state = {
            # for ABAQUS subprocess execution
            "jobname_head": '',
            "ABAQUS_input_filename": 'test_eig_fem_elastic_VUMAT_control',
            "num_cpus": 1,
            "umatname": 'test_eig_fem_elastic_VUMAT_control.for',
            "foldername_ABAQUS": 'ABAQUS',
            # "cmd_env": 'setenv.bat',

            # used defined parameters for DDM run
            "burger": 2.86e-10,
            "L_elem": 2.86e-10 * 500,
            "density": 2700,
            "E": 70e9,
            "rmax_DD": 50 * 2.86e-10,
            "sync_dd_fem_dt": 0,
            "LoadType": 1,
            "num_max_cycle": 100000,
            "DD_Coupling_All_Element": 1,
            "enabled_temp_disp": 1,
            "plot_freq": 1,
            "critical_eps_rate": 3.0e5,
            "FEM_Volume": 6.7348e10,
            "dd_size": 8000.,
            "DD_Volume": 8000.**3.,
            "Vol_dd_fem_ratio": 8000.**3./6.7348e10,

            "def_ParaDis_DEBUG": 1,
            "def_ABAQUS_DEBUG": 0,

            "foldername_DDM": "DDM_Linux_v0",

            "Element_type": 'C3D8R',
        }

    # (09/01/2025, kyeongmi) Setting up remaining parameters
    # DDM Job name = /'name of this file'
    state["Jobname_fullpath"] = os.path.abspath(__file__)
    state["Jobname"] = os.path.splitext(os.path.basename(__file__))[0]
    
    set_params(state)
    
    # (09/01/2025, kyeongmi) For initiating ABAQUS subprocess. ABAQUS starts and waits until "ABAQUS_running" flag is shown.
    run_abaqus(state)

    main()

    # explore the network after simulation
    G  = net.get_disnet()

    os.makedirs('output', exist_ok=True)
    net.write_json('output/frank_read_src_pydis_final.json')
    
    # (09/01/2025, kyeongmi) For time measure
    time_elapsed = time.time()-time_elapsed0
    print('Total computing time is %e [sec]\n' % time_elapsed)