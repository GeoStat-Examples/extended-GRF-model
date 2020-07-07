"""Generate a TPL ensemble of drawdowns with ogs5py and GSTools."""
import os
import sys
import numpy as np
from ogs5py import OGS, specialrange, generate_time, by_id
from gstools import SRF, TPLGaussian, transform as tf
from mpi4py import MPI

# rank is the actual core-number, size is total number of cores
rank = MPI.COMM_WORLD.Get_rank()
size = MPI.COMM_WORLD.Get_size()

# size of the ensembles
ens_size = 500
# pumping rate
prate = -1e-4
# parameter lists to generate the para_set (single one)
TG = [1e-4]  # mu = log(TG)
var = [2.25]
len_scale = [10]
S = [1e-4]
hurst = [0.5, 0.99]
para_set = np.array(
    [[s, t, v, l, h] for t in TG for v in var for l in len_scale for s in S for h in hurst]
)

RES = os.path.join("..", "..", "results") if len(sys.argv) < 2 else sys.argv[1]
# ogs configuration
task_root = os.path.abspath(os.path.join(RES, "sim2d"))
pcs_type_flow = "GROUNDWATER_FLOW"
var_name_flow = "HEAD"
model = OGS(task_root=task_root, task_id="model")

# spatio-temporal configuration
# define the time stepping: 2 h with 32 steps and increasing stepsize
time = specialrange(0, 7200, 32, typ="cub")
# radial discretization: 1000 m with 100 steps and increasing stepsize
rad = specialrange(0, 1000, 100, typ="cub")
# 64 angles for discretization
angles = 64

# generate mesh and gli
model.msh.generate("radial", dim=2, angles=angles, rad=rad)
model.gli.generate("radial", dim=2, angles=angles, rad_out=rad[-1])
# add the pumping well
model.gli.add_points(points=[0.0, 0.0, 0.0], names="pwell")

# --------------generate different ogs input classes------------------------- #

model.pcs.add_block(  # set the process type
    PCS_TYPE=pcs_type_flow, NUM_TYPE="NEW"
)
model.mpd.add(name="conductivity")
model.mpd.add_block(  # edit recent mpd file
    MSH_TYPE=pcs_type_flow, MMP_TYPE="PERMEABILITY", DIS_TYPE="ELEMENT",
)
model.mmp.add_block(  # permeability, storage and porosity
    GEOMETRY_DIMENSION=2,
    STORAGE=[1, 1.0e-04],
    PERMEABILITY_TENSOR=["ISOTROPIC", 1.0],
    PERMEABILITY_DISTRIBUTION=model.mpd.file_name,
)
model.bc.add_block(  # set boundary condition
    PCS_TYPE=pcs_type_flow,
    PRIMARY_VARIABLE=var_name_flow,
    GEO_TYPE=["POLYLINE", "boundary"],
    DIS_TYPE=["CONSTANT", 0.0],
)
model.ic.add_block(  # set the initial condition
    PCS_TYPE=pcs_type_flow,
    PRIMARY_VARIABLE=var_name_flow,
    GEO_TYPE="DOMAIN",
    DIS_TYPE=["CONSTANT", 0.0],
)
model.st.add_block(  # set pumping condition at the pumpingwell
    PCS_TYPE=pcs_type_flow,
    PRIMARY_VARIABLE=var_name_flow,
    GEO_TYPE=["POINT", "pwell"],
    DIS_TYPE=["CONSTANT_NEUMANN", prate],
)
model.num.add_block(  # set the parameters for the solver
    PCS_TYPE=pcs_type_flow, LINEAR_SOLVER=[2, 5, 1.0e-14, 1000, 1.0, 100, 4],
)
model.tim.add_block(  # set the TIMESTEPS
    PCS_TYPE=pcs_type_flow, **generate_time(time)
)
model.out.add_block(  # set the outputformat for the whole domain
    PCS_TYPE=pcs_type_flow,
    NOD_VALUES=var_name_flow,
    GEO_TYPE="DOMAIN",
    DAT_TYPE="PVD",
    TIM_TYPE=["STEPS", 1],
)

# --------------run OGS simulation------------------------------------------- #

print("write files")
model.write_input()
np.savetxt(os.path.join(model.task_root, "time.txt"), time)
np.savetxt(os.path.join(model.task_root, "rad.txt"), rad)
np.savetxt(os.path.join(model.task_root, "angles.txt"), [angles])

FAIL = []

for para_no, para in enumerate(para_set):
    print("PARA_SET {:04}".format(para_no))
    model.mmp.update_block(STORAGE=[1, para[0]])
    cov_mod = TPLGaussian(dim=2, var=para[2], len_scale=para[3], hurst=para[4])
    srf = SRF(cov_mod, mean=np.log(para[1]), upscaling="coarse_graining")
    # run the ensemble
    # for i in [158, 116]:
    for i in range(ens_size):
        # parallel running the right jobs on each core
        if (para_no * ens_size + i) % size != rank:
            continue
        # generate new transmissivity field
        srf.mesh(model.msh, seed=i, point_volumes=model.msh.volumes_flat)
        # transfrom to log-normal field
        tf.normal_to_lognormal(srf)
        # add the transmissivity to the ogs project
        model.mpd.update_block(DATA=by_id(srf.field))
        # write the new mpd file
        model.mpd.write_file()
        # set the new output-directory
        model.output_dir = os.path.join(
            model.task_root, "para{:04}".format(para_no), "seed{:04}".format(i)
        )
        print("  run model {:04}".format(i), end=" ")
        success = model.run_model(print_log=False)
        print("  ...success") if success else print("  ...error!")
        if not success:
            FAIL.append(str(para_no) + "_" + str(i))
        # export the generated transmissivity field as vtk
        model.msh.export_mesh(
            os.path.join(model.output_dir, "field.vtu"),
            file_format="vtk",
            cell_data_by_id={"transmissivity": srf.field},
        )
    # save the current parameter set
    np.savetxt(
        os.path.join(model.task_root, "para{:04}".format(para_no), "para.txt"),
        para,
        header="storage, trans_gmean, var, len_scale, hurst",
    )

print("FAILED:", FAIL)
