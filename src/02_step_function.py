import os
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from anaflow.tools.coarse_graining import TPL_CG, TPL_CG_error, T_CG, T_CG_error
from anaflow.tools.mean import annular_hmean
from anaflow.tools.special import specialrange_cut, specialrange, step_f

plt.style.use('default')
mpl.rc("text", usetex=True)
mpl.rc('lines', linewidth=3.5)
plt.close("all")
fig = plt.figure(figsize=[5, 3])
ax = fig.add_subplot(1, 1, 1)

cond_gmean = 1e-4            # the geometric mean of the transmissivity
len_scale = 5.0              # correlation length of the log-transmissivity
hurst = 0.5                  # hurst coefficient
var = 1.0                    # variance of the log-transmissivity
rate = -1e-4                 # pumping rate
dim = 2.0                    # using a fractional dimension

far_err = 0.01               # absolute error for cut-off
parts = 20                   # number of partitions

# genearte rlast from a given relativ-error to farfield-conductivity
r_last_tpl = TPL_CG_error(far_err, cond_gmean, len_scale, hurst, var, dim=dim)
# generate the partition points (with cut-off)
R_part_tpl = specialrange_cut(0, np.inf, parts + 1, r_last_tpl)
# calculate the harmonic mean conductivity values within each partition
K_part_tpl = annular_hmean(
    TPL_CG,
    R_part_tpl,
    ann_dim=dim,
    cond_gmean=cond_gmean,
    len_scale=len_scale,
    hurst=hurst,
    var=var,
    dim=dim,
)

# genearte rlast from a given relativ-error to farfield-conductivity
r_last = T_CG_error(far_err, cond_gmean, var, len_scale)
# generate the partition points
R_part = specialrange_cut(0, np.inf, parts + 1, r_last)
# calculate the harmonic mean conductivity values within each partition
K_part = annular_hmean(
    T_CG,
    R_part,
    ann_dim=dim,
    trans_gmean=cond_gmean,
    var=var,
    len_scale=len_scale,
)

rad = specialrange(0, 3 * len_scale, 1000)

ax.plot(
    rad,
    step_f(rad, R_part_tpl, K_part_tpl),
    color="C0",
    linewidth=3.5,
    solid_joinstyle="miter",
    label="$T_{\mathrm{TPL}}$ - truncated power law",
)
ax.plot(
    rad,
    TPL_CG(rad, cond_gmean, len_scale, hurst, var, dim=dim),
    color="k",
    linewidth=2,
    alpha=0.6,
)

ax.plot(
    rad,
    step_f(rad, R_part, K_part),
    color="C1",
    linewidth=3.5,
    solid_joinstyle="miter",
    label="$T_{\mathrm{RCG}}$ - gaussian",
)
ax.plot(
    rad,
    T_CG(rad, cond_gmean, var, len_scale),
    color="k",
    linewidth=2,
    alpha=0.6,
)

ax.ticklabel_format(axis='y', style='scientific', scilimits=(-4, -4))
ax.set_xlim([-0.06 * len_scale, 3.06 * len_scale])
ax.set_xticks([0, len_scale, 2 * len_scale, 3 * len_scale])
ax.set_xticklabels(["$0$", "$\ell$", "$2\ell$", "$3\ell$"])

ax.set_xlabel("$r$ in $[m]$")
ax.set_ylabel("$T(r)$ in $[m^2 / s]$")
ax.grid(linestyle=':')
ax.legend()

fig.tight_layout()
fig.savefig(
    os.path.join("..", "results", "02_step_functions.pdf"), dpi=300
)
fig.show()
