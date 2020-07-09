import os
import numpy as np
from matplotlib import pyplot as plt
from anaflow.tools.coarse_graining import TPL_CG, TPL_CG_error
from anaflow.tools.mean import annular_hmean
from anaflow.tools.special import specialrange_cut, specialrange, step_f

plt.close('all')
plt.style.use('ggplot')

cond_gmean = 1e-4            # the geometric mean of the transmissivity
len_scale = 5.0              # correlation length of the log-transmissivity
hurst = 0.5                  # hurst coefficient
var = 1.25                   # variance of the log-transmissivity
rate = -1e-4                 # pumping rate
dim = 2.5                    # using a fractional dimension

far_err = 0.01               # absolute error for cut-off
parts = 20                   # number of partitions

# genearte rlast from a given relativ-error to farfield-conductivity
r_last = TPL_CG_error(far_err, cond_gmean, len_scale, hurst, var, dim=dim)
# generate the partition points (with cut-off)
R_part = specialrange_cut(0, np.inf, parts + 1, r_last)
# calculate the harmonic mean conductivity values within each partition
K_part = annular_hmean(
    TPL_CG,
    R_part,
    ann_dim=dim,
    cond_gmean=cond_gmean,
    len_scale=len_scale,
    hurst=hurst,
    var=var,
    dim=dim,
)

print(r_last)

rad = specialrange(0, np.ceil(r_last + 5), 1000)

plt.plot(rad, step_f(rad, R_part, K_part))
plt.plot(rad, TPL_CG(rad, cond_gmean, len_scale, hurst, var, dim=dim))
plt.legend()
plt.show()
