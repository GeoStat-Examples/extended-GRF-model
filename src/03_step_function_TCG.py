import os
import numpy as np
from matplotlib import pyplot as plt
from anaflow.tools.coarse_graining import T_CG, T_CG_error
from anaflow.tools.mean import annular_hmean
from anaflow.tools.special import specialrange_cut, specialrange, step_f

plt.style.use('ggplot')

cond_gmean = 1e-4            # the geometric mean of the transmissivity
len_scale = 5.0              # correlation length of the log-transmissivity
hurst = 0.5                  # hurst coefficient
var = 5.5                    # variance of the log-transmissivity
rate = -1e-4                 # pumping rate
dim = 1.5                    # using a fractional dimension

far_err = 0.01
parts = 30

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

rad = specialrange(0, r_last + 1, 1000)

plt.plot(rad, step_f(rad, R_part, K_part))
plt.plot(rad, T_CG(rad, cond_gmean, var, len_scale))
plt.show()
