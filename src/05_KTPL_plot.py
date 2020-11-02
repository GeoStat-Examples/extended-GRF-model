# -*- coding: utf-8 -*-
"""
Plot showing the effective well flow conductivity K_TPL(r) with
truncated power law (TPL) variogram in a fractal dimension (d=2.5).

It is compared to its approximation as step function, which serves the
numerical solution of the groundwater flow equation
under radial flow conditions according to the extended Theis TPL method.

K_TPL(r) was derived making use of the upscaling Method Coarse Graining as
extension of the effective well flow conductivity T_CG(r) a with
Gaussian variogram.
It interpolates between the representative conductivity values for pumping
tests flow setting:
    - K_H as near well average representative (considering constant head BC)
    - K_efu (of uniform flow) as far field average representative, which is
        depending on the considered dimension:
        K_efu (2D) = K_G (geometric mean)
        K_efu (3D) = K_G exp(1/6 sigma^2) [for isotropic media]
        K_efu (d)  = K_G exp((1/2 - 1/d) sigma^2) [for isotropic media],
            representing also fractional dimension such as d = 2.5 chosen here
The transition is a function of radial distance to the well r and is controlled
by the stastisitical parameters of the underlying log-normal conducitivity and
its variogram properties: variance $\sigma^2$, correlation length $\all$ and
Hurst coefficient $H$.
"""

import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
from anaflow.tools.coarse_graining import TPL_CG

# from anaflow.tools.mean import annular_hmean
# from anaflow.tools.special import specialrange_cut, specialrange, step_f

###############################################################################
### Plot Settings
plt.style.use("default")
mpl.rc("text", usetex=True)
mpl.rc("lines", linewidth=2.5)
plt.close("all")
fig = plt.figure(figsize=[5, 3])
ax = fig.add_subplot(1, 1, 1)

###############################################################################
### Parameter Settings
cond_gmean = 1e-4  # the geometric mean of the transmissivity
len_scale = 5.0  # correlation length of the log-transmissivity
hurst = 0.5  # hurst coefficient
var = 1.0  # variance of the log-transmissivity

rad = np.linspace(0, 3 * len_scale, 1000)
K_H = np.exp(-var / 2)

###############################################################################
### Plotting
for dim in [3, 2.5, 2]:
    ax.plot(
        rad,
        TPL_CG(rad, cond_gmean, len_scale, hurst, var, dim=dim) / cond_gmean,
        label="$K_{{TPL}}(r)$ (d = {})".format(dim),
    )
    K_efu = np.exp(var * (0.5 - 1 / dim))

    ax.plot(rad, K_efu * np.ones_like(rad), "k--", lw=1)
    ax.text(
        rad[-1],
        1.01 * K_efu,
        "$K_{{\mathrm{{efu}}}}/K_G$ ({}D)".format(dim),
        horizontalalignment="right",
    )

ax.plot(rad, K_H * np.ones_like(rad), "k--", lw=1)
ax.text(
    rad[-1],
    1.01 * K_H,
    "$K_{\mathrm{well}}/K_G$",
    horizontalalignment="right",
)

ax.set_ylim(0.58, 1.25)
ax.set_ylabel("$K(r)/K_G$")
ax.set_xlim([-0.06 * len_scale, 3.06 * len_scale])
ax.set_xticks([0, len_scale, 2 * len_scale, 3 * len_scale])
ax.set_xticklabels(["$0$", "$\ell$", "$2\ell$", "$3\ell$"])
ax.set_xlabel("$r$ in $[m]$")

ax.grid(linestyle=":")
ax.legend(loc=(0.59, 0.21))

fig.tight_layout()
fig.savefig(os.path.join("..", "results", "05_KTPL.pdf"), dpi=300)
fig.show()
