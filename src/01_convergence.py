import os
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from anaflow import ext_theis_tpl, ext_thiem_tpl

plt.style.use('default')
mpl.rc("text", usetex=True)
mpl.rc('lines', linewidth=3.5)
plt.close("all")
fig = plt.figure(figsize=[5, 3])
ax = fig.add_subplot(1, 1, 1)

time = 1e4                   # time point for steady state
rad = np.geomspace(0.1, 9)  # radius from the pumping well in [0, 4]
r_ref = 9.0                 # reference radius
KG = 1e-4                    # the geometric mean of the transmissivity
len_scale = 3.0              # correlation length of the log-transmissivity
hurst = 0.5                  # hurst coefficient
var = 0.5                    # variance of the log-transmissivity
rate = -1e-4                 # pumping rate
dim = 1.5                    # using a fractional dimension

head1 = ext_thiem_tpl(rad, r_ref, KG, len_scale, hurst, var, dim=dim, rate=rate)
head2 = ext_theis_tpl(time, rad, 1e-4, KG, len_scale, hurst, var, dim=dim, rate=rate, r_bound=r_ref)
head3 = ext_theis_tpl(time, rad, 1e-4, KG, len_scale, hurst, var, dim=dim, rate=rate)
head3 -= head3[-1]  # quasi-steady

ax.plot(rad, head1, label="Ext. Thiem TPL (steady state)")
ax.plot(rad, head2, label="Ext. Theis TPL bounded (t={})".format(time), linestyle="--")
ax.plot(rad, head3, label="Ext. Theis TPL quasi-steady (t={})".format(time), linestyle=":", color="k")
ax.set_xticks([0, len_scale, 2 * len_scale, rad[-1]])
ax.set_xticklabels(["$0$", "$\ell$", "$2\ell$", "$3\ell$"])
ax.set_ylim([-1.7, 0.2])
ax.set_xlim([-0.3, 9.3])
ax.set_yticks([-1.5, -1, -0.5, 0])
ax.set_xlabel("$r$ in $[m]$")
ax.set_ylabel("$h(r)$ in $[m]$")
ax.grid(linestyle=':')
ax.legend()
fig.tight_layout()
fig.savefig(os.path.join("..", "results", "01_ext_theis_tpl_conv.pdf"), dpi=300)
fig.show()
