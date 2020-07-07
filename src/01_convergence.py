import os
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from anaflow import ext_theis_tpl, ext_thiem_tpl

plt.style.use('ggplot')
mpl.rc('lines', linewidth=3)

time = 1e4                   # time point for steady state
rad = np.geomspace(0.1, 10)  # radius from the pumping well in [0, 4]
r_ref = 10.0                 # reference radius
KG = 1e-4                    # the geometric mean of the transmissivity
len_scale = 5.0              # correlation length of the log-transmissivity
hurst = 0.5                  # hurst coefficient
var = 0.5                    # variance of the log-transmissivity
rate = -1e-4                 # pumping rate
dim = 1.5                    # using a fractional dimension

head1 = ext_thiem_tpl(rad, r_ref, KG, len_scale, hurst, var, dim=dim, rate=rate)
head2 = ext_theis_tpl(time, rad, 1e-4, KG, len_scale, hurst, var, dim=dim, rate=rate, r_bound=r_ref)
head3 = ext_theis_tpl(time, rad, 1e-4, KG, len_scale, hurst, var, dim=dim, rate=rate)
head3 -= head3[-1]  # quasi-steady

plt.plot(rad, head1, label="Ext Thiem TPL (steady state)")
plt.plot(rad, head2, label="Ext Theis TPL bounded (t={})".format(time), linestyle="--")
plt.plot(rad, head3, label="Ext Theis TPL quasi-steady (t={})".format(time), linestyle=":", color="k")

plt.xlabel("r in [m]")
plt.ylabel("h in [m]")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join("..", "results", "01_ext_theis_tpl_conv.png"), dpi=150)
plt.show()
