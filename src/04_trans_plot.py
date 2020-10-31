"""Generate a TPL field with GSTools."""
import os
import numpy as np
from ogs5py import MSH, specialrange
import gstools as gs
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.tri as tri
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

plt.style.use('default')
mpl.rc("text", usetex=True)
mpl.rc('lines', linewidth=3.5)
plt.close("all")
fig = plt.figure(figsize=[5, 4])
ax = fig.add_subplot(1, 1, 1)

# parameter lists to generate the para_set (single one)
TG = 1e-4  # mu = log(TG)
var = 1.0
len_scale = 10
hurst = 0.5

# radial discretization: 1000 m with 100 steps and increasing stepsize
rad = specialrange(0, 1000, 100, typ="cub")
# 64 angles for discretization
angles = 64
msh = MSH()
msh.generate("radial", dim=2, angles=angles, rad=rad)

# use a rank dependently seeded pseudo-random number generator
seed = 1001
# init cov model (truncated power law with gaussian modes)
cov = gs.TPLGaussian(dim=2, var=var, len_scale=len_scale, hurst=hurst)
# init spatial random field class
srf = gs.SRF(cov, mean=np.log(TG), upscaling="coarse_graining")
# generate new transmissivity field
srf.mesh(msh, seed=seed, point_volumes=msh.volumes_flat)
triang = tri.Triangulation(srf.pos[0], srf.pos[1])
field = srf.field.ravel()

ax.tricontour(triang, field, zorder=-10, levels=16) # anti-alias
cont = ax.tricontourf(triang, field, levels=16)
circle = plt.Circle((0, 0), 990, linewidth=1.5, color='k', fill=False)
ax.add_artist(circle)
cbar = fig.colorbar(cont, ticks=[-12, -10, -8, -6])
cbar.ax.set_ylabel("log-transmissivity")

ax.set_aspect("equal")
ax.set_xticks([-1000, -500, 0, 500, 1000])
ax.set_yticks([-1000, -500, 0, 500, 1000])
ax.set_xlabel("$x$ in $[m]$")
ax.set_ylabel("$y$ in $[m]$")

axins = zoomed_inset_axes(ax, zoom=8.5, loc=1)
axins.tricontour(triang, field, zorder=-10, levels=16) # anti-alias
cont2 = axins.tricontourf(triang, field, levels=16)
axins.scatter(0, 0, 5, color="k")
axins.annotate("well", (5, -2))

axins.set_xlim(-50, 50)
axins.set_ylim(-50, 50)
axins.set_xticks([])
axins.set_yticks([])

mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="k")

# anti-alias
for c in cont.collections:
    c.set_edgecolor("face")
for c in cont2.collections:
    c.set_edgecolor("face")

fig.tight_layout()
fig.savefig(os.path.join("..", "results", "04_trans_plot.pdf"), dpi=300)
fig.show()
