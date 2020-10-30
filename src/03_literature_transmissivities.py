# -*- coding: utf-8 -*-
import os
import numpy as np
from scipy import stats
from matplotlib import pyplot as plt
from matplotlib import rc
from matplotlib.offsetbox import AnchoredText
from anaflow import ext_grf, ext_grf_steady
from anaflow.tools import neuman2004_trans, K_CG, T_CG
from anaflow.tools import specialrange_cut, annular_hmean, step_f


plt.style.use('default')
rc("text", usetex=True)


def dashes(i=1, max_n=12, width=1):
    """Dashes for matplotlib"""
    return i * [width, width] + [max_n * 2 * width - 2 * i * width, width]

class sol_container:
    """Container to store all information about an eGRF solution."""
    def __init__(self, r_part, t_part, storage, label, title, t_w=None, trans=None):
        self.step_trans = lambda r: step_f(r, r_part, t_part)
        self.step = trans is None
        if self.step:
            self.trans = self.step_trans
        else:
            self.trans = trans
        self.t_w = self.trans(0) if t_w is None else t_w
        self.r_part = r_part
        self.t_part = t_part
        self.s_part = np.full_like(t_part, storage)
        self.label = label
        self.title = title


time_labels = [r"$t=10$s", r"$t=10$min", r"$t=10$h"]
time = [10, 600, 36000]      # 10s, 10min, 10h
rad = np.geomspace(0.1, 9)
rad_lin = np.linspace(0, rad[-1], 1000)
S = 1e-3
K_well = 2e-5
K_far = 1e-4
K_diff = K_far - K_well
K_mean = stats.hmean((K_well, K_far))  # harmonic mean
# variance to result in K_well in 2D
var = -2 * np.log(K_well / K_far)
# K_gmean to result in K_far in 3D with var=1 and anis=1
K_g3d = K_far * np.exp(- 1. / 6.)
len_scale = 3.0
rate = -1e-4
dim = 2

cut_off = 2 * len_scale
parts = 20
r_well = 0.0
r_bound = 50.0
R_part = specialrange_cut(r_well, r_bound, parts, rad[-1])
solutions = []
# SOL 1: GRF with constant conductivity and differing well conductivity
sol = sol_container(
    r_part=[r_well, r_bound],
    t_part=[K_far],
    t_w=K_far,
    storage=S,
    label=r"\textit{Barker (1988)}",
    title="A",
)
solutions.append(sol)
# SOL 2: Butler 1988: two zone aquifer
sol = sol_container(
    r_part=[r_well, len_scale, r_bound],
    t_part=[K_well, K_far],
    t_w=K_well,
    storage=S,
    label=r"\textit{Butler (1988)}",
    title="B",
)
solutions.append(sol)
# SOL 3: Avci 2014: three (multi) zone aquifer
sol = sol_container(
    r_part=[r_well, len_scale, cut_off, r_bound],
    t_part=[K_well, K_mean, K_far],
    t_w=K_well,
    storage=S,
    label=r"\textit{Avci \& Sahin (2014)}",
    title="C",
)
solutions.append(sol)
# SOL 4: Neuman 2004: apparent transmissivity
sol = sol_container(
    trans=lambda r: neuman2004_trans(r, trans_gmean=K_far, var=var, len_scale=len_scale),
    r_part=R_part,
    t_part=annular_hmean(neuman2004_trans, R_part, ann_dim=2, trans_gmean=K_far, var=var, len_scale=len_scale),
    storage=S,
    label=r"\textit{Neuman et al. (2004)}",
    title="D",
)
solutions.append(sol)
# SOL 5: Schneider 2008: extended theis 2D
sol = sol_container(
    trans=lambda r: T_CG(r, trans_gmean=K_far, var=var, T_well=K_well, len_scale=len_scale),
    r_part=R_part,
    t_part=annular_hmean(T_CG, R_part, ann_dim=2, trans_gmean=K_far, var=var, T_well=K_well, len_scale=len_scale),
    storage=S,
    label=r"\textit{Schneider \& Attinger (2008)}",
    title="E",
)
solutions.append(sol)
# SOL 6: Zech 2012: extended theis 3D
sol = sol_container(
    trans=lambda r: K_CG(r, cond_gmean=K_g3d, var=1, anis=1, K_well=K_well, len_scale=len_scale),
    r_part=R_part,
    t_part=annular_hmean(K_CG, R_part, ann_dim=2, cond_gmean=K_g3d, var=1, anis=1, K_well=K_well, len_scale=len_scale),
    storage=S,
    label=r"\textit{Zech et al. (2012)}",
    title="F",
)
solutions.append(sol)

# PLOTTING
plt.close("all")
fig = plt.figure(figsize=[9, 5], constrained_layout=False)
gs = fig.add_gridspec(2, 1)  # two rows
rows = []
rows.append(gs[0].subgridspec(2, 3, height_ratios=[1, 3], hspace=0, wspace=0.1))
rows.append(gs[1].subgridspec(2, 3, height_ratios=[1, 3], hspace=0, wspace=0.1))
# resulting axis pairs as flat list
axis = []
# generate 6 axis with two rows
for j in range(2):
    for i in range(3):
        ax = []
        ax.append(fig.add_subplot(rows[j][0, i]))
        ax.append(fig.add_subplot(rows[j][1, i], sharex=ax[0]))
        plt.setp(ax[0].get_xticklabels(), visible=False)
        ax[0].set_yticks([K_well, K_far])
        if i > 0:  # only show y axis tick-labels for left axis
            plt.setp(ax[0].get_yticklabels(), visible=False)
            plt.setp(ax[1].get_yticklabels(), visible=False)
        else:
            ax[0].set_yticklabels(["$T_{\mathrm{well}}$", "$T_{\mathrm{far}}$"])
            ax[1].set_ylabel(r"$h(r,t)$ in $[m]$")
        # if j == 1:
        #     ax[1].set_xlabel(r"$r$ in $[m]$")
        if j == 0:  # hide x axis tick-labels for the top row
            plt.setp(ax[1].get_xticklabels(), visible=False)
        ax[1].set_xticks([0, len_scale, cut_off, rad[-1]])
        ax[1].set_xticklabels(["$0$", "$\ell$", "$2\ell$", "$3\ell$"])
        axis.append(ax)
# plot given solutions
for i, sol in enumerate(solutions):
    axis[i][0].plot(rad_lin, sol.step_trans(rad_lin), label="transmissivity $T(r)$", linewidth=3.5, solid_joinstyle="miter")
    if not sol.step:
        axis[i][0].plot(rad_lin, sol.trans(rad_lin), linewidth=1, color="w", alpha=0.6)
    axis[i][0].set_ylim([K_well - 0.15 * K_diff, K_far + 0.15 * K_diff])
    axis[i][0].set_title(sol.label)
    # timesteps
    head1 = ext_grf(time, rad, sol.s_part, sol.t_part, sol.r_part, K_well=sol.t_w, dim=dim, rate=rate)
    head2 = ext_grf_steady(rad, r_bound, sol.trans, dim=dim, rate=rate)
    for k, step in enumerate(time):
        axis[i][1].plot(rad, head1[k], dashes=dashes(i=k+1, max_n=9), alpha=1., label=time_labels[k], linewidth=3, color="k")
    axis[i][1].plot(rad, head2, label="steady state", linewidth=2.5, color="C1", alpha=1, zorder=-10)
    axis[i][1].set_ylim([-2.2, 0.2])
    text_box = AnchoredText(sol.title, frameon=True, loc=4, pad=0.5)
    plt.setp(text_box.patch, facecolor='white', alpha=0.5)
    axis[i][1].add_artist(text_box)
    # axis[i][0].grid(axis="x", linestyle=':')
    axis[i][1].grid(linestyle=':')

    # create legend
    if i == 0:
        handles1, labels1 = axis[i][0].get_legend_handles_labels()
        handles2, labels2 = axis[i][1].get_legend_handles_labels()
# add legend
handles = handles1 + handles2
labels = labels1 + labels2
fig.legend(handles, labels, loc='lower center', ncol=5, handlelength=4.0)
# tighten layout
plt.tight_layout()
fig.subplots_adjust(bottom=0.13)
fig.savefig(
    os.path.join("..", "results", "03_literature_comparison.pdf"), dpi=300
)
plt.show()
