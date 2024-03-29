"""Post processing the TPL ensembles."""
import os
import glob
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm, rc
import anaflow as ana
from PyPDF2 import PdfFileMerger


rc("text", usetex=True)
CWD = os.path.abspath(os.path.join("..", "..", "results"))


def calc_ensemble_mean(base="eGRF_TPL_2D", p_min=0, p_max=np.inf):
    """Generate mean for all simulations in ensemble from single means."""
    para_sets = sorted(glob.glob(os.path.join(CWD, base, "para*")))
    time = np.loadtxt(os.path.join(CWD, base, "time.txt"))
    rad = np.loadtxt(os.path.join(CWD, base, "rad.txt"))
    # iterate over all parameter sets
    for para_no, para_set in enumerate(para_sets):
        if para_no < p_min or para_no > p_max:
            continue
        print(para_no, "PARA_SET: ensemble mean calculation")
        ensemble = sorted(glob.glob(os.path.join(para_set, "seed*")))
        rt_head = np.zeros(time.shape + rad.shape, dtype=float)
        # collect all ensemble members
        for cnt, single in enumerate(ensemble, start=1):
            sgl_head = np.loadtxt(os.path.join(single, "rad_mean_head.txt"))
            if cnt % 50 == 0:
                print("   para-Set", para_no, ": member", cnt)
            rt_head += sgl_head
        # ensemble mean
        if ensemble:  # skip if ensemble is empty
            print(" --> write file 'rad_mean_head.txt'")
            rt_head /= len(ensemble)
            np.savetxt(os.path.join(para_set, "rad_mean_head.txt"), rt_head)


def compare(base="eGRF_TPL_2D", p_min=0, p_max=np.inf):
    """Compare ensemble mean to effective drawdown solution."""
    path = os.path.join(CWD, base)
    para_sets = sorted(glob.glob(os.path.join(path, "para*")))
    time = np.loadtxt(os.path.join(path, "time.txt"))
    rad = np.loadtxt(os.path.join(path, "rad.txt"))
    time_range = time > 60
    rad_range = np.logical_and(rad > 0.2, rad < 40)
    time_select = time[time_range]
    rad_select = rad[rad_range]
    # iterate over all parameter sets
    for para_no, para_set in enumerate(para_sets):
        if para_no < p_min or para_no > p_max:
            continue
        para = np.loadtxt(os.path.join(para_set, "para.txt"))
        rt_head = np.loadtxt(os.path.join(para_set, "rad_mean_head.txt"))
        print(para_no, "PARA_SET")
        rt_head = rt_head[time_range]
        rt_head = rt_head[:, rad_range]
        et_head = ana.ext_theis_tpl(  # effective transmissivity
            time=time_select,
            rad=rad_select,
            storage=para[0],
            cond_gmean=para[1],
            var=para[2],
            len_scale=para[3],
            hurst=para[4],
            rate=-1e-4,
            parts=30,
            # gaussian covmodel in GSTools normalized to integral scale
            prop=np.sqrt(np.pi * 2),
        )
        # print relative errors
        abs_diff = np.abs(rt_head - et_head)
        abs_mean = 0.5 * (np.abs(rt_head) + np.abs(et_head))
        rel_diff = abs_diff / np.max(abs_mean)
        print("  max rel. diff:", np.max(rel_diff))
        # creat plots
        plot_diff(
            time_select, rad_select, rt_head, et_head, para_no, path, para
        )
    # merge pdfs
    pdfs = sorted(glob.glob(os.path.join(path, "*_diff.pdf")))
    merger = PdfFileMerger()
    for pdf in pdfs:
        merger.append(pdf)
    merger.write(os.path.join(path, "diff.pdf"))
    # for pdf in pdfs:
    #     os.remove(pdf)


def plot_diff(time, rad, rt_head, et_head, para_no, path, para):
    """Plot the comparisson between effective head and ensemble mean."""
    plt.close("all")
    fig = plt.figure(figsize=[10, 3.4])
    ax0 = plt.subplot2grid((1, 3), (0, 0), fig=fig, projection=Axes3D.name)
    ax1 = plt.subplot2grid((1, 3), (0, 2), fig=fig, projection=Axes3D.name)
    ax2 = plt.subplot2grid((1, 3), (0, 1), fig=fig, projection=Axes3D.name)

    time_m, rad_m = np.meshgrid(time, rad, indexing="ij")
    diff = np.abs(rt_head - et_head)
    z_max = 0.1
    z_min = -2.1

    ax0.plot_surface(
        rad_m,
        time_m,
        rt_head,
        rstride=1,
        cstride=1,
        cmap=cm.RdBu,
        linewidth=0.3,
        antialiased=True,
        edgecolors="k",
    )
    ax1.plot_surface(
        rad_m,
        time_m,
        et_head,
        rstride=1,
        cstride=1,
        cmap=cm.RdBu,
        linewidth=0.3,
        antialiased=True,
        edgecolors="k",
    )
    ax2.plot_surface(
        rad_m,
        time_m,
        diff,
        rstride=1,
        cstride=1,
        cmap=cm.RdBu_r,
        linewidth=0.3,
        antialiased=True,
        edgecolors="k",
        vmin=0, vmax=1,
    )
    ax2.plot_surface(
        rad_m,
        time_m,
        np.zeros_like(diff),
        rstride=1,
        cstride=1,
        color="k",
        alpha=0.4,
        antialiased=True,
    )
    ax0.view_init(elev=15, azim=-150)
    ax1.view_init(elev=15, azim=-150)
    ax2.view_init(elev=15, azim=-150)
    fig.suptitle(
        r"Parameter set P{}: ".format(para_no)
        + r"$S={:.1e}".format(para[0])
        + r"$, $T_G={:.1e}".format(para[1])
        + r"$, $\sigma^2={}".format(para[2])
        + r"$, $\ell={}".format(para[3])
        + r"$, $H={}".format(para[4])
        + r"$"
    )
    xlab = r"$r$ in $\mathrm{[m]}$"
    ylab = r"$t$ in $\mathrm{[s]}$"
    ax0.set_title("Ensemble mean drawdown in $\mathrm{[m]}$", pad=-5)
    ax1.set_title("Effective drawdown in $\mathrm{[m]}$", pad=-5)
    ax2.set_title("Absolute difference in $\mathrm{[m]}$", pad=-5)
    # change axes labels
    ax0.set_xlabel(xlab)
    ax0.set_ylabel(ylab)
    ax1.set_xlabel(xlab)
    ax1.set_ylabel(ylab)
    ax2.set_xlabel(xlab)
    ax2.set_ylabel(ylab)
    ax0.set_zlim((z_min, z_max))
    ax1.set_zlim((z_min, z_max))
    ax2.set_zlim((-1, 1))
    fig.tight_layout()
    plt.savefig(os.path.join(path, "{:04}_diff.pdf".format(para_no)), dpi=300)
    plt.close("all")


calc_ensemble_mean()
compare()
