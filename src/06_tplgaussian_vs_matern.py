"""TPL-Gaussian model vs. Matern model family."""
import os
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import gstools as gs

plt.style.use("default")
mpl.rc("text", usetex=True)
mpl.rc("lines", linewidth=3.5)
plt.close("all")


def dashes(i=1, max_d=12, space=1):
    """Dashes for matplotlib."""
    return i * [space, space] + [max_d - 2 * i * space, space]


def format_ax(axis):
    """Format axis."""
    axis.set_xlim([-0.06, 3.06])
    axis.set_xticks(range(4))
    axis.set_xticklabels(["$0$", r"$\ell$", r"$2\ell$", r"$3\ell$"])
    axis.set_xlabel("distance")
    axis.set_ylabel("semivariance")
    axis.set_ylim([0, 1])
    axis.grid(linestyle=":")
    axis.legend()


save = True
x = np.geomspace(0.01, 3, 10)
grid = np.linspace(0, 10, 100)

# Matern(nu=1.5) vs TPL-Gaussian

fig, ax = plt.subplots(figsize=[5, 3])

m1 = gs.Matern(dim=2, integral_scale=1, nu=1.5)
fit_m1 = gs.TPLGaussian(dim=2)
fit_m1.fit_variogram(x, m1.variogram(x), len_low=0, nugget=0)
m1.plot(ax=ax, x_max=3, label="Matern(nu=1.5)", color="k", linewidth=2)
fit_m1.plot(
    ax=ax, x_max=3, label="TPL-Gaussian(hurst=1.0)", linestyle=":", color="C0"
)

print(m1)
print(fit_m1)

format_ax(ax)
fig.tight_layout()
if save:
    fig.savefig(
        os.path.join("..", "results", "06_matern_tpl_1-5.pdf"), dpi=300
    )
    fig.show()

# Matern(nu=0.5) vs TPL-Gaussian

fig, ax = plt.subplots(figsize=[5, 3])

m2 = gs.Matern(dim=2, integral_scale=1, nu=0.5)
fit_m2 = gs.TPLGaussian(dim=2)
fit_m2.fit_variogram(x, m2.variogram(x), len_low=0, nugget=0)
m2.plot(ax=ax, x_max=3, label="Matern(nu=0.5)", color="k", linewidth=2)
fit_m2.plot(
    ax=ax, x_max=3, label="TPL-Gaussian(hurst=0.45)", linestyle=":", color="C0"
)

print(m2)
print(fit_m2)

format_ax(ax)
fig.tight_layout()
if save:
    fig.savefig(
        os.path.join("..", "results", "07_matern_tpl_0-5.pdf"), dpi=300
    )
    fig.show()

# model families

fig, ax = plt.subplots(figsize=[5, 3])

gau = gs.Gaussian(integral_scale=1)
exp = gs.Exponential(integral_scale=1)
nus = [0.5, 0.75, 1.0, 1.5, 2.0]
# ax = exp.plot(ax=ax, x_max=3, linestyle="-.", color="k")
ax = gau.plot(ax=ax, x_max=3, linestyle="-", color="k")
for i, nu in enumerate(nus):
    m = gs.Matern(integral_scale=1, nu=nu)
    ax = m.plot(
        ax=ax,
        x_max=3,
        label=f"Matern(nu={nu:.2})",
        linewidth=2,
        dashes=dashes(i),
        color="C0",
    )
format_ax(ax)
fig.tight_layout()
if save:
    fig.savefig(os.path.join("..", "results", "08_matern_family.pdf"), dpi=300)
fig.show()

fig, ax = plt.subplots(figsize=[5, 3])

hursts = [0.45, 0.5, 0.6, 0.8, 0.999]
# ax = exp.plot(ax=ax, x_max=3, linestyle="-.", color="k")
ax = gau.plot(ax=ax, x_max=3, linestyle="-", color="k")
for i, hurst in enumerate(hursts):
    m = gs.TPLGaussian(integral_scale=1, hurst=hurst)
    ax = m.plot(
        ax=ax,
        x_max=3,
        label=f"TPL-Gaussian(hurst={hurst:.2})",
        linewidth=2,
        dashes=dashes(i),
        color="C0",
    )
format_ax(ax)
fig.tight_layout()
if save:
    fig.savefig(os.path.join("..", "results", "09_tpl_family.pdf"), dpi=300)
fig.show()

# Fields

fig, ax = plt.subplots(figsize=[5, 4])
srf = gs.SRF(m1, seed=1234)
field1a = srf.structured((grid, grid))
srf.plot(ax=ax, contour_plot=False)
ax.set_title("Matern(nu=1.5)")
fig.tight_layout()
if save:
    fig.savefig(
        os.path.join("..", "results", "10_field_matern_1-5.pdf"), dpi=300
    )
fig.show()

fig, ax = plt.subplots(figsize=[5, 4])
srf = gs.SRF(m2, seed=1234)
srf.structured((grid, grid))
srf.plot(ax=ax, contour_plot=False)
ax.set_title("Matern(nu=0.5)")
fig.tight_layout()
if save:
    fig.savefig(
        os.path.join("..", "results", "11_field_matern_0-5.pdf"), dpi=300
    )
fig.show()
