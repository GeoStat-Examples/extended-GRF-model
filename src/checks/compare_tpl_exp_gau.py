import numpy as np
import gstools as gs

x = np.linspace(0, 20)

mod_tpl = gs.TPLGaussian(dim=2, var=2.25, len_scale=10.0)
mod_exp = gs.Exponential(dim=2)
mod_gau = gs.Gaussian(dim=2)
mod_mat = gs.Matern(dim=2)

mod_tpl.hurst = 0.5
y = mod_tpl.variogram(x)
mod_exp.fit_variogram(x, y, nugget=False)
mod_mat.fit_variogram(x, y, nugget=False)

ax1 = mod_tpl.plot(x_max=30)
mod_exp.plot(ax=ax1, x_max=30)
mod_mat.plot(ax=ax1, x_max=30)
print(mod_gau)
print(mod_mat)

mod_tpl.hurst = 0.99
mod_tpl.var = 2.25
y = mod_tpl.variogram(x)
mod_gau.fit_variogram(x, y, nugget=False)
mod_mat.fit_variogram(x, y, nugget=False)

ax2 = mod_tpl.plot(x_max=30)
mod_gau.plot(ax=ax2, x_max=30)
mod_mat.plot(ax=ax2, x_max=30)
print(mod_gau)
print(mod_mat)
