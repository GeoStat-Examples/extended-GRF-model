[![GS-Frame](https://img.shields.io/badge/github-GeoStat_Framework-468a88?logo=github&style=flat)](https://github.com/GeoStat-Framework)
[![Gitter](https://badges.gitter.im/GeoStat-Examples/community.svg)](https://gitter.im/GeoStat-Examples/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4246460.svg)](https://doi.org/10.5281/zenodo.4246460)

# The eGRF model and its application to effective conductivity for TPL variograms


## Description

The *extended* generalized radial flow (eGRF) model is an extension to the GRF model
by allowing radial variable transmissivity and storativity values.
The GRF model was derived by:

> Barker, J.A., 1988.
> A generalized radial flow model for hydraulic tests infractured rock.
> Water Resources Research 24, 1796–1804. https://doi.org/10.1029/WR024i010p01796

In this workflow, we demonstrate the abilities of the eGRF model and numerically
prove, that the effective transmissivity for truncated power law (TPL) variograms
reproduces the ensemble mean drawdown of pumping tests on synthetic aquifers.


## Structure

The workflow is organized by the following structure:
- `src/` - here you should place your python scripts
  - `00_ext_theis_tpl.py` - plotting the effective head for TPL variograms
  - `01_est_run.sh` - bash file running `02_para_estimation.py` in parallel
  - `01_convergence.py` - demonstating the convergence of the effective TPL solution
  - `02_step_function.py` - plot different step function approximations
  - `03_literature_transmissivities.py` - comparision of drawdowns for different
    transimissivites from literature
  - `04_trans_plot.py` - plot a realization of a TPL transmissivity field
  - `05_KTPL_plot.py` - plot K_TPL for different dimensions
  - `06_tplgaussian_vs_matern.py` - comparison of TPL-Gaussian and Matern models
  - `comparison/` - scripts for the comparison of ensemble mean to effective TPL heads
    - `00_run_sim_mpi.sh` - bash file running `01_run_sim.py` in parallel
    - `01_run_sim.py` - run all ensemble simulations for pumping tests on TPL aquifers
    - `02_compare_mean.py` - generate comparision plots for the ensemble means
- `results/` - all produced results


## Python environment

Main Python dependencies are stored in `requirements.txt`:

```
gstools==1.3.0
anaflow==1.0.1
ogs5py==1.1.1
matplotlib
```

You can install them with `pip` (potentially in a virtual environment):

```bash
pip install -r requirements.txt
```


## Contact

You can contact us via <info@geostat-framework.org>.


## License

MIT © 2021
