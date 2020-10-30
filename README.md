[![GS-Frame](https://img.shields.io/badge/github-GeoStat_Framework-468a88?logo=github&style=flat)](https://github.com/GeoStat-Framework)
[![Gitter](https://badges.gitter.im/GeoStat-Examples/community.svg)](https://gitter.im/GeoStat-Examples/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

# the eGRF model and its application to effective conductivity for TPL variograms


## Description

The *extended* generalized radial flow (eGRF) model is an extension to the GRF model
derived by:

> Barker, J.A., 1988.
> A generalized radial flow model for hydraulic tests infractured rock.
> Water Resources Research 24, 1796–1804. https://doi.org/10.1029/WR024i010p01796

In this workflow, we demonstrate the abilities of the eGRF model and numerically
prove, that the effective transmissivity for truncated power law (TPL) variograms
reproduce the ensemble mean drawdown of pumping tests on synthetic aquifers.


## Structure

Please organize your example in the given Structure
- `src/` - here you should place your python scripts
  - `00_ext_theis_tpl.py` - plotting the effective head for TPL variograms
  - `01_est_run.sh` - bash file running `02_para_estimation.py` in parallel
  - `01_convergence.py` - demonstating the convergence of the effective TPL solution
  - `02_step_function.py` - plotting different step function approximations
  - `03_literature_transmissivities.py` - comparision of drawdowns for different
    transimissivites from literature
  - `comparison`:
    - `00_run_sim_mpi.sh` - bash file running `01_run_sim.py` in parallel
    - `01_run_sim.py` - run all ensemble simulations for pumping tests on TPL aquifers
    - `02_compare_mean.py` - generate comparision plots for the ensemble means
- `results/` - all produced results


## Python environment

Main Python dependencies are stored in `requirements.txt`:

```
gstools==1.2.1
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

MIT © 2020
