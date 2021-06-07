# Documentation

* naima webpage with many examples: https://naima.readthedocs.io/en/latest/index.html 
* naima github repo with more example scripts and example data: https://github.com/zblz/naima
* Documentation and more info about posterior sampling using emcee: http://dfm.io/emcee/current/

# Prerequisites

Naima can be installed using `pip` or `conda`:
    
    pip install naima

or

    conda install naima

Recent-ish versions of matplotlib and astropy are needed.

# Examples provided here

`exploring_pion_decay_spectra.ipynb` and `exploring_electron_spectra.ipynb` showcase how to produce emission spectra given an underlying population of protons or electrons. Feel free to play around with the model parameters to see how the emission spectra change.

`fitting_IC_spectrum.ipynb` is a very minimal example of how to fit parameters of an underlying electron population to VHE gamma-ray measurements.

All scripts are also provided as standalone python scripts that can be run outside of  jupyter notebooks.

# For Fermibottle users

As of summer 2021, `naima 0.8.4` is installed as part of the fermit bottle. Unfortunately, there is a conflict in one of the depdencencies, `emcee`. Naima needs an older version: `emcee<3.0,>=2.2.0`, whereas `pint-pulsar 0.7` (also part of the fermi bottle) requires `emcee>=3.0`. If you just want to run the "exploring ..." scripts, this does not matter. If you want to run the `fitting_IC_spectrum` script or notebook, you will have to downgrade `emcee`. One way to do that is to run `pip install --upgrade naima`. This will install `naima 0.9.1` and `emcee 2.2.1`. Unfortunately, this might cause issues with `threeML` and `pint-pulsar`, both of which require newer versions of `emcee`. If you want to use `emcee` with either of those packages, you will have to create separate conda environments with the old and new versions of `naima`. 
