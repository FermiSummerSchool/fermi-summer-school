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
