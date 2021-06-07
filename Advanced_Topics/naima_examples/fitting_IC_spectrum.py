
# coding: utf-8

# In[1]:

import numpy as np
import naima
import astropy.units as u
from astropy.io import ascii
from naima.models import Synchrotron, InverseCompton, Bremsstrahlung, ExponentialCutoffPowerLaw
import matplotlib as mpl
import matplotlib.pyplot as plt

#Full example can be found at 
#https://naima.readthedocs.io/en/latest/examples.html#fitting-a-minimal-radiative-model
#copied here for convenience.

#see also the examples here: https://github.com/zblz/naima/tree/master/examples


# In[2]:

#read in data. See https://naima.readthedocs.io/en/latest/dataformat.html#dataformat

data = ascii.read("RXJ1713_HESS_2007.dat")

print(data)


# In[3]:

#model definition
#pars is a list/array of parameters to be fitted.
#the 'data' parameter will contain the energies of the data points,
#which are the energies at which the model will be evaluated.
#return value is the flux (dN/dE/dA/dt) at the energies at which we have data points.

def ElectronIC(pars, data):
    """
    Define particle distribution model, radiative model, and return model flux
    at data energy values
    """

    ECPL = ExponentialCutoffPowerLaw(
        pars[0] / u.eV, 10.0 * u.TeV, pars[1], 10 ** pars[2] * u.TeV
    )
    IC = InverseCompton(ECPL, seed_photon_fields=["CMB"])

    return IC.flux(data, distance=1.0 * u.kpc)

#Priors for the bayesian sampler.
#In this case, set un-informative priors.
#We could add limits or even gaussian priors for the other parameters here.
def lnprior(pars):
    # Limit amplitude to positive domain
    logprob = naima.uniform_prior(pars[0], 0.0, np.inf)
    return logprob


# In[4]:

#initial parameters, sampling, and saving of output.
#This will take a few minutes to run.
#several output files will be created.


## Set initial parameters and labels
#parameters here: 
#electron normalization in eV**-1, electron index, electron cutoff energy in TeV.
p0 = np.array((1e30, 3.0, np.log10(30))) 
labels = ["norm", "index", "log10(cutoff)"]

## Run sampler
sampler, pos = naima.run_sampler(
        data_table=data,
        p0=p0,
        labels=labels,
        model=ElectronIC,
        prior=lnprior,
        nwalkers=32,
        nburn=100,
        nrun=20,
        threads=4,
        prefit=True,
        interactive=False,
)
#Try running with more walkers or increasing the number of runs.

## Save run results
out_root = "RXJ1713_IC_minimal"
naima.save_run(out_root, sampler)

## Save diagnostic plots and results table
naima.save_diagnostic_plots(out_root, sampler, sed=True)
naima.save_results_table(out_root, sampler)


# In[51]:

#Plot the fit results.
#Try confs=None, or try e_range=None

f = naima.plot_fit(sampler, modelidx=0, label=None, sed=True, last_step=False, 
               n_samples=100, confs=[3,2,1], ML_info=True, 
               figure=None, plotdata=True, plotresiduals=True, 
               e_unit=None, e_range=(300*u.GeV,200*u.TeV), e_npoints=200, threads=3, 
               xlabel=None, ylabel=None, ulim_opts={}, errorbar_opts={})
f.savefig("RXJ1713_SED.png" )



#corner plot shows 1-D and 2-D projections of posterior probability.
#Useful to check correlations between parameters, see if the sampling was run
#with enough statistics etc.

#Chain plots show the 1-D projection of the posterior distribution for a given parameter
#as well as the evolution of the parameter during the sampling.




