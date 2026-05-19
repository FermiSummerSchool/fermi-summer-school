# CTAO Instrument Response Functions - prod5 version v0.1

**Please check always the CTAO webpage ([https://www.cta-observatory.org/science/cta-performance/](https://www.cta-observatory.org/science/cta-performance/)) for the most recent instrument response functions.**

The CTA Observatory (CTAO) will provide very wide energy range and excellent angular resolution and sensitivity in comparison to any existing gamma-ray detector. Energies down to 20 GeV will allow CTAO to study the most distant objects. Energies up to 300 TeV will push CTAO beyond the edge of the known electromagnetic spectrum, providing a completely new view of the sky. 
This data repository provides access to performance evaluation and instrument response functions (IRFs) for CTA.

- IRF version: prod5 v0.1
- Publication date: Sep 2021
- Archived webpage with performance figures included: [CTAO Performance Description (file Website.md)](Website.md)
- Licence: his work is licensed under a [Creative Commons Attribution 4.0 International License](LICENSE).

## Citation and Acknowledgements:

In cases for which the CTA instrument response functions are used in
a research project, we ask to add the following acknowledgement in
any resulting publication:

“This research has made use of the CTA instrument response functions
provided by the CTA Observatory and Consortium, 
see https://www.cta-observatory.org/science/cta-performance/ (version prod5 v0.1; https://doi.org/10.5281/zenodo.5499840) for more details.”

## Description

### Monte Carlo Simulations:

The performance values are derived from detailed Monte Carlo (MC)
simulations of the CTA instrument based on the CORSIKA air shower
code (v7.71, with the hadronic interaction models QGSjet-II-04 and
URQMD, [1]) and telescope simulation tool sim\_telarray [2]. A power-
law gamma-ray spectrum with photon index 2.62 was assumed in the
calculations, although none of the instrument response functions
(e.g. differential flux sensitivities, effective areas, angular or
energy resolutions) depends on the assumed spectral shape of the
gamma-ray source. Background cosmic-ray spectra of proton and
electron/positron particle types are modelled according to recent
measurements from cosmic-ray instruments.

Nominal telescope pointing is assumed, with all telescopes pointing
directions parallel to each other (performance estimation for other
pointing modes, e.g. divergent pointing will be provided in the
future). Performance estimations are available for three zenith angles
(20 deg, 40 deg, and 60 deg), and for each zenith angle for two different
azimuth angles (corresponding to pointing towards the magnetic North
and South). There are significant performance differences found
between the two azimuthal pointing directions  (especially for the
Northern site) as the impact of the geomagnetic field is large
enough to influence notably the air shower development. For general
studies, the  use of the azimuth-averaged instrument response
functions is recommended.

### Instrument Response Functions (IRFs):

The analysis has been tuned to maximize the performance in terms of
flux sensitivity. The optimal analysis cuts depend on the duration
of the observation, therefore the IRFs are provided for 3 different
observation times, from 0.5 to 50 h. IRFs are provided as binned
histogram or FITS tables. It should be stressed, that the full
potential of CTA in terms of angular and energy resolution is not
revealed by these IRFs, due to the focus on the optimisation for
best flux sensitivity.

In general all histograms are binned with a 0.2-binning on the
logarithmic energy axis (5 bins per decade); some selected
histograms (e.g. effective areas or energy migration matrices) are
provided with a finer binning. Effective area and energy migration
matrix are available in a double version: one for the case in which
there is no a priori knowledge of the true direction of incoming
gamma rays (e.g. for the observation of diffuse sources), and
another for observations of point-like objects (including among the
analysis cuts one on the angle between the true and the
reconstructed gamma-ray direction).

IRFs are provided in ROOT format and as FITS tables.
The FITS tables can be used directly as input to science analysis tools.
The values of the IRFs are identical for the different file format, with one exception: the
angular point-spread function is approximated by a Gaussian function
for the FITS tables, while the ROOT files contain the full distribution.

Telescope layouts are preliminary and subject to change. The following
array layouts (Alpha configuration) have been assumed:

- CTA South with 14 MSTs and 37 SSTs (see [figure](figures/CTA-Performance-prod5-v0.1-South-Alpha-Layout.png))
- CTA North with 4 LSTs and 9 MSTs (see [figure](figures/CTA-Performance-prod5-v0.1-North-Alpha-Layout.png))

File Naming (examples):

- Prod5-North-40deg-AverageAz-4LSTs09MSTs.18000s-v0.1.root: IRF for CTA
Northern site on La Palma, 40 deg zenith angle, azimuth-averaged
pointing, optimised for 50 hours of observation time
- Prod5-South-20deg-AverageAz-14MSTs37SSTs.180000s-v0.1.fits.gz: IRF for CTA
Southern site in Paranal, 20 deg zenith angle, azimuth-averaged
pointing, optimised for 50 hours of observation time

List of files:

FITS format:
- fits/CTA-Performance-prod5-v0.1-North-20deg.FITS.tar.gz
- fits/CTA-Performance-prod5-v0.1-North-40deg.FITS.tar.gz
- fits/CTA-Performance-prod5-v0.1-North-60deg.FITS.tar.gz
- fits/CTA-Performance-prod5-v0.1-South-20deg.FITS.tar.gz
- fits/CTA-Performance-prod5-v0.1-South-40deg.FITS.tar.gz
- fits/CTA-Performance-prod5-v0.1-South-60deg.FITS.tar.gz

ROOT format:
- root/CTA-Performance-prod5-v0.1-North-20deg.tar.gz
- root/CTA-Performance-prod5-v0.1-North-40deg.tar.gz
- root/CTA-Performance-prod5-v0.1-North-60deg.tar.gz
- root/CTA-Performance-prod5-v0.1-South-20deg.tar.gz
- root/CTA-Performance-prod5-v0.1-South-40deg.tar.gz
- root/CTA-Performance-prod5-v0.1-South-60deg.tar.gz

IRFs for subarrays of e.g., MSTs only are in the files named MSTSubArray (similar for all other telescope types.

## References

- [1] https://www.ikp.kit.edu/corsika/
- [2] [Bernloehr, K. 2008, Astroparticle Physics, 30, 149](https://ui.adsabs.harvard.edu/abs/2008APh....30..149B/abstract)

## Acknowledgements

We would like to thank the computing centres that provided resources for the generation of the Prod 5 Instrument Response Functions (IRFs):

- CAMK, Nicolaus Copernicus Astronomical Center, Warsaw, Poland
- CIEMAT-LCG2, CIEMAT, Madrid, Spain
- CYFRONET-LCG2, ACC CYFRONET AGH, Cracow, Poland
- DESY-ZN, Deutsches Elektronen-Synchrotron, Standort Zeuthen, Germany
- GRIF, Grille de Recherche d’Ile de France, Paris, France
- IN2P3-CC, Centre de Calcul de l’IN2P3, Villeurbanne, France
- IN2P3-CPPM, Centre de Physique des Particules de Marseille, Marseille, France
- IN2P3-LAPP, Laboratoire d Annecy de Physique des Particules, Annecy, France
- INFN-FRASCATI, INFN Frascati, Frascati, Italy
- INFN-T1, CNAF INFN, Bologna, Italy
- INFN-TORINO, INFN Torino, Torino, Italy
- MPIK, Heidelberg, Germany
- OBSPM, Observatoire de Paris Meudon, Paris, France
- PIC, port d’informacio cientifica, Bellaterra, Spain
- prague_cesnet_lcg2, CESNET, Prague, Czech Republic
- praguelcg2, FZU Prague, Prague, Czech Republic
- UKI-NORTHGRID-LANCS-HEP, Lancaster University, United Kingdom
