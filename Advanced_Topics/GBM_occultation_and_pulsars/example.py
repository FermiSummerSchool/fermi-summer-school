from OccView import OccView


ov = OccView("VELAX-1_daily_lc.fits")


# Inputs for full Plots:
# tstart
# tstop
# channel(s) = 1 ([0,3,4])
# quality="good/bad"  (not needed if you just want good)

ov.PlotFluxes(56100,56500,0,quality="good",save="fullFlux.pdf")

# Inputs for binned Plots:
# tstart
# tstop
# dt
# channel(s) = 1 ([0,3,4])
# quality="good/bad"  (not needed if you just want good)


ov.PlotBinnedFluxes(56100,56500,20,[1,2,5],quality="good",save = "binnedFlux.pdf")



# For Extracting data :
# Binned data
# tstart
# tstop
# dt
# channel
# quality="good/bad"  (not needed if you just want good)

f=ov.GetBinnedFluxes(56100,56500,20,1)

#to get the time
t = f["time"]

#to get the flux

flux = f["flux"]

#to get the error

err = f["error"]

print t
print
print flux
print
print err


# For Extracting data :
# Full data
# tstart
# tstop
# channel
# quality="good/bad"  (not needed if you just want good)

f=ov.GetFluxes(56100,56500,1)

#to get the time
t = f["time"]

#to get the flux

flux = f["flux"]

#to get the error

err = f["error"]
print t
print
print flux
print
print err
