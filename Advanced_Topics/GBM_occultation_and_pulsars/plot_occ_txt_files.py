import numpy as np
import matplotlib.pyplot as plt
import sys
print sys.argv
filename = sys.argv[1]

from numpy import *
tstart, tstop, flux1, err1, flux2, err2, flux3, err3, flux4, err4, nocc = loadtxt(filename, unpack=True,skiprows=9)

avg1 = sum(flux1/err1**2)/sum(1.0/err1**2)
avg2 = sum(flux2/err2**2)/sum(1.0/err2**2)
avg3 = sum(flux3/err3**2)/sum(1.0/err3**2)
avg4 = sum(flux4/err4**2)/sum(1.0/err4**2)

tmid = (tstart+tstop)/2.0
terr=(tstop-tstart)/2.0
if len(sys.argv) > 2:
        xmin = float(sys.argv[2])
        xmax = float(sys.argv[3])
elif len(sys.argv) <= 2:
        xmin = min(tstart)
        xmax = max(tstop)

print min(tmid), max(tmid), min(flux1), max(flux1)

plt.figure(1)
plt.errorbar(tmid,flux1,xerr=terr,yerr=err1,fmt='o')
plt.plot([xmin,xmax],[0,0])
plt.plot([xmin,xmax],[avg1,avg1],color='m')
plt.xlim(xmin=xmin, xmax=xmax)
plt.ylabel('Flux (mcrab, 12-25 keV)')
plt.title(filename)
plt.xlabel('Time (Modified Julian Date)')

plt.figure(2)
plt.errorbar(tmid,flux2,xerr=terr,yerr=err2,fmt='o')
plt.plot([xmin,xmax],[0,0])
plt.plot([xmin,xmax],[avg2,avg2],color='m')
plt.xlim(xmin=xmin, xmax=xmax)
plt.ylabel('Flux (mcrab, 25-50 keV)')
plt.title(filename)
plt.xlabel('Time (Modified Julian Date)')

plt.figure(3)
plt.errorbar(tmid,flux3,xerr=terr,yerr=err3,fmt='o')
plt.plot([xmin,xmax],[0,0])
plt.plot([xmin,xmax],[avg3,avg3],color='m')
plt.xlim(xmin=xmin, xmax=xmax)
plt.ylabel('Flux (mcrab, 50-100 keV)')
plt.title(filename)
plt.xlabel('Time (Modified Julian Date)')

plt.figure(4)
plt.errorbar(tmid,flux4,xerr=terr,yerr=err4,fmt='o')
plt.plot([xmin,xmax],[0,0])
plt.plot([xmin,xmax],[avg4,avg4],color='m')
plt.xlim(xmin=xmin, xmax=xmax)
plt.ylabel('Flux (mcrab, 100-300 keV)')
plt.title(filename)
plt.xlabel('Time (Modified Julian Date)')
plt.show()
