import netCDF4 as nc
import numpy as np
import pylab as pl

start = 20
end   = 30
step  = 1

hc = 10.

stats = nc.Dataset("canopy.default.0000000.nc","r")
t  = stats.variables["t"][start:end]
z  = stats.variables["z"][:]
zh = stats.variables["zh"][:]
u  = stats.variables["u"][start:end,:]
ustd = stats.variables["u2"][start:end,:]**.5
vstd = stats.variables["v2"][start:end,:]**.5
wstd = stats.variables["w2"][start:end,:]**.5

uflux = stats.variables["uflux"][start:end,:]
vflux = stats.variables["vflux"][start:end,:]

ih = (np.abs(zh-hc)).argmin()
ustar = (uflux[:,ih]**2 + vflux[:,ih]**2)**.25

stats.close()

pl.rc('font',**{'family':'serif','serif':['Palatino']})
pl.rc('text', usetex=True)

pl.figure()
pl.subplot(221)
for n in range(t.size):
    pl.plot(ustd[n,:]/ustar[n], z/hc, color='#cccccc')
pl.xlabel(r'$\sigma_u / u_*$')
pl.ylabel(r'$z / h_c$')
pl.ylim(0, 6)

pl.subplot(222)
for n in range(t.size):
    pl.plot(vstd[n,:]/ustar[n], z/hc, color='#cccccc')
pl.xlabel(r'$\sigma_v / u_*$')
pl.ylabel(r'$z / h_c$')
pl.ylim(0, 6)

pl.subplot(223)
for n in range(t.size):
    pl.plot(wstd[n,:]/ustar[n], zh/hc, color='#cccccc')
pl.xlabel(r'$\sigma_w / u_*$')
pl.ylabel(r'$z / h_c$')
pl.ylim(0, 6)

pl.subplot(224)
for n in range(t.size):
    pl.plot(uflux[n,:]/ustar[n]**2, zh/hc, color='#cccccc')
pl.xlabel(r'$\overline{u^\prime w^\prime} / u_*^2$')
pl.ylabel(r'$z / h_c$')
pl.ylim(0, 6)

