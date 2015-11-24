import netCDF4 as nc
import numpy as np
import pylab as pl

end = 2400
hc = 10.

stats = nc.Dataset("canopy.default.0000000.nc","r")
start = end-900

t  = stats.variables["t"][start:end]
z  = stats.variables["z"][:]
zh = stats.variables["zh"][:]
u  = stats.variables["u"][start:end,:]
ustd = stats.variables["u2"][start:end,:]**.5
vstd = stats.variables["v2"][start:end,:]**.5
wstd = stats.variables["w2"][start:end,:]**.5

uw = stats.variables["uw"][start:end,:]
udiff = stats.variables["udiff"][start:end,:]
uflux = stats.variables["uflux"][start:end,:]
vflux = stats.variables["vflux"][start:end,:]
ugrad = stats.variables["ugrad"][start:end,:]
stats.close()

#ih = (np.abs(zh-hc)).argmin()
#ustar = (uflux[:,ih]**2 + vflux[:,ih]**2)**.25

ustar = np.zeros(t.size)
for n in range(t.size):
    ustar[n] = np.interp(hc, zh, (uflux[n,:]**2 + vflux[n,:]**2)**.25)

u_avg = np.mean(u, 0)
ustd_avg = np.mean(ustd, 0)
vstd_avg = np.mean(vstd, 0)
wstd_avg = np.mean(wstd, 0)
uflux_avg = np.mean(uflux, 0)
uw_avg = np.mean(uw, 0)
udiff_avg = np.mean(udiff, 0)
ugrad_avg = np.mean(ugrad, 0)
ustar_avg = np.mean(ustar)

pl.rc('font',**{'family':'serif','serif':['Palatino']})
pl.rc('text', usetex=True)

pl.figure()
for n in range(t.size):
    pl.plot(u[n,:]/ustar[n], z/hc, color='#cccccc')
pl.plot(u_avg/ustar_avg, z/hc)
pl.xlabel(r'$u / u_*$')
pl.ylabel(r'$z / h_c$')
pl.xlim(0, 8)
pl.ylim(0, 3)

pl.figure()
for n in range(t.size):
    pl.plot(ustd[n,:]/ustar[n], z/hc, color='#cccccc')
pl.plot(ustd_avg/ustar_avg, z/hc)
pl.xlabel(r'$\sigma_u / u_*$')
pl.ylabel(r'$z / h_c$')
pl.ylim(0, 6)

pl.figure()
for n in range(t.size):
    pl.plot(vstd[n,:]/ustar[n], z/hc, color='#cccccc')
pl.plot(vstd_avg/ustar_avg, z/hc)
pl.xlabel(r'$\sigma_v / u_*$')
pl.ylabel(r'$z / h_c$')
pl.ylim(0, 6)

pl.figure()
for n in range(t.size):
    pl.plot(wstd[n,:]/ustar[n], zh/hc, color='#cccccc')
pl.plot(wstd_avg/ustar_avg, zh/hc)
pl.xlabel(r'$\sigma_w / u_*$')
pl.ylabel(r'$z / h_c$')
pl.ylim(0, 6)

pl.figure()
for n in range(t.size):
    pl.plot(uflux[n,:]/ustar[n]**2, zh/hc, color='#cccccc')
pl.plot(uflux_avg/ustar_avg**2, zh/hc, 'b-' )
pl.plot(uw_avg   /ustar_avg**2, zh/hc, 'b--')
pl.plot(udiff_avg/ustar_avg**2, zh/hc, 'b:' )
pl.xlabel(r'$\overline{u^\prime w^\prime} / u_*^2$')
pl.ylabel(r'$z / h_c$')
pl.ylim(0, 6)

pl.figure()
for n in range(t.size):
    pl.plot(-uw[n,:]*ugrad[n,:]*hc/ustar[n]**3, zh/hc, color='#cccccc')
pl.plot(-uw_avg*ugrad_avg*hc/ustar_avg**3, zh/hc)
pl.xlabel(r'$-\overline{u^\prime w^\prime} \partial u / \partial z h_c / u_*^3$')
pl.ylabel(r'$z / h_c$')
pl.ylim(0, 6)
