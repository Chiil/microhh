import numpy as np
import matplotlib.pyplot as pl
import netCDF4 as nc

def smooth_prof(prof, n):
    for i in range(n):
        prof[1:-1] = 0.5*(prof[0:-2] + prof[2::])

    return prof

ncfile = nc.Dataset("gabls4_les.nc","r")

z_ref     = ncfile.variables["height"][:]
theta_ref = ncfile.variables["theta"][:]
u_ref     = ncfile.variables["u"][:]
v_ref     = ncfile.variables["v"][:]
ug_ref    = ncfile.variables["Ug"][0,:]
vg_ref    = ncfile.variables["Vg"][0,:]

# set the height (ktot = 512)
kmax = 512
dn   = 1./kmax

"""
n = np.linspace(dn, 1.-dn, kmax)

nloc1 = 80.*dn
nbuf1 = 16.*dn

nloc2 = 512.*dn
nbuf2 = 72.*dn

dz1 = 0.001
dz2 = 0.002
dz3 = 0.016

dzdn1 = dz1/dn
dzdn2 = dz2/dn
dzdn3 = dz3/dn

dzdn = dzdn1 \
     + 0.5*(dzdn2-dzdn1)*(1. + np.tanh((n-nloc1)/nbuf1)) \
     + 0.5*(dzdn3-dzdn2)*(1. + np.tanh((n-nloc2)/nbuf2))

dz = dzdn*dn

z       = np.zeros(np.size(dz))
stretch = np.zeros(np.size(dz))

z      [0] = 0.5*dz[0]
stretch[0] = 1.

for k in range(1,kmax):
  z      [k] = z[k-1] + 0.5*(dz[k-1]+dz[k])
  stretch[k] = dz[k]/dz[k-1]

zsize = z[kmax-1] + 0.5*dz[kmax-1]
print('zsize = ', zsize)
"""

zsize = 400.
dz = zsize / kmax

z = np.arange(dz/2, 400., dz)

theta = np.zeros(kmax)
theta = np.where(z<=270, 277.7, theta)
theta = np.where(np.logical_and(z>270,z<320), 277.7+0.018*(z-270), theta)
theta = np.where(z>=320, 278.6+0.0075*(z-320), theta)

z_theta0 = np.array([0, 270, 320, 400])
theta0 = np.array([277.7, 277.7, 278.6, 279.2])
theta = np.interp(z, z_theta0, theta0)

ug0 = 1.25
vg0 = 4.5

ug = ug0*np.ones(kmax)
vg = vg0*np.ones(kmax)

z_u0 = np.array([0, 270, 320, 400])
u0 = np.array([2.25, 2.25, 1.25, 1.25])
u = np.interp(z, z_u0, u0)

v0 = np.array([5., 5, 4.5, 4.5])
v = np.interp(z, z_u0, v0)

nsmooth = 128
theta = smooth_prof(theta, nsmooth)
u = smooth_prof(u, nsmooth)
v = smooth_prof(v, nsmooth)

b = 9.81/theta[0] * (theta - theta[0])

# write the data to a file
proffile = open('gabls4_dns.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s} {4:^20s} {5:^20s}\n'.format('z','b','u','v','ug','vg'))
for k in range(kmax):
    proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E} {4:1.14E} {5:1.14E}\n'.format(z[k], b[k], u[k], v[k], ug[k], vg[k]))
proffile.close()

#plot the grid
#pl.figure()
#pl.subplot(131)
#pl.plot(n,z)
#pl.subplot(132)
#pl.plot(n,dz)
#pl.subplot(133)
#pl.plot(n,stretch)

ylim = 400

pl.figure()
pl.plot(theta, z, 'b-')
pl.plot(theta_ref, z_ref, 'k:')
pl.xlim(276, 280)
pl.ylim(0  , 400)

pl.figure()
pl.subplot(121)
pl.plot(u, z, 'b-')
pl.plot(u_ref , z_ref, 'k:')
pl.plot(ug, z, 'r-')
pl.plot(ug_ref, z_ref, 'k:')
pl.xlim(0, 6  )
pl.ylim(0, 400)
pl.subplot(122)
pl.plot(v, z, 'b-')
pl.plot(v_ref , z_ref, 'k:')
pl.plot(vg, z, 'r-')
pl.plot(vg_ref, z_ref, 'k:')
pl.xlim(0, 6  )
pl.ylim(0, 400)

pl.show()
