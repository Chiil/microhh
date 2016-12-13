import numpy as np
import matplotlib.pyplot as pl
import netCDF4 as nc

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

z = np.arange(0., 400., dz)

theta = np.zeros(kmax)
theta = np.where(z<=270, 277.7, theta)
theta = np.where(np.logical_and(z>270,z<320), 277.7+0.018*(z-270), theta)
theta = np.where(z>=320, 278.6+0.0075*(z-320), theta)

ug0 = 1.25
vg0 = 4.5

ug = ug0*np.ones(kmax)
vg = vg0*np.ones(kmax)



# write the data to a file
#proffile = open('drycbl.prof','w')
#proffile.write('{0:^20s} {1:^20s}\n'.format('z','b'))
#for k in range(kmax):
#    proffile.write('{0:1.14E} {1:1.14E}\n'.format(z[k], b[k]))
#proffile.close()

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
pl.plot(u_ref , z_ref, 'k:')
pl.plot(ug, z, 'b-')
pl.plot(ug_ref, z_ref, 'k:')
pl.xlim(0, 6  )
pl.ylim(0, 400)
pl.subplot(122)
pl.plot(v_ref , z_ref, 'k-')
pl.plot(vg, z, 'b-')
pl.plot(vg_ref, z_ref, 'k:')
pl.xlim(0, 6  )
pl.ylim(0, 400)

pl.show()
