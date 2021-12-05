import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

float_type = 'f8'

# Get number of vertical levels and size from .ini file
with open('hossein.ini') as f:
    for line in f:
        if (line.split('=')[0]=='itot'):
            itot = int(line.split('=')[1])
        if (line.split('=')[0]=='jtot'):
            jtot = int(line.split('=')[1])
        if (line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if (line.split('=')[0]=='xsize'):
            xsize = float(line.split('=')[1])
        if (line.split('=')[0]=='ysize'):
            ysize = float(line.split('=')[1])
        if (line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

# define the variables
z = np.zeros(kmax)
u = np.zeros(kmax)

# create non-equidistant grid
alpha = 0.967
for k in range(kmax):
  eta  = -1. + 2.*((k+1)-0.5) / kmax
  z[k] = zsize / (2.*alpha) * np.tanh(eta*0.5*(np.log(1.+alpha) - np.log(1.-alpha))) + 0.5*zsize

# create initial parabolic shape
dpdxls = -3.0e-6
visc   =  1.0e-5
for k in range(kmax):
  u[k] = 1./(2.*visc)*dpdxls*(z[k]**2. - zsize*z[k])


# write the data to a file
nc_file = nc.Dataset("hossein_input.nc", mode="w", datamodel="NETCDF4", clobber=True)

nc_file.createDimension("z", kmax)
nc_z  = nc_file.createVariable("z" , float_type, ("z"))

nc_group_init = nc_file.createGroup("init");
nc_u = nc_group_init.createVariable("u", float_type, ("z"))

nc_z[:] = z[:]
nc_u[:] = u[:]

nc_file.close()

ni = 8
nj = 2
h = 0.35
dx = xsize/itot
dy = ysize/jtot

nblock = int(np.rint(h/dx/2))

stepi = itot//ni
stepj = jtot//nj

print("nblock = {0}, step = {1}\n".format(nblock, stepi))

dem = np.zeros((jtot, itot))
for j in range(nj):
    for i in range(ni):
        iw = stepi//2 - nblock + i*stepi
        ie = stepi//2 + nblock + i*stepi
        js = stepj//2 - nblock + j*stepj
        jn = stepj//2 + nblock + j*stepj

        dem[js:jn, iw:ie] = h

dem.tofile('dem.0000000')

sbot = np.zeros((jtot, itot))
soff = (3*nblock)//2
for j in range(0, 1):
    for i in range(0, 2):
        iw = stepi//2 - 1 + i*stepi
        ie = stepi//2 + 1 + i*stepi
        js = stepj//2 - 1 + j*stepj
        jn = stepj//2 + 1 + j*stepj

        sbot[js+soff:jn+soff, iw:ie] = 0.1

for j in range(1, 2):
    for i in range(0, 2):
        iw = stepi//2 - 1 + i*stepi
        ie = stepi//2 + 1 + i*stepi
        js = stepj//2 - 1 + j*stepj
        jn = stepj//2 + 1 + j*stepj

        sbot[js-soff:jn-soff, iw:ie] = 0.1

sbot.tofile('s_sbot.0000000')
sbot.tofile('s_bot.0000000')

plt.figure()
plt.pcolormesh(dem)
plt.colorbar()

plt.figure()
plt.pcolormesh(dem + sbot)
plt.colorbar()
plt.show()
