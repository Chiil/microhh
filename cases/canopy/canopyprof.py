import numpy

u_val = 6.

# Get number of vertical levels and size from .ini file
with open('canopy.ini') as f:
    for line in f:
        if (line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if (line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / kmax

# set the height
z = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)
u = u_val * numpy.ones(zsize)

# write the data to a file
proffile = open('canopy.prof','w')
proffile.write('{0:^20s} {1:^20s}\n'.format('z','u'))
for k in range(kmax):
    proffile.write('{0:1.14E} {1:1.14E}\n'.format(z[k], u[k]))
proffile.close()
