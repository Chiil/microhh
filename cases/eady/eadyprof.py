import numpy

# Get number of vertical levels and size from .ini file
with open('eady.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / kmax

N2 = 1.
shear = 1e-4

# set the height
z = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)
b = numpy.zeros(numpy.size(z))
ug = numpy.zeros(numpy.size(z))

# linearly stratified profile
for k in range(kmax):
    b[k] = N2*z[k] # to be filled in
    ug[k] = shear*z[k]
u=ug.copy()

# write the data to a file
proffile = open('eady.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s}\n'.format('z','b','u','ug'))
for k in range(kmax):
    proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E}\n'.format(z[k], b[k], u[k], ug[k]))
proffile.close()
