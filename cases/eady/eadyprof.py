import numpy

# Get number of vertical levels and size from .ini file
with open('eady.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / kmax

N2 = 99999999.

# set the height
z = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)
b = numpy.zeros(numpy.size(z))
ug = numpy.zeros(numpy.size(z))

# linearly stratified profile
for k in range(kmax):
    b[k] = 666. # to be filled in
    ug[k] = 43434.

# write the data to a file
proffile = open('eady.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s}\n'.format('z','b','ug'))
for k in range(kmax):
    proffile.write('{0:1.14E} {1:1.14E} {2:1.14E}\n'.format(z[k], b[k], ug[k]))
proffile.close()
