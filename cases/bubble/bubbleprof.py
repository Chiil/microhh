import numpy

# set the height
kmax  = 384
zsize = 2400.
dz = zsize / kmax

dthetadz = 0.003

z = numpy.linspace(0.5*dz, zsize-0.5*dz, kmax)
thl = numpy.zeros(numpy.size(z))

# linearly stratified profile
for k in range(kmax):
  thl[k] = 300. + dthetadz*z[k]

# write the data to a file
proffile = open('bubble.prof','w')
proffile.write('{0:^20s} {1:^20s}\n'.format('z','thl'))
for k in range(kmax):
  proffile.write('{0:1.14E} {1:1.14E}\n'.format(z[k], thl[k]))
proffile.close()
