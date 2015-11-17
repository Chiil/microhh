import numpy as np

u_val = 6.     # Wind velocity in m s-1.
h_canopy = 10. # Canopy height in m.
delta = 3.     # Region of smoothing in m.

# Get number of vertical levels and size from .ini file
with open('canopy.ini') as f:
    for line in f:
        if (line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if (line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])

dz = zsize / kmax

# Set the arrays 
z  = np.linspace(0.5*dz, zsize-0.5*dz, kmax)
u  = u_val * np.ones(kmax)
Cd = np.zeros(kmax)

# Generate canopy profile with a sudden step at z = h_canopy
#relative_pad = np.where(z < h_canopy, 1., 0.)

# Generate canopy profile with a smooth step at z = h_canopy
relative_pad = np.zeros(kmax)
for k in range(kmax):
    relative_pad[k] = 0.5 - 0.5*np.math.erf(4.*(z[k] - h_canopy) / delta)
relative_pad /= sum(relative_pad)

print('Sum of relative_pad = {0}'.format(sum(relative_pad)))

# Write the data to a file
proffile = open('canopy.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s}\n'.format('z','u','Cd','relative_pad'))
for k in range(kmax):
    proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E}\n'.format(z[k], u[k], Cd[k], relative_pad[k]))
proffile.close()
