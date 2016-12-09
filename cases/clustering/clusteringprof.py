#==============================================================================
# Import necessary modules
#==============================================================================

import numpy
#from scipy.special import erf
from pylab import *

#==============================================================================
# User settings
#==============================================================================

case = 'sessions' # Based on Sessions et al(2014)
case_rad = True
#case_rad_mod = 'dryminwet'

# Wet = 'wet' cooling profile of (Muller, 2014)
# Dry = 'dry' cooling profile of (Muller, 2014)
# dryminwet = 'dry' minus 'wet' cooling profile of (Muller, 2014). This profile only cools in the lowest 2000 m of the model.
#case = 'ss08' # Moist RICO from Stevens/Seifert & Seifert/Heus
#case = 'test' # More moist mixed-layer for testing


#==============================================================================
# Define grid and stretch grid spacing change between two neighbouring grids in the z-direction.
#==============================================================================

# set the height (ktot = 512)
kmax = 192
dn   = 1./kmax

#Create series starting from 1/512 all the way to 1-(1/512) in 512 steps.

n  = numpy.linspace(dn, 1.-dn, kmax)

nloc1 = 24.*dn
nbuf1 = 10.*dn

nloc2 = 192.*dn
nbuf2 = 45.*dn

dz1 = 40 #z0 is calculated as 7.37e-4
dz2 = 80
dz3 = 900

dzdn1 = dz1/dn
dzdn2 = dz2/dn
dzdn3 = dz3/dn

dzdn = dzdn1 + 0.5*(dzdn2-dzdn1)*(1. + numpy.tanh((n-nloc1)/nbuf1)) + 0.5*(dzdn3-dzdn2)*(1. + numpy.tanh((n-nloc2)/nbuf2))

dz = dzdn*dn

z       = numpy.zeros(numpy.size(dz))
stretch = numpy.zeros(numpy.size(dz))

z      [0] = 0.5*dz[0]
stretch[0] = 1.

for k in range(1,kmax):
  z      [k] = z[k-1] + 0.5*(dz[k-1]+dz[k])
  stretch[k] = dz[k]/dz[k-1]

zsize = z[kmax-1] + 0.5*dz[kmax-1]
print('zsize = ', zsize)

b0    = 1.
delta = 4.407731e-3
N2    = 3.

b = numpy.zeros(numpy.size(z))

for k in range(kmax):
  #b[k] = N2*z[k] + b0*erf(-0.5*z[k]/delta) + b0
  b[k] = N2*z[k]

# write the data to a file
proffile = open('rico.prof','w')
proffile.write('{0:^20s} {1:^20s}\n'.format('z','b'))
for k in range(kmax):
  proffile.write('{0:1.14E} {1:1.14E}\n'.format(z[k], b[k]))
proffile.close()

"""
#plot the grid
plt.close('all')
plt.figure()
plt.subplot(131)
plt.plot(n,z,'bo-')
plt.subplot(132)
plt.plot(n,dz,'bo-')
plt.subplot(133)
plt.plot(n,stretch,'bo-')
plt.show()
"""



# Get number of vertical levels and size from .ini file
""" 
with open('rico.ini') as f:
    for line in f:
        if(line.split('=')[0]=='ktot'):
            kmax = int(line.split('=')[1])
        if(line.split('=')[0]=='zsize'):
            zsize = float(line.split('=')[1])
"""
#==============================================================================
# Set variables
#==============================================================================
dz = zsize / kmax

# set the height
#zz     = np.linspace(0.5*dz, zsize-0.5*dz, kmax)
thl     = np.zeros(z.size)
qt      = np.zeros(z.size)
radmean = np.zeros(z.size)
radflex = np.zeros(z.size)
#print(kmax)

print('Setup = %s'%case)

#==============================================================================
# Define profiles
#==============================================================================

for k in range(kmax):

    # Liquid water potential temperature: same in GCSS and SS08
    if(case == 'sessions'):
        if(z[k] < 200.):
            thl[k] = 300
        elif z[k] >= 200 and z[k] < 8000:
            thl[k] = 0.00513 * z[k] + 298.97
        elif z[k] >= 8000 and z[k] < 13000:
            thl[k] = 0.0018 * z[k] + 325.6
        elif z[k] >= 13000 and z[k] < 16000:
            thl[k] = 3.33*10**-4 * z[k] + 344.67
        elif z[k] > 16000:
            thl[k] = 0.0125 * z[k] + 150
    else: 
        if(z[k] < 740.):
            thl[k] = 297.9
        else:
            thl[k] = 297.9 + (317.0 - 297.9)/(4000. - 740.) * (z[k] - 740.) 

    #moisture profile
    if(case == 'gcss'):
        if(z[k] < 740.):
            qt[k] = 16.0 + (13.8 - 16.0) / 740. * z[k]
        elif(z[k] < 3260.):
            qt[k] = 13.8 + (2.4 - 13.8) / (3260. - 740.) * (z[k] - 740.) 
        else:
            qt[k] = 2.4 + (1.8 - 2.4)/(4000. - 3260.) * (z[k] - 3260.) 

    elif(case == 'ss08'):
        if(z[k] < 740.):
            qt[k] = 16.0 + (13.8 - 16.0) / 740. * z[k]
        elif(z[k] < 3260.):
            qt[k] = 13.8 + (4.4 - 13.8) / (3260. - 740.) * (z[k] - 740.) 
        else:
            qt[k] = 4.4 + (3.6 - 4.4)/(4000. - 3260.) * (z[k] - 3260.) 

    elif(case == 'test'):
        q0 = 18.
        q1 = 15.8
        if(z[k] < 740.):
            qt[k] = q0 + (q1 - q0) / 740. * z[k]
        elif(z[k] < 3260.):
            qt[k] = q1 + (2.4 - q1) / (3260. - 740.) * (z[k] - 740.) 
        else:
            qt[k] = 2.4 + (1.8 - 2.4)/(4000. - 3260.) * (z[k] - 3260.) 
            
    elif(case == 'sessions'): 
        if(z[k]<200):
            qt[k] = 0.019 - 1.5*10**-6 * z[k]
        elif(z[k] >= 200 and z[k] < 2500):
            qt[k] = 0.0194 - 3.78*10**-6 * z[k]
        elif(z[k] >= 2500 and z[k] < 7500):
            qt[k] = 0.015 - 2*10**-6 * z[k]
        else: 
            qt[k] = 0
    
    if(case != 'sessions'):
        qt[k] = qt[k]/1000
	
	

    # Advective and radiative tendency thl

    if case_rad == True:
        if z[k] < 7000 :
            radmean[k] = -2 / 86400.
        elif z[k] >= 7000 and z[k] < 12500:
            radmean[k] = (2/5500*z[k]-(2+2/5500*7000)) / 86400
        elif z[k] >= 12500:
            radmean[k] = 0 / 86400

        if z[k]<2000:
            radflex[k] = (-6+0.003*z[k]) / 86400  
        elif z[k]>= 2000:
            radflex[k] = 0.

    else:
        radmean[k] = -0.5 / 86400.



#==============================================================================
# #Smoothen lines at points where profile makes a sudden bend ('knikpunten').
#==============================================================================

#Function definition
def filter_profile(prof, n):
    new_prof = np.copy(prof)
    for i in range(n):
        new_prof[1:-1] = (new_prof[0:-2] + new_prof[2:]) / 2.
    return new_prof

if(case == 'sessions'):

#Temperature 

    hgt = 11200 #defines height above which stronger filtering needs to be applied
    thl_high=np.searchsorted(z[:], hgt)
    thl[thl_high:] = filter_profile(thl[thl_high:], 8)
    thl= filter_profile(thl, 15)
    
    #Moisture
    qt = filter_profile(qt, 15)

    #Radiative forcing
    radmean = filter_profile(radmean, 15)
    radflex = filter_profile(radflex, 15)

#==============================================================================
# Write the data to a file
#==============================================================================
proffile = open('rico.prof','w')
proffile.write('{0:^20s} {1:^20s} {2:^20s} {3:^20s} {4:^20s}\n'.format('z','thl','qt','radmean','radflex'))
for k in range(kmax):
    proffile.write('{0:1.14E} {1:1.14E} {2:1.14E} {3:1.14E} {4:1.14E}\n'.format(z[k], thl[k], qt[k], radmean[k], radflex[k]))
proffile.close()

#==============================================================================
# Define bottom constants
#==============================================================================

ep = 287.04 / 461.5 

# Surface settings
def esat(T):
    c0 = 0.6105851e+03; c1 = 0.4440316e+02; c2 = 0.1430341e+01; c3 = 0.2641412e-01 
    c4 = 0.2995057e-03; c5 = 0.2031998e-05; c6 = 0.6936113e-08; c7 = 0.2564861e-11 
    c8 = -.3704404e-13 
    x  = max(-80.,T-273.15)
    return c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*c8)))))))

def qsat(p, T):
    return ep*esat(T)/(p-(1.-ep)*esat(T))

ps  = 101540.
SST = 303
ths = SST / (ps/1.e5)**(287.04/1005.)
qs  = qsat(ps, SST) 
print('sbot[thl]=%f, sbot[qt]=%f'%(ths, qs))

