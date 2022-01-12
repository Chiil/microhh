import netCDF4 as nc
import numpy as np

dem = np.fromfile('dem.0000000').reshape(1, 64, 128)

"""
data_u = nc.Dataset('u.xy.nc', 'a')
data_u['u'][:, :, :] = np.where(dem > 0, np.nan, data_u['u'][:, :, :])
data_u.close()

data_v = nc.Dataset('v.xy.nc', 'a')
data_v['v'][:, :, :] = np.where(dem > 0, np.nan, data_v['v'][:, :, :])
data_v.close()
"""

for i in range(11):
    var = 's{0}'.format(i)
    data_s = nc.Dataset('{0}.xy.nc'.format(var), 'a')
    data_s[var][:, :, :] = np.where(dem > 0, np.nan, data_s[var][:, :, :])
    data_s.close()
