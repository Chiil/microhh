import netCDF4 as nc
import numpy as np

dem = np.fromfile('dem.0000000').reshape(1, 128, 256)

data_u = nc.Dataset('u.xy.nc', 'a')
data_u['u'][:, :, :] = np.where(dem > 0, np.nan, data_u['u'][:, :, :])
data_u.close()

data_v = nc.Dataset('v.xy.nc', 'a')
data_v['v'][:, :, :] = np.where(dem > 0, np.nan, data_v['v'][:, :, :])
data_v.close()

data_s = nc.Dataset('s.xy.nc', 'a')
data_s['s'][:, :, :] = np.where(dem > 0, np.nan, data_s['s'][:, :, :])
data_s.close()
