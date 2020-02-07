import statistics as stats
import sys
import numpy as np

filename = sys.argv[1]
load = sys.argv[2] #reads in as string, not True or False
lag = int(sys.argv[3]) 
#print(type(lag))

meani = stats.iter_mean(filename)
print('mean computed iteratively')
covi_spatial = stats.iter_spatial_cov(filename)
print('spatial covariance computed iteratively')
covi_st = stats.iter_spacetime_cov(filename,covi_spatial,lag)
print('spacetime covariance computed iteratively')

if load=='False' or load=='false':
	pass

else:
	mean = stats.full_mean(filename)
	print('full mean computed')
	if np.allclose(mean,meani)==True:
		print('means match up!')
	else:
		print('ALERT: mean calculations do NOT match')

	cov_spatial = stats.full_spatial_cov(filename)
	print('full cov computed')
	if np.allclose(cov_spatial,covi_spatial)==True:
		print('cov matches up!')
	else:
		print('ALERT: cov calculations do NOT match')
		print('full spatial max is',np.max(cov_spatial),'and iterative max is',np.max(covi_spatial))

	cov_st = stats.full_spacetime_cov(filename,cov_spatial,lag)
	print('full spacetime cov computed with lag',lag)
	if np.allclose(cov_st,covi_st,atol=0.001,rtol=0.001)==True:
		print('spacetime covs match up (enough)!')
	else:
		print('ALERT: spacetime covs calculations do NOT match')
		print('full ST max is',np.max(cov_st),'and iterative max is',np.max(covi_st))

import matplotlib.pyplot as plt

#if np.allclose(cov_st,covi_st)==False and load!=False:
	#plt.imshow(cov_st)
	#plt.show()

plt.imshow(cov_st[0,0,:,:])
plt.show()

#print(np.max(cov_spatial),np.min(cov_spatial))

##want to save files
from astropy.io import fits
hdu = fits.PrimaryHDU(covi_st)
hdu.writeto('/Volumes/Backup-Plus/fields_1000_0.001_Vortex_contrast5e3-summed-ST-COV.fits',overwrite=True)
print('full ST max is',np.max(cov_st),'and iterative max is',np.max(covi_st))
print('full ST min is',np.min(cov_st),'and iterative min is',np.min(covi_st))
