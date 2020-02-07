import statistics as stats
import sys
import numpy as np
from astropy.io import fits

filename = sys.argv[1]
prefix = filename.split('.h')[0]
print(prefix)
full = sys.argv[2]

##need to decide what lags we care about. placeholder here for now
lags = [0,1,2,3,4,5]

##only run full version if needed for comparison
if full=='True' or full=='true':
	for lag in lags:
		if lag==0:
			cov_st = stats.full_spacetime_cov(filename,lag)
			covs = cov_st[:,:,:,:,np.newaxis]
		else:
			print(np.shape(covs),np.shape(cov_st))
			cov_st = stats.full_spacetime_cov(filename,lag)
			cov_st = cov_st[:,:,:,:,np.newaxis]
			covs = np.concatenate([covs,cov_st],axis=-1)

	print('cov array made in full with shape',np.shape(covs),'(nx,ny,nx,ny,n_lag)')

	#save file
	hdu = fits.PrimaryHDU(covs)
	hdu.writeto(prefix+'-covariances-full.fits',overwrite=True)

	print('file saved with name',prefix+'-covariances-full.fits')

#always run iterative version!
for lag in lags:
	if lag==0:
		covi_st = stats.iter_spacetime_cov(filename,lag)
		covs = covi_st[:,:,:,:,np.newaxis]
	else:
		covi_st = stats.iter_spacetime_cov(filename,lag)
		covi_st = covi_st[:,:,:,:,np.newaxis]
		covs = np.concatenate([covs,covi_st],axis=-1)

print('cov array made iteratively with shape',np.shape(covs),'(nx,ny,nx,ny,n_lag)')

#save file
hdu = fits.PrimaryHDU(covs)
hdu.writeto(prefix+'-covariances-iterative.fits',overwrite=True)

print('file saved with name',prefix+'-covariances-iterative.fits')



