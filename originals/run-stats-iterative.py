import statistics as stats
import sys
import numpy as np
from astropy.io import fits

#file info
filename = sys.argv[1]

filename_parts = filename.split("/")
nameonly = filename_parts[-1]
prefix = nameonly.split('.')[0]
savedir = '/Volumes/Backup-Plus/'

#set parameters for region of interest
start = 75
end = 175

##first, let's do the spatial covariances
print('starting spatial covariance')
cov_spatial = stats.iter_spatial_cov(filename,start=start,end=end)
print('spatial cov computed iteratively')

##then, we need to run space-time cov for multiple lags
##need to decide what lags we care about. placeholder here for now
lags = [0,1,2,3,4,5,10,20,30,40,50,70,100]
for lag in lags:
	print('starting covariance for lag',lag)
	if lag==0:
		cov_st = cov_spatial
		covs = cov_st[np.newaxis,:,:]
		print('lag 0 complete')
	else:
		cov_st = stats.iter_spacetime_cov(filename,cov_spatial,lag,start=start,end=end)
		cov_st = cov_st[np.newaxis,:,:]
		covs = np.concatenate([covs,cov_st],axis=0)
		print('lag',lag,'complete')

print('cov array made iteratively with shape',np.shape(covs),'(n_lag,nx*ny,nx*ny)')

#save file with header info
hdr = fits.Header()
hdr['RANGE'] = str([start,end])
hdr['METHOD'] = 'iterative'
hdr['INPUT'] = filename
hdr['OUTPUT'] = str(np.shape(covs))
hdr['LAGS'] = str(lags)
empty_primary = fits.PrimaryHDU(header=hdr)
hdu = fits.ImageHDU(covs)
hdul = fits.HDUList([empty_primary,hdu])
hdul.writeto(savedir+prefix+'-covariances.fits',overwrite=True)


print('file saved at',savedir+prefix+'-covariances.fits')


