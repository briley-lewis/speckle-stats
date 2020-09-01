import math
from astropy.io import fits
import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
import sys

savedir='/Users/blewis/Desktop/grad-school/Research/Speckle-Stats/Outputs/new-stKLIP/'
datadir='/Users/blewis/Desktop/grad-school/Research/Speckle-Stats/'

##### basic functions
#####no guarantees the legacy version is fully functional!! it was for an older version of MEDIS files, kind of a workaround

def load_full(filename,legacy=False,window=[0,256]):
	"""loads in full file - do NOT use this for large files!"""
	#import file
	f_in = h5py.File(datadir+filename, 'r')

	if legacy==True:
	##FOR OLD FILE FORMAT
	#make list of keys in correct order
		n_screen = len(f_in.keys())
		nums = np.arange(0,n_screen)
		keys = ['t' + str(n) for n in nums]

		#read in file and sum along wavefront axis, and square to get intensity
		intensities=[]
		start = window[0]
		end = window[1]
		for key in keys:
			frame = f_in[key]
			frame = np.sum(frame,axis=2)
			intensity = np.abs(frame[1,0,start:end,start:end])**2
			intensities.append(intensity)

		print('data imported')

	if legacy==False:
		start = window[0]
		end = window[1]
		data = f_in['data']
		frames = np.sum(data,axis=3) #sum along wavefront axis, so that you get both planet and star signal 
		focal_plane = frames[:,1,0,:,:] #show focal plane (change 1 to 0 for pupil)
		intensities = np.abs(focal_plane)**2

		print('data imported')

	else:
		print('data import error')

def iter_mean(f_in,start,end,secondary=False):
	"""get the mean image from an hdf5 dataset, iterative approach"""
	#open file
	if secondary==True:
		f_in = h5py.File(savedir+f_in, 'r')
		data = f_in['data']
		#print(np.shape(data))
		arr_start=1
		shape=np.shape(data)
	if secondary==False:
		f_in = h5py.File(datadir+f_in, 'r')
		data = f_in['data']
		#print(np.shape(data))
		arr_start=0
		shape=np.shape(data[0,0,0,0,start:end,start:end])
	size = shape[-1]
	
	#iterate through to compute mean
	mn_im = np.zeros((size,size)) #initializing array with same spatial dimensions
	k = 0
	for num in np.arange(arr_start,np.shape(data)[0]): ##exclude first frame of zeros from file initialization
		if secondary==False:
			frame = data[num,1,0,:,start:end,start:end]
			frame = np.sum(frame,axis=0)
			intensity = np.abs(frame)**2
		else:
			intensity = data[num,:,:]
		k = k+1
		mn_im = ((k-1)*mn_im + intensity)/k

	f_in.close()

	return mn_im

def cov2d(intensities,lag=0):
	###untested!
	"""2d covariance at a given lag, using fully loaded file. do NOT use for large files! takes in intensities, so you need to use load_full first"""
	steps, nx, ny = np.shape(intensities)
	mean = np.mean(intensities,axis=0)
	mean_subtracted = np.reshape(intensities-mean[np.newaxis,:,:], (steps, ny*nx))
	cov = np.dot(mean_subtracted[0:steps-lag,:].T, mean_subtracted[lag:,:])/(steps-1)
	if lag!=0:
		cov = (cov+cov.T)/2
	return cov

def iter_cov2d(f_in,lag,start,end,legacy=False,return_mean=False,verbose=True):
	if legacy==True:
		#open file
		f_in = h5py.File(datadir+f_in, 'r')
		#print('file opened')
		
		#make list of keys in correct order
		n_screens = np.arange(0,len(f_in.keys()))
		keys = ['t' + str(n) for n in n_screens]
		nx,ny=[end-start,end-start]

		#iterate through to compute mean and covariance
		n = nx*ny
		mn_im = np.zeros((n,)) # initialize
		cov = np.empty((n,n))
		vv = np.empty((n,n))
		intensities = []
		k = 0

		for k in range(lag,len(n_screens)-lag):
			if k>=len(keys) or (k-lag)<0: ##catching for array out of bounds errors
				pass
			else:
				framek = f_in[keys[k]]
				ik = np.abs(framek[1,0,start:end,start:end])**2
				ik = np.reshape(ik,n) #image at k step
				intensities.append(ik)
				framekm1 = f_in[keys[k-lag]]
				ikm1 = np.abs(framekm1[1,0,start:end,start:end])**2
				ikm1 = np.reshape(ikm1,n)
				k = k+1
				######THING TO CONSIDER: what version of mean should we be using? for kth step, or k-lag step? 
				##currently using mean at step k
				##also, maybe our pre-factors should be different k-l instead of k? need to check
				mn_im = ((k-1)*mn_im + ik)/k 
				if k == 1:
					#start with regular covariance
					vk = ik-mn_im
					vkm1 = ikm1-mn_im
					cov = np.outer(vkm1,vk)
				else:
					#for all other timesteps, use iterative calculation
					vk = ik-mn_im
					vkm1 = ikm1-mn_im
					vv = ((k)/(k-1)**2)*np.outer(vkm1,vk)
					cov *= (k-2)/(k-1)
					cov += vv
					if verbose==True:
						print('step',k,'done')
		f_in.close()

	else:
		#open file
		f_in = h5py.File(datadir+f_in, 'r')
		data = f_in['data']
		shape = np.shape(data[:,0,0,0,start:end,start:end])

		#iterate through to compute mean and covariance
		n = (shape[-1])**2
		mn_im = np.zeros((n,)) # initialize
		cov = np.empty((n,n))
		vv = np.empty((n,n))
		k = 0

		for k in range(lag,shape[0]-lag):
			if k>=shape[0] or (k-lag)<0: ##catching for array out of bounds errors
				print('error!')
				pass
			else:
				framek_orig = data[k,1,0,:,start:end,start:end]
				#print(np.shape(framek_orig))
				framek = np.sum(framek_orig,axis=0)
				ik = np.abs(framek)**2
				ik = np.reshape(ik,n) #image at k step
				framekm1 = data[k-lag,1,0,:,start:end,start:end]
				framekm1 = np.sum(framekm1,axis=0)
				ikm1 = np.abs(framekm1)**2
				ikm1 = np.reshape(ikm1,n)
				k = k+1
				######THING TO CONSIDER: what version of mean should we be using? for kth step, or k-lag step? 
				##currently using mean at step k
				##also, maybe our pre-factors should be different k-l instead of k? need to check
				mn_im = ((k-1)*mn_im + ik)/k 
				if k == 1:
					#start with regular covariance
					vk = ik-mn_im
					vkm1 = ikm1-mn_im
					cov = np.outer(vkm1,vk)
				else:
					#for all other timesteps, use iterative calculation
					vk = ik-mn_im
					vkm1 = ikm1-mn_im
					vv = ((k)/(k-1)**2)*np.outer(vkm1,vk)
					cov *= (k-2)/(k-1)
					cov += vv
					if verbose==True:
						print('step',k,'done out of',shape[0]-lag)
		f_in.close()

	#return outputs
	if return_mean:
		return mn_im, cov
	else:
		return cov

def stcov_matrix(img,nlag):
	####untested!
	"""creates larger block diagonal covariance matrix out of smaller single lag covariance matrices"""
	npix = (np.shape(img)[1])**2
	cov = np.full((nlag*npix,nlag*npix),np.nan)
	i = 0
	j = 1
	#start with diagonals
	while j<nlag+1:
		print('covering indices',i*npix,'to',j*npix)
		cov[i*npix:j*npix,i*npix:j*npix] = get_cov2d(img,lag=0)
		i = i+1
		j = j+1
	#then, off diagonals!
	total = npix*nlag
	for k in range(0,nlag):
		j=1
		print('k',k,'completed')
		while ((j+k)*npix)<(npix*nlag):
			print('covering indices',(j-1)*npix,'to',(j+k)*npix,'and',(j+k)*npix,'to',(j+k+1)*npix)
			#print(np.shape(cov[(j-1)*npix:j*npix,(j+k)*npix:(j+k+1)*npix]))
			cov_now = get_cov2d(img,lag=k+1)
			cov[(j-1)*npix:j*npix,(j+k)*npix:(j+k+1)*npix] = cov_now
			cov[(j+k)*npix:(j+k+1)*npix,(j-1)*npix:j*npix] = cov_now
			j = j+1

	return cov

def iter_stcov_matrix(filename,nlag,start,end,legacy=False):
	"""creates larger block diagonal covariance matrix out of smaller single lag covariance matrices"""
	f_in = h5py.File(datadir+filename, 'r')
	data = f_in['data']
	shape = np.shape(data[0,0,0,0,start:end,start:end])
	f_in.close()
	npix = (shape[-1])**2 
	#print('file loaded')
	print(nlag,npix)
	#if ((nlag*npix*32)**2)>5e11:
		#print('WARNING: SIZE TOO BIG CANNOT HOLD IN MEMORY')
	cov = np.full((nlag*npix,nlag*npix),np.nan) ###consider writing to hdf5 file, would enable larger windows and larger nlag. not needed though since eigendecom is bigger bottleneck
	print('cov matrix shell made')
	i = 0
	j = 1
	#start with diagonals
	while j<nlag+1:
		print('covering indices',i*npix,'to',j*npix)
		cov[i*npix:j*npix,i*npix:j*npix] = iter_cov2d(filename,0,start,end,legacy=legacy,return_mean=False,verbose=False)
		i = i+1
		j = j+1
	#then, off diagonals!
	total = npix*nlag
	for k in range(0,nlag):
		j=1
		print('k',k,'completed')
		while ((j+k)*npix)<(npix*nlag):
			print('covering indices',(j-1)*npix,'to',(j+k)*npix,'and',(j+k)*npix,'to',(j+k+1)*npix)
			#print(np.shape(cov[(j-1)*npix:j*npix,(j+k)*npix:(j+k+1)*npix]))
			cov_now = iter_cov2d(filename,(k+1),start,end,legacy=legacy,return_mean=False,verbose=False)
			cov[(j-1)*npix:j*npix,(j+k)*npix:(j+k+1)*npix] = cov_now
			cov[(j+k)*npix:(j+k+1)*npix,(j-1)*npix:j*npix] = cov_now
			j = j+1

	return cov

def iterative_eigendecomp(filename,legacy=False):
	"""considered implemented an iterative algorithm for eigendecomposition / PCA since that is also computationally heavy; considered NIPALS, but did not implement due to need for manually created stcov matrix"""
	print('not implemented - could use NIPALS algorithm?')

def eigendecomp(cov,max_ev):
	"""eigendecomposition using scipy eigh - max_ev sets the maximum number of eigenvalues computed (e.g. max_ev = 100 means we only keep the 100 largest ev)"""
	if max_ev=='all':
		max_ev = np.shape(cov)[0]
	first_ev = np.shape(cov)[0]-1
	last_ev = np.shape(cov)[0]-max_ev
	print('calculating eigenvectors',first_ev,'through',last_ev)
	ev0, P0 = eigh(cov,subset_by_index=[last_ev,first_ev])
	print('eigendecomposition complete')

	#reverse arrays so they're in descending order
	ev0 = ev0[::-1]
	P0 = P0[:,::-1]
	P0 = P0.T
	#P0.shape = max_ev,nl,nx,ny
	#print(np.shape(P0))
	print('rearrangement complete')

	return ev0,P0

def nipals(n_real,n_var):
	"Nonlinear iterative partial least-squares approach to PCA --- NOT USED"
	pass
	"""# following http://folk.uio.no/henninri/pca_module/
	##taken from pyfitz misc.py
	####NOT CURRENTLY USED
	
	n_real, n_var = x.shape

	if skip_mean_subtract:
		xms = x.copy()
	else:	
		# average across realizations
		x_mn = x.mean(axis=0)
		# remove mean
		xms = x-x_mn[n.newaxis,:]

	# number of principal components
	n_pc = n_pc_max if ((n_pc_max is not None) and (n_var >= n_pc_max)) else n_var


	# NIPALS algorithm
	E = xms
	t = xms[:,0].copy()
	P = n.empty((n_var,n_pc), dtype=n.float)
	ev = n.empty(n_pc, dtype=n.float)
	for i in range(n_pc):
		#_log.debug("NIPALS principal component %d/%d" % (i, n_pc))
		tau_old = 0.
		for j in range(tau_iter_max):
			p = n.dot(E.T,t) / n.dot(t.T,t)
			p /= n.sqrt(n.dot(p.T,p))
			t = n.dot(E,p) #/ n.dot(p.T,p)
			tau = n.dot(t.T, t)
			if n.abs(1.-tau_old/tau) < thresh:
				break
			else:
				tau_old = tau
		E -= n.outer(t, p)
		# store
		P[:,i] = p
		ev[i] = tau
		# ev thresh check
		if ev_thresh is not None:
			if i==0: ev0 = tau
			elif tau/ev0 < ev_thresh: break
	if i != n_pc-1: # result of ev_thresh
		P = P[:,0:i+1] #eigenvectors
		ev = ev[0:i+1] #eigenvalues

	return ev, P""",

def stKLIP(ev0,P0,f_in,num_ev=10,seq_len=5,window=[0,256],iterative=True,return_all=False,**kwargs):
	"""given a set of eigenvalues and eigenvectors, this runs the rest of KLIP over the *whole* image sequence (not just one image subsequence.
	returns the averaged residuals, and also the whole set of subtracted images and model images. IF USING ITERATIVE, PLEASE SUPPLY MEAN FROM ITERATIVE COV CALC AS KWARG 'mean_img'"""
	start = window[0]
	end = window[1]
	mean_img = kwargs['mean_img']

	if iterative==False:
		filename = f_in.split('.')[0]
		f_in = h5py.File(datadir+f_in, 'r')
		print('file opened for stKLIP')
		data = f_in['data']
		full_seq = data[:,1,0,:,start:end,start:end]
		full_seq=np.sum(full_seq,axis=1)
		full_seq = np.abs(full_seq)**2
		mean_img = np.mean(full_seq,axis=0) ###won't work for iterative
		central_index = int(np.median(np.arange(0,seq_len)))
		#print('for subsequence length',seq_len,'central image has index',central_index)
		
		models=[]
		subtracteds=[]
		i=0
		
		while (i+seq_len-1)<np.shape(full_seq)[0]:
			
			target_seq = full_seq[i:i+seq_len,:,:]

			#choose target image - we're going with one, the "central" image (2 in a series of 5 lags)
			ms_target_seq = target_seq - mean_img
			seq_len,img_shape,img_shape2 = np.shape(target_seq)

			#choose modes
			chosen_ev = np.reshape(P0[0:num_ev],[num_ev,seq_len*img_shape*img_shape]) 
			chosen_evals = ev0[0:num_ev]

			#get image coefficients
			coeffs = np.dot(chosen_ev,np.reshape(target_seq,seq_len*img_shape**2)) 

			#create model
			psf_model = np.reshape(np.dot(coeffs,chosen_ev),[seq_len,img_shape,img_shape]) 
			central_img = target_seq[central_index]
			central_model = psf_model[central_index]
			models.append(central_model)
			subtracted = (central_img-central_model)
			subtracteds.append(subtracted)
			#print('central image',(i+central_index),'done')
			
			i=i+1
			
		averaged = np.average(subtracteds,axis=0)

		print('nothing saved on disk!')
		if return_all==True:
			return averaged, subtracteds, models
		else:
			return averaged

	if iterative==True:
		
		filename = f_in.split('.')[0]
		f_in = h5py.File(datadir+f_in, 'r')
		print('file opened for stKLIP')
		full_seq = f_in['data']

		central_index = int(np.median(np.arange(0,seq_len)))
		#print('for subsequence length',seq_len,'central image has index',central_index)
		
		##initialize save files
		test_seq = full_seq[0,1,0,0,start:end,start:end]
		models = h5py.File(savedir+'{}_models-ev{}-seq{}.h5'.format(filename,num_ev,seq_len), 'w')
		subtracteds = h5py.File(savedir+'{}_subtracteds-ev{}-seq{}.h5'.format(filename,num_ev,seq_len), 'w')
		models.create_dataset('data', data=np.zeros((0,np.shape(test_seq)[0],np.shape(test_seq)[1])), compression="gzip", chunks=True,maxshape=(None, None, None))
		models.close()
		subtracteds.create_dataset('data', data=np.zeros((0,np.shape(test_seq)[0],np.shape(test_seq)[1])), compression="gzip", chunks=True,maxshape=(None, None, None))
		subtracteds.close()

		#reopen for adding on
		models = h5py.File(savedir+'{}_models-ev{}-seq{}.h5'.format(filename,num_ev,seq_len), 'a')
		subtracteds = h5py.File(savedir+'{}_subtracteds-ev{}-seq{}.h5'.format(filename,num_ev,seq_len), 'a')

		#add attributes from original hdf5 simulation
		att_count = 0
		for key in f_in.attrs:
			models.attrs[key] = f_in.attrs[key]
			subtracteds.attrs[key] = f_in.attrs[key]
			att_count = att_count+1
		print(att_count,'attributes written to file')

		i=0
		
		while (i+seq_len-1)<np.shape(full_seq)[0]:
			
			#print(np.shape(full_seq))
			target_seq = full_seq[i:i+seq_len,1,0,:,start:end,start:end]
			target_seq = np.sum(target_seq,axis=1)
			target_seq = np.abs(target_seq)**2
			#print(np.shape(target_seq))
			#print(np.shape(mean_img))

			#choose target image - we're going with one, the "central" image (2 in a series of 5 lags)
			ms_target_seq = target_seq - mean_img
			seq_len,img_shape,img_shape2 = np.shape(target_seq)

			#choose modes
			chosen_ev = np.reshape(P0[0:num_ev],[num_ev,seq_len*img_shape*img_shape]) 
			chosen_evals = ev0[0:num_ev]

			#get image coefficients
			coeffs = np.dot(chosen_ev,np.reshape(target_seq,seq_len*img_shape**2)) 

			#create model
			psf_model = np.reshape(np.dot(coeffs,chosen_ev),[seq_len,img_shape,img_shape]) 
			central_img = target_seq[central_index]
			central_model = psf_model[central_index]

			####okay, this will get out of hand in memory VERY fast! need to write out iteratively to hdf5
			index = models["data"].shape[0]-1
			models["data"].resize((models["data"].shape[0] + 1), axis = 0)
			models["data"][index,:,:] = central_model
			subtracted = (central_img-central_model)
			subtracteds["data"].resize((subtracteds["data"].shape[0] + 1), axis = 0)
			subtracteds["data"][index,:,:] = subtracted
			print('central image',(i+central_index),'done')
			
			i=i+1
	   
		print('models saved at {}_models-ev{}-seq{}.h5'.format(filename,num_ev,seq_len))
		models.close()
		print('subtracteds saved at {}_subtracteds-ev{}-seq{}.h5'.format(filename,num_ev,seq_len))
		subtracteds.close()
		f_in.close()

		print('stklip done, beginning averaging process')
		##compute iterative mean on subtracteds
		averaged = iter_mean('{}_subtracteds-ev{}-seq{}.h5'.format(filename,num_ev,seq_len),start,end,secondary=True) ##start and end should be global variables when running this~ not the best implementation?
		print('average residuals calculated -- not yet saved on disk')

		return averaged




##### wrappers / higher level functions

def run_stKLIP(filename,nlag,nmode,window=[0,256],legacy=False):
	"runs stKLIP once for one mode one lag -- not tested, use model_grid instead"
	name = filename.split('.')[0]
	f_in = h5py.File(datadir+filename,'r')
	start = window[0]
	end = window[1]

	print('creating stcov matrix for {} lags'.format(nlag))
	stcov_matrix = iter_stcov_matrix(filename,nlag,start,end)

	print('doing eigendecomposition')
	ev0, P0 = eigendecomp(stcov_matrix,max_ev=(nmode+1))

	print('computing overall image mean')
	mean = iter_mean(filename,start,end)

	print('data residuals will be saved to {}_avg-res_{}-lags_{}-modes.fits'.format(name,nlag,nmode))
	print('running stKLIP for {} modes'.format(nmode))
	avg_res = stKLIP(ev0,P0,filename,num_ev=nmode,seq_len=nlag,mean_img=mean)
	print('stKLIP completed for {} modes'.format(nmode))

	###SAVE TO FILE
	modes = h5py.File(savedir+'{}_avg-res_{}-lags.h5'.format(name,nlag), 'w')
	modes.create_dataset('data', avg_res, compression="gzip")
	list_modes=nmode.copy()
	list_modes.insert(0,np.nan)
	modes.attrs['modes']= list_modes
	modes.attrs['lag']=nlag
	if window==[0,256]:
	  	modes.attrs['window']='full simulation'
	else:
	   	modes.attrs['window']=window
	modes.attrs['input_file']=filename
	modes.attrs['method']='Load Full'
	att_count = 5
	for key in f_in.attrs:
		modes.attrs[key] = f_in.attrs[key]
		att_count = att_count+1
	print(att_count,'attributes written to file')
	print('LAG {} FINISHED - averaged residuals saved to file at {}_avg-res_{}-lags_single-mode.h5'.format(name,nlag))

	#save to fits cube for easy viewing
	attributes = modes.attrs
	hdr = fits.Header(attributes)
	fits.writeto(savedir+'{}_avg-res_{}-lags_single-mode.fits'.format(name,nlag),avg_res,header=hdr,overwrite=True)
	print('also saved in fits format as {}_avg-res_{}-lags.fits'.format(name,nlag))
	f_in.close()
	modes.close()


def model_grid(filename,nlags,nmodes,window=[0,256],legacy=False):
	"""future function / wrapper for testing over multiple lags + modes. nlags is an array of lags to test, nmodes is a number of modes to test at each lag"""
	name = filename.split('.')[0]
	f_in = h5py.File(datadir+filename,'r')
	start = window[0]
	end = window[1]

	#running over multiple lags, multiple modes
	for nlag in nlags:
		print('creating stcov matrix for {} lags'.format(nlag))
		stcov_matrix = iter_stcov_matrix(filename,nlag,start,end)

		print('doing eigendecomposition')
		ev0, P0 = eigendecomp(stcov_matrix,max_ev=(nmodes[-1]+5)) ##weird nmodes thing here ensures that there are more ev calculated than modes we want to use

		print('computing overall image mean')
		mean = iter_mean(filename,start,end,secondary=False)
		shape = np.shape(mean)


		##initialize save file
		modes = h5py.File(savedir+'{}_avg-res_{}-lags.h5'.format(name,nlag), 'w')
		modes.create_dataset('data', data=np.zeros((0,shape[0],shape[1])), compression="gzip", chunks=True,maxshape=(None, None, None))
		list_modes=nmodes.copy()
		list_modes.insert(0,np.nan)
		modes.attrs['modes']= list_modes
		modes.attrs['lag']=nlag
		if window==[0,256]:
			modes.attrs['window']='full simulation'
		else:
			modes.attrs['window']=window
		#modes.attrs['input_file']=filename
		modes.attrs['method']='Iterative'
		att_count = 5
		for key in f_in.attrs:
			modes.attrs[key] = f_in.attrs[key]
			att_count = att_count+1
		print(att_count,'attributes written to file')
		modes.close()
		print('data residuals will be saved to {}_avg-res_{}-lags.h5'.format(name,nlag))
		#reopen for addition
		modes = h5py.File(savedir+'{}_avg-res_{}-lags.h5'.format(name,nlag), 'a')
		averaged = []

		for mode in nmodes:
			if mode=='nan' or np.isnan(mode)==True or mode=='np.nan':
				pass
			else:
				print('running stKLIP for {} modes'.format(mode))
				avg_res = stKLIP(ev0,P0,filename,num_ev=mode,seq_len=nlag,mean_img=mean,window=[start,end])
				print('stKLIP completed for {} modes'.format(mode))

				###add slice to file; somehow include KL mode in header?
				index = modes["data"].shape[0]-1
				modes["data"].resize((modes["data"].shape[0] + 1), axis = 0)
				modes["data"][index,:,:] = avg_res
				print('slice for {} modes saved to disk'.format(mode))
				averaged.append(avg_res)

		print('LAG {} FINISHED - averaged residuals all saved to file at {}_avg-res_{}-lags.h5'.format(nlag,name,nlag))

		#save to fits cube for easy viewing 
		hdr = fits.Header()
		for key in modes.attrs:
			hdr[key] = str(modes.attrs[key])
			#print(key,str(modes.attrs[key]))
		averaged=np.asarray(averaged)
		fits.writeto(savedir+'{}_avg-res_{}-lags.fits'.format(name,nlag),averaged,header=hdr,overwrite=True)
		print('also saved in fits format as {}_avg-res_{}-lags.fits'.format(name,nlag))
		modes.close()

	f_in.close()

	print('run complete')



###later analysis functions!

def SNR():
	"""computes signal to noise ratio"""
	pass

def contrast_curve():
	"""computes contrast curve for a given residual"""
	pass


###for running from command line

def stringtolist(string):
	splits = string.split(',')
	#print(splits)
	values=[]
	for i in splits:
		values.append(int(i))
	return values

if __name__ == '__main__':
	input_file = sys.argv[1]
	print('running stKLIP for',datadir+input_file)
	print('outputs will write to',savedir)
	lag_string = sys.argv[2]
	input_lags = stringtolist(lag_string)
	print("lags:",input_lags)
	#print(type(input_lags),type(input_lags[0]))
	mode_string = sys.argv[3]
	input_modes = stringtolist(mode_string)
	print("modes:",input_modes)
	window_string = sys.argv[4]
	input_window = stringtolist(window_string)
	print("window:",input_window)
	model_grid(input_file,input_lags,input_modes,window=input_window)



