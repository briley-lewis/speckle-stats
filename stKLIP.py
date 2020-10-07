import math
from astropy.io import fits
import h5py
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
import sys
from photutils import CircularAperture
from photutils import aperture_photometry
sys.path.insert(0,'/u/scratch/b/blewis34/medis_data/observations/complex-fields/')
sys.path.insert(0,'/u/home/b/blewis34/speckle-stats/')

#savedir='/u/scratch/b/blewis34/stKLIP/'
#datadir='/u/scratch/b/blewis34/medis_data/observations/complex-fields/'
datadir='../Outputs/'
savedir='../Outputs/'

##### basic functions

def load_full(filename,window=[0,256,0,256]):
	"""loads in full file - do NOT use this for large files!"""
	#import file
	f_in = h5py.File(datadir+filename, 'r')
	starty = window[0]
	endy = window[1]
	startx = window[2]
	endx = window[3]
	data = f_in['data'] 
	focal_plane = frames[:,1,0,:,starty:endy,startx:endx] #show focal plane (change 1 to 0 for pupil)
	intensities = np.abs(focal_plane)**2
	intensities = np.sum(intensities,axis=1) #sum along wavefront axis, so that you get both planet and star signal 

	print('data imported')
	return intensities

def iter_mean(f_in,starty,endy,startx,endx,secondary=False):
	"""get the mean image from an hdf5 dataset, iterative approach"""

	#open file
	if secondary==True:
		f_in = h5py.File(savedir+f_in, 'r')
		data = f_in['data']
		arr_start=1
		shape=np.shape(data)
	if secondary==False:
		f_in = h5py.File(datadir+f_in, 'r')
		data = f_in['data']
		#print(np.shape(data))
		arr_start=0
		shape=np.shape(data[0,0,0,0,starty:endy,startx:endx])
	size = shape[-1]
	
	#iterate through to compute mean
	mn_im = np.zeros((size,size)) #initializing array with same spatial dimensions
	k = 0
	for num in np.arange(arr_start,np.shape(data)[0]): ##exclude first frame of zeros from file initialization
		if secondary==False:
			frame = data[num,1,0,:,starty:endy,startx:endx]
			intensity = np.abs(frame)**2
			intensity = np.sum(intensity,axis=0)
		else:
			intensity = data[num,:,:]
		k = k+1
		mn_im = ((k-1)*mn_im + intensity)/k

	f_in.close()

	return mn_im

def get_cov2d(intensities,lag=0):
	"""2d covariance at a given lag, using fully loaded file. do NOT use for large files! takes in intensities, so you need to use load_full first"""
	steps, nx, ny = np.shape(intensities)
	mean = np.mean(intensities,axis=0)
	mean_subtracted = np.reshape(intensities-mean[np.newaxis,:,:], (steps, ny*nx))
	cov = np.dot(mean_subtracted[0:steps-lag,:].T, mean_subtracted[lag:,:])/(steps-1)
	if lag!=0:
		cov = (cov+cov.T)/2
	return cov

def iter_cov2d(f_in,lag,starty,endy,startx,endx,return_mean=False,verbose=True):
	#open file
	f_in = h5py.File(datadir+f_in, 'r')
	data = f_in['data']
	shape = np.shape(data[:,0,0,0,starty:endy,startx:endx])

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
			framek = data[k,1,0,:,starty:endy,startx:endx]
			#print(np.shape(framek_orig))
			ik = np.abs(framek)**2
			ik = np.sum(ik,axis=0)
			ik = np.reshape(ik,n) #image at k step
			framekm1 = data[k-lag,1,0,:,starty:endy,startx:endx]
			ikm1 = np.abs(framekm1)**2
			ikm1 = np.sum(ikm1,axis=0)
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

def iter_stcov_matrix(filename,nlag,starty,endy,startx,endx):
	"""creates larger block diagonal covariance matrix out of smaller single lag covariance matrices"""
	f_in = h5py.File(datadir+filename, 'r')
	data = f_in['data']
	shape = np.shape(data[0,0,0,0,starty:endy,startx:endx])
	f_in.close()
	npix = (shape[-1])**2 
	cov = np.full((nlag*npix,nlag*npix),np.nan) ###consider writing to hdf5 file, would enable larger windows and larger nlag. only worthwhile if using dask for eigh though
	print('cov matrix shell made')
	i = 0
	j = 1
	#start with diagonals
	while j<nlag+1:
		print('covering indices',i*npix,'to',j*npix)
		cov[i*npix:j*npix,i*npix:j*npix] = iter_cov2d(filename,0,starty,endy,startx,endx,return_mean=False,verbose=False)
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
			cov_now = iter_cov2d(filename,(k+1),starty,endy,startx,endx,return_mean=False,verbose=False)
			cov[(j-1)*npix:j*npix,(j+k)*npix:(j+k+1)*npix] = cov_now
			cov[(j+k)*npix:(j+k+1)*npix,(j-1)*npix:j*npix] = cov_now
			j = j+1

	return cov

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

def stKLIP(ev0,P0,f_in,num_ev=10,seq_len=5,window=[0,256,0,256],iterative=True,return_all=False,**kwargs):
	"""given a set of eigenvalues and eigenvectors, this runs the rest of KLIP over the *whole* image sequence (not just one image subsequence.
	returns the averaged residuals, and also the whole set of subtracted images and model images. IF USING ITERATIVE, PLEASE SUPPLY MEAN FROM ITERATIVE COV CALC AS KWARG 'mean_img'"""
	starty = window[0]
	endy = window[1]
	startx = window[2]
	endx = window[3]

	if iterative==False:
		filename = f_in.split('.')[0]
		f_in = h5py.File(datadir+f_in, 'r')
		print('file opened for stKLIP')
		data = f_in['data']
		full_seq = data[:,1,0,:,starty:endy,startx:endx]
		full_seq = np.abs(full_seq)**2
		full_seq = np.sum(full_seq,axis=1)
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
		mean_img = kwargs['mean_img']

		central_index = int(np.median(np.arange(0,seq_len)))
		#print('for subsequence length',seq_len,'central image has index',central_index)
		
		##initialize save files
		test_seq = full_seq[0,1,0,0,starty:endy,startx:endx]
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
			target_seq = full_seq[i:i+seq_len,1,0,:,starty:endy,startx:endx]
			target_seq = np.abs(target_seq)**2
			target_seq = np.sum(target_seq,axis=1)
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
			#print('central image',(i+central_index),'done')
			
			i=i+1
	   
		print('models saved at {}_models-ev{}-seq{}.h5'.format(filename,num_ev,seq_len))
		models.close()
		print('subtracteds saved at {}_subtracteds-ev{}-seq{}.h5'.format(filename,num_ev,seq_len))
		subtracteds.close()
		f_in.close()

		print('stklip done, beginning averaging process')
		##compute iterative mean on subtracteds
		averaged = iter_mean('{}_subtracteds-ev{}-seq{}.h5'.format(filename,num_ev,seq_len),starty,endy,startx,endx,secondary=True) ##start and end should be global variables when running this~ not the best implementation?
		print('average residuals calculated -- not yet saved on disk')

		return averaged




##### wrappers / higher level functions

def model_grid(filename,nlags,nmodes,window=[0,256,0,256],legacy=False):
	"""future function / wrapper for testing over multiple lags + modes. nlags is an array of lags to test, nmodes is a number of modes to test at each lag"""
	name = filename.split('.')[0]
	f_in = h5py.File(datadir+filename,'r')
	starty = window[0]
	endy = window[1]
	startx = window[2]
	endx = window[3]

	#running over multiple lags, multiple modes
	for nlag in nlags:
		print('creating stcov matrix for {} lags'.format(nlag))
		stcov_matrix = iter_stcov_matrix(filename,nlag,starty,endy,startx,endx)

		print('doing eigendecomposition')
		ev0, P0 = eigendecomp(stcov_matrix,max_ev=(nmodes[-1]+5)) ##weird nmodes thing here ensures that there are more ev calculated than modes we want to use

		print('computing overall image mean')
		mean = iter_mean(filename,starty,endy,startx,endx,secondary=False)
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
				avg_res = stKLIP(ev0,P0,filename,num_ev=mode,seq_len=nlag,mean_img=mean,window=[starty,endy,startx,endx])
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
		hdr['modes']= str(list_modes)
		hdr['lag']=str(nlag)
		if window==[0,256]:
			hdr['window']='full simulation'
		else:
			hdr['window']=str(window)
		#modes.attrs['input_file']=filename
		hdr['method']='Iterative'
		#for key in modes.attrs:
			#hdr[key] = str(modes.attrs[key])
			#print(key,str(modes.attrs[key]))
		averaged=np.asarray(averaged)
		fits.writeto(savedir+'{}_avg-res_{}-lags.fits'.format(name,nlag),averaged,header=hdr,overwrite=True)
		print('also saved in fits format as {}_avg-res_{}-lags.fits'.format(name,nlag))
		modes.close()

	f_in.close()

	print('run complete')



###later analysis functions!

def SNR(data,vx,vy,w,window,display=True):

	tol=w*3
	
	r = np.zeros((256,256))
	for i in np.arange(0,256):
		for j in np.arange(0,256):
			r[i,j] = np.sqrt((i-128)**2 + (j-128)**2)

	##planet and image parameters
	px = window[2]+vx
	py = window[0]+vy
	rp = np.sqrt((px-128)**2 + (py-128)**2)
	#print(px,py,rp)

	apertures = []
	
	for i in np.arange(0,256):
		for j in np.arange(0,256):
			#print(r[i,j])
			if r[i,j]>rp-0.5 and r[i,j]<rp+0.5:
				if i<px+tol and i>px-tol and j<py+tol and j>py-tol:
					pass
				else:
					if i-w<window[2] or i+w>window[3] or j-w<window[0] or j+w>window[1]:
						pass
					else:
						planet_mark = plt.Circle((i-window[2],j-window[0]), w,facecolor='None',edgecolor='red')
						ax = plt.gca()
						ax.add_artist(planet_mark)
						aper = CircularAperture((i-window[2],j-window[0]), w)
						flux = aperture_photometry(data,aper)['aperture_sum']
						apertures.append(flux)
	if display==True:
		plt.imshow(data,vmax=0.1)
		planet_mark = plt.Circle((vx,vy), w,facecolor='None',edgecolor='yellow')
		ax = plt.gca()
		ax.add_artist(planet_mark)
		plt.title('SNR Apertures')
		plt.colorbar()
		plt.show()
	
	noise = np.asarray(np.average(apertures))
	pa = CircularAperture((vx,vy), w)
	signal = np.asarray(aperture_photometry(data,pa)['aperture_sum'])
	SNR = signal/noise
	print('SNR is',SNR)
	
	return np.asarray([noise,signal,SNR])

def contrast_curve(data,window,w):
	i=15
	contrasts=[]
	pixel_val=[]
	centerx = 128-window[2]
	centery = 128-window[0]
	while centerx+i+5<np.shape(data)[0]:
		ap = CircularAperture((centerx+i,centery), w)
		fl = aperture_photometry(data,ap)['aperture_sum']
		planet_mark = plt.Circle((centerx+i,centery), w,facecolor='None',edgecolor='red')
		ax = plt.gca()
		ax.add_artist(planet_mark)
		mn = fl/ap.area
		contrast = mn/mean_flux
		contrasts.append(contrast)
		pixel_val.append(i)
		i=i+5
	return pixel_val,contrasts


###for running from command line

def stringtolist(string):
	splits = string.split(',')
	#print(splits)
	values=[]
	for i in splits:
		values.append(int(i))
	return values

if __name__ == '__main__':
	print('starting now!')
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
	print("windows (y,y,x,x):",input_window)
	print('moving on to model grid function')
	model_grid(input_file,input_lags,input_modes,window=input_window)



####################deprecated functions
def stcov_matrix(img,nlag):
	####untested!
	"""creates larger block diagonal covariance matrix out of smaller single lag covariance matrices"""
	shape = np.shape(img)
	npix = (shape[1])**2
	cov = np.full((nlag*npix,nlag*npix),np.nan)
	i = 0
	j = 1
	#start with diagonals
	while j<nlag+1:
		#print('covering indices',i*npix,'to',j*npix)
		cov[i*npix:j*npix,i*npix:j*npix] = get_cov2d(img,lag=0)
		i = i+1
		j = j+1
	#then, off diagonals!
	total = npix*nlag
	for k in range(0,nlag):
		j=1
		#print('k',k,'completed')
		while ((j+k)*npix)<(npix*nlag):
			#print('covering indices',(j-1)*npix,'to',(j+k)*npix,'and',(j+k)*npix,'to',(j+k+1)*npix)
			#print(np.shape(cov[(j-1)*npix:j*npix,(j+k)*npix:(j+k+1)*npix]))
			cov_now = get_cov2d(img,lag=k+1)
			cov[(j-1)*npix:j*npix,(j+k)*npix:(j+k+1)*npix] = cov_now
			cov[(j+k)*npix:(j+k+1)*npix,(j-1)*npix:j*npix] = cov_now
			j = j+1

	return cov

def iter_stcov_matrix2(filename,nlag,starty,endy,startx,endx,legacy=False):
	"""creates larger block diagonal covariance matrix out of smaller single lag covariance matrices. this one saves it to an h5 file so that it doesn't overwhelm memory, but takes longer!
	note: will still cause memory errors unless eigendecomp is changed, since you've gotta load in the whole matrix for that. also, this still isn't fully functional so DONT USE IT OK.
	as of 9/3/2020 this for some reason fails np.allclose with the output of iter_stcov_matrix function."""
	print('starting stcov matrix!')
	f_in = h5py.File(datadir+filename, 'r')
	splits = filename.split('.')
	name = splits[0]
	data = f_in['data']
	shape = np.shape(data[0,0,0,0,starty:endy,startx:endx])
	f_in.close()
	print('file loaded')
	npix = (shape[-1])**2
	hf = h5py.File('{}_cov-matrix_{}-lags.h5'.format(name,nlag),'w')
	hf.create_dataset('cov',(nlag*npix,nlag*npix))
	hf.close()
	hf = h5py.File('{}_cov-matrix_{}-lags.h5'.format(name,nlag), 'a')
	print('saving cov matrix at {}_cov-matrix_{}-lags.h5'.format(name,nlag))
	i = 0
	j = 1
	#start with diagonals
	while j<nlag+1:
		print('covering indices',i*npix,'to',j*npix)
		hf['cov'][i*npix:j*npix,i*npix:j*npix] = iter_cov2d(filename,0,starty,endy,startx,endx,legacy=legacy,return_mean=False,verbose=False)
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
			cov_now = iter_cov2d(filename,(k+1),starty,endy,startx,endx,legacy=legacy,return_mean=False,verbose=False)
			hf['cov'][(j-1)*npix:j*npix,(j+k)*npix:(j+k+1)*npix] = cov_now
			hf['cov'][(j+k)*npix:(j+k+1)*npix,(j-1)*npix:j*npix] = cov_now
			j = j+1

	#return np.asarray(hf['cov']) ###THIS WILL CAUSE MEMORY ERRORS. REMOVE IF IMPLEMENTING ITERATIVE EIGNENDECOMP
	print('matrix created at {}_cov-matrix_{}-lags.h5, closing file'.format(name,nlag))
	hf.close()

def iterative_eigendecomp(filename,legacy=False):
	"""considered implemented an iterative algorithm for eigendecomposition / PCA since that is also computationally heavy; considered NIPALS, but did not implement due to need for manually created stcov matrix"""
	print('not implemented - could use NIPALS algorithm?')

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

def run_stKLIP(filename,nlag,nmode,window=[0,256],legacy=False):
	"runs stKLIP once for one mode one lag -- not tested, use model_grid instead"
	name = filename.split('.')[0]
	f_in = h5py.File(datadir+filename,'r')
	start = window[0]
	end = window[1]

	print('creating stcov matrix for {} lags'.format(nlag))
	stcov_matrix = iter_stcov_matrix(filename,nlag,starty,endy,startx,endx)

	print('doing eigendecomposition')
	ev0, P0 = eigendecomp(stcov_matrix,max_ev=(nmode+1))

	print('computing overall image mean')
	mean = iter_mean(filename,starty,endy,startx,endx)

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


