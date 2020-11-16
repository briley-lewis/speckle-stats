import math
import numpy as np
import h5py
import matplotlib.pyplot as plt
from astropy.io import fits
from photutils import CircularAperture
from photutils import aperture_photometry
import glob
import sys

filestring = 'MEDIS_short_Aug2020_WITHATTRS_avg-res*.fits'
outname = 'MEDIS-short-cn14-SNR'

####params used for stKLIP:
#modes = [0,1,2,3,4,5,6,7,10,15,20,30,50,60,75,85,100,110,125,140,150,160,175,185,200,300]
#lags = [0,1,2,3,4,5,6,8,10]
#window = [75,175,100,200]
modes = [0,1,2,3,4,5,6,7,10,15,20,30,50,75,100,150,200,300]
lags = [0,2,4,6,8]
window = [75,175,100,200]

def run_plot_SNR(filestring,outname,modes,lags,window):

	##get files
	files = glob.glob('/Volumes/Backup-Plus/stklip-results/{}'.format(filestring))
	print(files)
	data = []
	for f in files:
		d = fits.getdata(f)
		data.append(d)
	data = np.asarray(data)

	#data quicklook
	plt.imshow(data[1,4,:,:],vmax=0.2)
	vy = 27
	vx = 54
	w = 5
	planet_mark = plt.Circle((vx,vy), w,facecolor='None',edgecolor='red')
	ax = plt.gca()
	ax.add_artist(planet_mark)
	plt.xlim(0,99)
	plt.ylim(0,99)
	plt.xlabel("x (pixels)")
	plt.ylabel("y (pixels)")
	plt.title('Input MEDIS Image')
	plt.colorbar()
	plt.savefig('{}_input-img.eps'.format(outname))
	plt.show()

	data_ok = query_yes_no('Does this data look ok (not weird, planet is centered in circle)?',default=None)
	if data_ok==False:
		sys.exit("Something's wrong with your data")
	else:
		print('...contuinuing...')

	#create annulus values for SNR calc
	r = np.zeros((256,256))
	for i in np.arange(0,256):
		for j in np.arange(0,256):
			r[i,j] = np.sqrt((i-128)**2 + (j-128)**2)

	#determine relevant planet values
	px = window[2]+vx
	py = window[0]+vy
	rp = np.sqrt((px-128)**2 + (py-128)**2)

	#compute apertures
	test = data[2,5,:,:]
	apertures=[]
	tol = 15
	w=5
	for i in np.arange(0,256):
		for j in np.arange(0,256):
			#print(r[i,j])
			if r[i,j]>rp-0.1 and r[i,j]<rp+0.1:
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
						flux = aperture_photometry(test,aper)['aperture_sum']
						apertures.append(flux)
					
	plt.imshow(test)
	planet_mark = plt.Circle((vx,vy), w,facecolor='None',edgecolor='yellow')
	ax = plt.gca()
	ax.add_artist(planet_mark)
	plt.colorbar()
	plt.show()

	data_ok = query_yes_no('Do the apertures look ok?',default=None)
	if data_ok==False:
		sys.exit("Error")
	else:
		print('...continuing to SNR calculation...')

	#calculate over all lags and modes
	results = []
	for l in range(0,len(lags)):
		lag = lags[l]
		im = []
		for m in range(0,len(modes)):
			mode = modes[m]
			sto = data[l,m,:,:]
			result = SNR(sto,vx,vy,window)
			im.append(result)
		im = np.asarray(im)
		results.append(im)
	results=np.asarray(results)
	#print(np.shape(results),"[lags,modes,[noise,signal,SNR]]")
	#[lag,mode,noise,signal,SNR]

	##make plots
	#SNR vs. lags
	averaged_over_lag = np.mean(results,axis=1)
	max_over_lag = np.max(results,axis=1)
	plt.scatter(lags,averaged_over_lag[:,2],label='Average SNR')
	plt.scatter(lags,max_over_lag[:,2],label='Peak SNR')
	plt.xlabel('Lags')
	plt.ylabel('SNR')
	plt.legend()
	plt.savefig('{}_lags-v-snr.png'.format(outname))
	plt.show()
	#modes vs. SNR
	cmap = plt.cm.get_cmap('plasma')
	for l in range(0,len(lags)):
		if l==1:
			pass
		else:
			lag = lags[l]
			plt.plot(modes,results[l,:,2],color=cmap(l*60),label="{} lags".format(lag))
	plt.plot(modes,np.full_like(results[0,:,2],results[1,0,2]),linestyle='dashed',color='black',label='Original')
	plt.xlabel('Number of KL Modes')
	plt.ylabel('SNR')
	#plt.ylim(0,60)
	plt.legend(loc=4)
	plt.savefig('{}_modes-v-snr.png'.format(outname))
	plt.show()
	#cropped modes v. SNR
	for l in range(0,len(lags)):
		if l==1:
			pass
		else:
			lag = lags[l]
			plt.plot(modes,results[l,:,2],color=cmap(l*60),label="{} lags".format(lag))
	plt.plot(modes,np.full_like(results[0,:,2],results[1,0,2]),linestyle='dashed',color='black',label='Original')
	plt.xlabel('Number of KL Modes')
	plt.ylabel('SNR')
	#plt.xlim(-1,50)
	#plt.ylim(10,50)
	plt.legend(loc=4)
	plt.savefig('{}_modes-v-snr-cropped.png'.format(outname))
	plt.show()
	#SNR comparison
	orig = SNR(data[1,0,:,:],vx,vy,window)
	for l in range(0,len(lags)):
		if l==1:
			pass
		else:
			lag = lags[l]
			plt.plot(modes,np.abs(results[l,:,2])/orig[2],color=cmap(l*60),label="{} lags".format(lag))

	plt.xlabel('Number of KL Modes')
	plt.ylabel('SNR Comparison')
	plt.legend()
	plt.savefig('{}_modes-v-snr-comparison.png'.format(outname))
	plt.show()

	##now for contrast curves
	###contrast curve
	##no coron for comparison
	noc = h5py.File('/Users/blewis/Desktop/grad-school/Research/Speckle-Stats/Outputs/MEDIS_1sec-nc_Aug2020.h5','r')
	img = noc['data']
	summed = np.sum(img,axis=3)
	nocoron = np.abs(summed[:,1,0,:,:])**2
	#find best aperture fit for star
	plt.imshow(nocoron[0,:,:])
	star_mark = plt.Circle((128,128), 5,facecolor='None',edgecolor='red')
	ax = plt.gca()
	ax.add_artist(star_mark)
	plt.colorbar()
	plt.show()
	data_ok = query_yes_no('Is the star centered in the aperture?',default=None)
	if data_ok==False:
		sys.exit("Need to refine params")
	else:
		print('...contuinuing...')
	##measure star flux
	avg_nocoron = np.mean(nocoron,axis=0)
	star = CircularAperture((128,128), 5)
	star_flux = aperture_photometry(avg_nocoron,star)['aperture_sum']
	mean_flux = star_flux/star.area
	#contrast curve for no-processing original
	pix,contrasts_orig = contrast_curve(data[1,0,:,:],mean_flux)
	#contrast curves for each lag
	for l in range(0,len(lags)):
		lag = lags[l]
		cons=[]
		for m in range(0,len(modes)):
			mode = modes[m]
			pix,con = contrast_curve(data[l,m,:,:],mean_flux)
			plt.semilogy(pix,con,color=cmap(m*15),label="{} Modes".format(mode))
		plt.title("{} Lags Contrast Curve".format(lag))
		plt.legend(bbox_to_anchor=(1, 1.10))
		plt.ylabel('Contrast')
		plt.xlabel('Pixels from center')
		plt.savefig('{}_contrast_{}-lag.png'.format(outname,lag))
		plt.show()
	#contrast curvey comparison for each lag
	for l in range(0,len(lags)):
		lag = lags[l]
		cons=[]
		for m in range(0,len(modes)):
			mode = modes[m]
			pix,con = contrast_curve(data[l,m,:,:],mean_flux)
			con=np.asarray(con)
			contrasts_orig=np.asarray(contrasts_orig)
			plt.plot(pix,con/contrasts_orig,color=cmap(m*15),label="{} Modes".format(mode))
		plt.title("{} Lags Contrast Curve".format(lag))
		plt.legend(bbox_to_anchor=(1, 1.10))
		plt.ylabel('Contrast Improvement')
		plt.xlabel('Pixels from center')
		plt.savefig('{}_contrast_comparison_{}-lag.png'.format(outname,lag))
		plt.show()
	#average contrast over lags
	for l in range(0,len(lags)):
		lag = lags[l]
		cons=[]
		for m in range(0,len(modes)):
			mode = modes[m]
			pix,con = contrast_curve(data[l,m,:,:],mean_flux)
			cons.append(con)
		con = np.mean(np.asarray(cons),axis=0)
		plt.semilogy(pix,con,color=cmap(l*60),label="{} lags".format(lag,mode))
	plt.legend()
	#	plt.ylim(1e-5,1e-3)
	plt.ylabel('Contrast')
	plt.xlabel('Pixels from center')
	plt.savefig('{}_averaged-over-modes-contrast_{}-lag.png'.format(outname,lag))
	plt.show()
	for l in range(0,len(lags)):
		lag = lags[l]
		cons=[]
		for m in range(0,len(modes)):
			mode = modes[m]
			pix,con = contrast_curve(data[l,m,:,:],mean_flux)
			cons.append(con)
		con = np.mean(np.asarray(cons),axis=0)
		plt.plot(pix,con/contrasts_orig,color=cmap(l*60),label="{} lags".format(lag,mode))
	plt.legend()
	plt.ylabel('Contrast Improvement')
	plt.xlabel('Pixels from center')
	plt.savefig('{}_averaged-over-modes-contrast-improvement_{}-lag.png'.format(outname,lag))
	plt.show()

	##let's look at the total image variance
	for l in range(0,len(lags)):
		if l==1:
			pass
		else:
			lag = lags[l]
			var = []
			for m in range(0,len(modes)):
				mode = modes[m]
				image = data[l,m,:,:]
				variance = np.var(image)
				var.append(variance)
			plt.semilogy(modes,var,color=cmap(l*50),label="Variance for {} Lags".format(lag))
	plt.legend()
	plt.ylabel('Variance')
	plt.xlabel('KL Modes')
	plt.savefig('{}_variance-v-modes.png'.format(outname))
	plt.show()

	print('SNR calculated and plots made with name {}'.format(outname))
	return [lags,modes,results]

def contrast_curve(data,mean_flux,window=[75,175,100,200],w=5):
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

#ok, making succint SNR definition to make this easy
def SNR(data,vx,vy,window,tol=15,w=5):
	
	r = np.zeros((256,256))
	for i in np.arange(0,256):
		for j in np.arange(0,256):
			r[i,j] = np.sqrt((i-128)**2 + (j-128)**2)

	apertures=[]
	
	##planet and image parameters
	px = window[2]+vx
	py = window[0]+vy
	rp = np.sqrt((px-128)**2 + (py-128)**2)
	#print(px,py,rp)
	
	for i in np.arange(0,256):
		for j in np.arange(0,256):
			#print(r[i,j])
			if r[i,j]>rp-0.1 and r[i,j]<rp+0.1:
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
						
	#plt.imshow(data,vmax=0.1)
	#planet_mark = plt.Circle((vx,vy), w,facecolor='None',edgecolor='yellow')
	#ax = plt.gca()
	#ax.add_artist(planet_mark)
	#plt.title('SNR Apertures for {} lags, {} modes'.format(lag,mode))
	#plt.colorbar()
	#plt.show()
	
	noise = np.asarray(np.average(apertures))
	pa = CircularAperture((vx,vy), w)
	signal = np.asarray(aperture_photometry(data,pa)['aperture_sum'])
	SNR = signal/noise
	#print(SNR)
	
	return np.asarray([noise,signal,SNR])

def query_yes_no(question, default="yes"):
	"""Ask a yes/no question via input() and return their answer.

	"question" is a string that is presented to the user.
	"default" is the presumed answer if the user just hits <Enter>.
		It must be "yes" (the default), "no" or None (meaning
		an answer is required of the user).

	The "answer" return value is True for "yes" or False for "no".
	"""
	valid = {"yes": True, "y": True, "ye": True,
			 "no": False, "n": False}
	if default is None:
		prompt = " [y/n] "
	elif default == "yes":
		prompt = " [Y/n] "
	elif default == "no":
		prompt = " [y/N] "
	else:
		raise ValueError("invalid default answer: '%s'" % default)

	while True:
		sys.stdout.write(question + prompt)
		choice = input().lower()
		if default is not None and choice == '':
			return valid[default]
		elif choice in valid:
			return valid[choice]
		else:
			sys.stdout.write("Please respond with 'yes' or 'no' "
							 "(or 'y' or 'n').\n")