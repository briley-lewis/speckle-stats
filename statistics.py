##step 1: import relevant modules
import os
import numpy as np
import matplotlib.pyplot as plt
import pylab
import h5py
import math
from astropy.io import fits
import pickle

def sum_wavefronts(f_in,to_intensity=False,load=True):
    """sums over # of wavefronts in each timestep for a given hdf5 file, returns hdf5 file.
    NOTE: attributes are currently lost in transition, will add in future update. in the meantime,
    make sure your f_in is descriptive"""
    #open file
    #f_in = h5py.File(f_in, 'r')

    #make list of keys in correct order
    #n_screens = np.arange(0,len(f_in.keys()))
    #keys = ['t' + str(n) for n in n_screens]
    #print('starting with',keys[0])

    #if to_intensity==True:
        #would be cool to implement functionality where I could save intensities, 
        #but then that wouldn't be compatible with my other definitions so maybe not
       #pass

    f_out = h5py.File('summed'+f_in, 'w')

    if load=True:
        f = h5py.File(f_in)
        data = f['data']
        summed = np.sum(data,axis=3)
        print(np.shape(summed))
    
    #iterative not fully implemented yet
    #with f = h5py.File(f_in):
        #data_set = f['data']
        #shape = np.shape(data_set)
            #for i in range(len(shape[0])):
                #frame=data_set[i,:,:,:,:,:]
                #frame=np.sum(frame,axis=3)
    
    f_out.create_dataset(key,data=frame)

    print('data saved with dimensions',summed)
    print('save name is','summed-'+f_in+'.h5')

    f_out.close()

def load_full(f_in,summed=True,start=100,end=150):
    """loads in full data set into an array -- do NOT use for large files! 
    returns array of intensities and number of timesteps"""
    #open file
    f_in = h5py.File(f_in, 'r')

    #make list of keys in correct order
    n_screens = np.arange(0,len(f_in.keys()))
    keys = ['t' + str(n) for n in n_screens]
    ##error catching
    if np.ndim(f_in[keys[0]])>4:
        summed=False

    #get intensities
    intensities=[]
    for key in keys:
       frame = f_in[key]
       if summed==False:
          frame = np.sum(frame,axis=2)
       intensity = np.abs(frame[1,0,start:end,start:end])**2
       intensities.append(intensity)
    intensities = np.asarray(intensities)

    return intensities,n_screens

def full_mean(f_in,summed=True,start=100,end=150):
    """requires input of .hdf5 file summed over wavefront axis, with option to use not-summed file"""
    ##error catching
    #if np.ndim(f_in[keys[0]])>4:
        #summed=False

    if "fits" in f_in:
        intensities = fits.getdata(f_in)

    else:
        nx,ny=[end-start,end-start]
        intensities, n_screens = load_full(f_in,summed,start,end)

    mean = np.mean(intensities,axis=0)

    return mean


def iter_mean(f_in,summed=True,start=100,end=150):
    """get the mean image from an hdf5 dataset, iterative approach, with option to use not-summed file"""
    #open file
    f_in = h5py.File(f_in, 'r')
    
    #make list of keys in correct order
    n_screens = np.arange(0,len(f_in.keys()))
    keys = ['t' + str(n) for n in n_screens]
    ##error catching
    if np.ndim(f_in[keys[0]])>4:
        summed=False
    
    #iterate through to compute mean
    mn_im = np.zeros((end-start,end-start)) #initializing array with same spatial dimensions
    k = 0
    for key in keys:
       frame = f_in[key]
       if summed==False:
           frame = np.sum(frame,axis=2)
       intensity = np.abs(frame[1,0,start:end,start:end])**2
       k = k+1
       mn_im = ((k-1)*mn_im + intensity)/k
    return mn_im

def full_spatial_cov(f_in,summed=True,start=100,end=150):
    """get spatial covariance (no time lags) from an hdf5 data set, with option to use file not summed over 
    number of wavefronts"""
    ##error catching
    #if np.ndim(f_in[keys[0]])>4:
        #summed=False

    intensities, n_screens = load_full(f_in,summed,start,end)
    nx,ny=[end-start,end-start]
    #print(n_screens,nx,ny)

    mean = np.mean(intensities,axis=0)
    mean_subtracted = np.reshape(intensities-mean[np.newaxis,:,:], (len(n_screens), ny*nx))
    cov = np.dot(mean_subtracted.T, mean_subtracted)/(len(n_screens)-1)

    stdev = np.sqrt(np.diagonal(cov))
    norm = np.outer(stdev,stdev)
    corr = cov/norm

    return corr

def iter_spatial_cov(f_in,summed=True,start=100,end=150,return_mean=False):
    """get spatial covariance (no time lags) from an hdf5 data set, iterative approach, must used summed file
    NOTE: this currently WILL break if you use a large enough file, because of the implementation of correlation
    normalization relying on the full array"""

    ##error catching
    #if np.ndim(f_in[keys[0]])>4:
        #summed=False

    if summed==False:
       print("ERROR: you must use a file summed over your # of wavefronts!")
    
    else:
        #open file
        f_in = h5py.File(f_in, 'r')
        
        #make list of keys in correct order
        n_screens = np.arange(0,len(f_in.keys()))
        keys = ['t' + str(n) for n in n_screens]
        nx,ny=[end-start,end-start]

        #iterate through to compute mean and covariance
        n = nx*ny
        mn_im = np.zeros((n,)) # initialize
        cov = np.empty((n,n))
        vv = np.empty((n,n))
        intensities=[]
        k = 0
        for key in keys:
           frame = f_in[key]
           intensity = np.abs(frame[1,0,start:end,start:end])**2
           intensity = np.reshape(intensity,n)
           intensities.append(intensity)
           k = k+1
           mn_im = ((k-1)*mn_im + intensity)/k
           if k == 1:
              #start with regular covariance
              v = intensity-mn_im
              cov = np.outer(v,v)
           else:
              #for all other timesteps, use iterative calculation
              v = intensity-mn_im
              vv = ((k)/(k-1)**2)*np.outer(v,v)
              cov *= (k-2)/(k-1)
              cov += vv

        stdev = np.sqrt(np.diagonal(cov))
        norm = np.outer(stdev,stdev)
        corr = cov/norm
              
        if return_mean==True:
           return np.reshape(mn_im,[nx,ny]), corr
        else:
           return corr

def full_spacetime_cov(f_in,spatial_cov,lag,summed=True,start=100,end=150):
    """get spatial covariance for a given lag from an hdf5 data set, with option to use file not summed over 
    number of wavefronts"""

    ##error catching
    #if np.ndim(f_in[key[0]])>4:
        #summed=False

    intensities, n_screens = load_full(f_in,summed,start,end)
    #print(np.shape(intensities))
    n_screen, nx, ny = np.shape(intensities)

    mean = np.mean(intensities,axis=0)
    mean_subtracted = np.reshape(intensities-mean[np.newaxis,:,:], (n_screen, ny*nx))

    cov = np.dot(mean_subtracted[0:n_screen-lag,:].T, mean_subtracted[lag:,:])/(n_screen-1)

    stdev = np.sqrt(np.diagonal(spatial_cov))
    norm = np.outer(stdev,stdev)
    corr = cov/norm
    print('maximum of full st cov (normalized) is',np.max(corr))
    corr = np.reshape(corr,[nx,ny,nx,ny])

    return corr


###unfinished below this point!!!!
def iter_spacetime_cov(f_in,spatial_cov,lag,summed=True,start=100,end=150,return_mean=False):
    """get spatial covariance for a given time lag from an hdf5 data set, iterative approach, 
    must used summed file     
    NOTE: this currently WILL break if you use a large enough file, because of the implementation of correlation
    normalization relying on the full array"""

    if summed==False:
        print("ERROR: you must use a file summed over your # of wavefronts!")

    else:
        #open file
        f_in = h5py.File(f_in, 'r')
        
        #make list of keys in correct order
        n_screens = np.arange(0,len(f_in.keys()))
        keys = ['t' + str(n) for n in n_screens]
        nx,ny=[end-start,end-start]

        ##error catching
        if np.ndim(f_in['t0'])>4:
            summed=False
            print('WRONG FILE USED -- ABORT')

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
                    #print('step done')

        stdev = np.sqrt(np.diagonal(spatial_cov))
        norm = np.outer(stdev,stdev)
        corr = cov/norm

        #return outputs
        if return_mean:
            return mn_im, corr
        else:
            print('maximum of iterative st cov (normalized) is',np.max(corr))
            corr = np.reshape(corr,[nx,ny,nx,ny])
            return corr

