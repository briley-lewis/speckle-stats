
import sys
sys.path.insert(0,'/u/home/b/blewis34/MEDIS/Proper')
sys.path.insert(0,'/u/home/b/blewis34/MEDIS/Proper/proper')
sys.path.insert(0,'/u/home/b/blewis34/MEDIS/medis')
sys.path.insert(0,'/u/home/b/blewis34/MEDIS')
import medis.save_photon_data as gpd
from medis.params import ap, tp, sp, iop, cp
import numpy as np
import matplotlib.pyplot as plt
import os

##set random seed
np.random.seed(seed=0)

##set astro parameters
ap.companion=True #is there a planet in there?
ap.sample_time = 0.005 #exposure time in [s]
ap.numframes = 200
ap.grid_size = 128*2
ap.contrast = [5e-2]
ap.star_spec = None
#ap.nwsamp = 1 #number of wavefronts in proper to sample from

##set telescope parameters
tp.use_spiders = False #spiders off/on
tp.obscure = False #no spiders, no secondary
tp.use_ao = False #AO off/on
#tp.servo_error = [0.005,1] ##only relevant if trying AO!
tp.piston_error = False #piston error off/on
tp.use_coron = True #coronagraph off/on
tp.occulter_type = 'Vortex' #occulter type - vortex, none, gaussian, 8th_Order, None (Lyot Stop)
tp.detector = 'ideal'
tp.beam_ratio = 25/64./2
tp.check_args()

##set simulation parameters
sp.show_wframe = False
sp.save_obs = False
sp.num_processes = 4
sp.save_locs = np.array([['add_atmos',], ['coronagraph',]])
sp.return_E = True
sp.show_cube = False

#atmos params -- all currently set to defaults
cp.model = 'single'  # single|hcipy_standard|evolving
cp.show_caosparams= True  # for control over all other variables
cp.cn = 1e-13 #lower cn2 should mean less turbulence - this value is equivalent to r0~30 (>>D)
cp.L0 = 10 #longer L0 (coherence length) means less turbulent, I think
cp.v = np.asarray([5,0])
cp.h = 100

iop.update("complex-fields/")
###SET FILE SAVE NAME
iop.fields = os.path.join(iop.testdir, 'MEDIS_1sec_Aug2020.h5')

if __name__ == '__main__':
	gpd.run_medis()
