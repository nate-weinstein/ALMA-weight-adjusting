def new_var_vis(file,collapse=False,realimag=False):
    ''' Calculate the weight based on the variance in a visibility map at each u,v point and each channel. The codes estimate the variance among the 50 closest uv-points in a limited range in uv-space. 

    :param file:
    Name of the visibility fits file for which the weights will be calculated. The file is assumed to contain a spectral line (ie. it contains a spectral dimension). Use var_vis_cont if you have continuum data that has been averaged along the spectral dimension.
    
    :param collapse:
    Calculates the average across the spectral windows, rather than calculating an average in each spectral window separatly. This is necessary if you are using line-free channels to calculate the dispersion.

    :param realimag:
    Estimates one weight for both the real and imaginary part, based exclusively on the real part of the visibilities. Otherwise a separate weight is calculated for the real and imaginary part of the visibilities. 
    ''' 
    from astropy.io import fits
    import numpy as np

    fits_file = fits.open(file)

    vis = (fits_file[0].data['data']).squeeze()
    fits_file.close()

    import time
    start=time.time()
    chi = 0
    for i in range(len(vis[0,:])):
        chi += ((vis[:,i,0,0]**2)*vis[:,i,0,2]).sum() + ((vis[:,i,0,1]**2)*vis[:,i,0,2]).sum() + ((vis[:,i,1,0]**2)*vis[:,i,1,2]).sum() + ((vis[:,i,1,1]**2)*vis[:,i,1,2]).sum()
    red_chi = chi/(len(vis[:])*2*2*len(vis[0,:,0,0]))
    
    file1 = open("new_var_vis_output.txt","w")
    file1.write('Chi square = '+str(chi))
    file1.write('Reduced chi square = '+str(red_chi))
    file1.write('Elapsed time (hrs): '+str((time.time()-start)/3600.))
    file1.close()

    return red_chi

def adjust_weights(file,red_chi):
    from astropy.io import fits
    import numpy as np

    fits_file = fits.open(file)
    fits_file[0].data['data'][:,0,0,0,:,:,2]/=red_chi
    fits_file.close()


red_chi=new_var_vis('f1_HD36546_data.uvfits')
#adjust_weights('f1_HD36546_data.uvfits',red_chi)
#Can edit this file name in vim on the cluster
#Write to a file instead of printing the chi square, reduced chi square, etc
#Use screen to run in background on cluster
#Run this all in ipython: %run new_var_vis
