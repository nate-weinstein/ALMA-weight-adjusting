def new_var_vis(file,collapse=False):
    ''' Calculate the weight based on the variance in a visibility map at each u,v point and each channel. The codes estimate the variance among the 50 closest uv-points in a limited range in uv-space. 

    :param file:
    Name of the visibility fits file for which the weights will be calculated. The file is assumed to contain a spectral line (ie. it contains a spectral dimension). Use var_vis_cont if you have continuum data that has been averaged along the spectral dimension.
    
    :param collapse:
    Calculates the average across the spectral windows, rather than calculating an average in each spectral window separately. This is necessary if you are using line-free channels to calculate the dispersion.
    ''' 
    from astropy.io import fits
    import numpy as np

    fits_file = fits.open(file)

    vis = (fits_file[0].data['data']).squeeze()
    fits_file.close()

    chi = 0
    if collapse:
        chi += ((vis[i,0,0]**2)*vis[:,0,2]).sum() + ((vis[:,0,1]**2)*vis[:,0,2]).sum() + ((vis[:,1,0]**2)*vis[:,1,2]).sum() + ((vis[:,1,1]**2)*vis[:,1,2]).sum()
        if np.shape(vis[:,:,2]) > np.shape(vis[:,:,2]!=0): #Will happen with weights of zero as well, which usually are due to u or v =0 (unphysical points)
            print("Alert: Non-positive weights in data set")
        n_els=len(np.ravel(vis[:,:,2]!=0))
    else:
        for i in range(len(vis[0,:])):
            chi += ((vis[:,i,0,0]**2)*vis[:,i,0,2]).sum() + ((vis[:,i,0,1]**2)*vis[:,i,0,2]).sum() + ((vis[:,i,1,0]**2)*vis[:,i,1,2]).sum() + ((vis[:,i,1,1]**2)*vis[:,i,1,2]).sum()
        if np.shape(vis[:,:,:,2]) > np.shape(vis[:,:,:,2]!=0):
            print("Alert: Non-positive weights in data set")
        n_els=len(np.ravel(vis[:,:,:,2]!=0))

    red_chi = chi/(2*n_els)

    return red_chi

def adjust_weights(file,red_chi,collapse=False):
    ''' Scales the weights in given file by a factor of red_chi

    :param file:
    Name of the visibility fits file for which the weights will be calculated. The file is assumed to contain a spectral line (ie. it contains a spectral dimension). Use var_vis_cont if you have continuum data that has been averaged along the spectral dimension.

    :param collapse:
    Calculates the average across the spectral windows, rather than calculating an average in each spectral window separately. This is necessary if you are using line-free channels to calculate the dispersion.
    '''
    from astropy.io import fits
    import numpy as np

    fits_file = fits.open(file)
    if collapse:
        fits_file[0].data['data'][:,0,0,0,:,2]/=red_chi
    else:
        fits_file[0].data['data'][:,0,0,0,:,:,2]/=red_chi
    fits_file.writeto(file,overwrite=True)
    fits_file.close()
    

red_chi=new_var_vis('file_name')
adjust_weights('file_name',red_chi)
