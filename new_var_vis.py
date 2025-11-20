def new_var_vis(collapse=False,realimag=False):
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
    file = input("Input file path: ")
    fits_file = fits.open(file)
    #u,v = fits_file[0].data['UU'],fits_file[0].data['VV']
    #freq0 = fits_file[0].header['crval4']
    #klam = freq0/1e3

    vis = (fits_file[0].data['data']).squeeze()
    fits_file.close()

    import time
    start=time.time()
    #weight = np.zeros((u.size(),(vis.shape)[1],2))
    chi = 0
    for i in range():
        chi += ((vis[:,i,0,0]**2)*vis[:,i,0,2]).sum() + ((vis[:,i,0,1]**2)*vis[:,i,0,2]).sum() + ((vis[:,i,1,0]**2)*vis[:,i,1,2]).sum() + ((vis[:,i,1,1]**2)*vis[:,i,1,2]).sum()
    red_chi = chi/(len(vis[:])*2*2*vis[0,:,0,0].size())
    print('Chi square = ',str(chi))
    print('Reduced chi square = ',str(red_chi))
    print('Elapsed time (hrs): ',(time.time()-start)/3600.)