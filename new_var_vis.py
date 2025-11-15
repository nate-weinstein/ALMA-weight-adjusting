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
    u,v = fits_file[0].data['UU'],fits_file[0].data['VV']
    freq0 = fits_file[0].header['crval4']
    klam = freq0/1e3
    
    vis = (fits_file[0].data['data']).squeeze()
    