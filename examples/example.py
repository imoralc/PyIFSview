from der_snr import der_snr
from astropy.io import fits
import numpy as np
from PyIFSview import test

#test(cube_fits_file = 'NGC2253.fits.gz')

cube_fits_file = 'NGC2253.fits.gz'

hdu_list = fits.open(cube_fits_file) 

image_data = hdu_list[0].data
image_data[image_data == 0] = np.nan # convert 0 to nan
naxis1 = hdu_list[0].header['NAXIS1']
naxis2 = hdu_list[0].header['NAXIS2']
naxis3 = hdu_list[0].header['NAXIS3']
crval3 = hdu_list[0].header['CRVAL3']

try:
    cdelt3 = hdu_list[0].header['CDELT3']
except:
    cdelt3 = hdu_list[0].header['CD3_3'] 

# Original:
wave2 = np.arange(crval3, crval3+cdelt3*naxis3, cdelt3) #wavelength

                #image_error = hdu_list[1].data

                # Read the error spectra if available. Otherwise estimate the errors with the der_snr algorithm
print("No error extension found. Estimating the error spectra with the der_snr algorithm")
image_error = np.zeros( image_data.shape )
for i in range( 0, image_data.shape[1] ):
    image_error[:,i] = der_snr( image_data[:,i] )
print(image_error)

test(cube_fits_file = 'manga-7443-12703-LOGCUBE.fits.gz', origin_cube='MANGA')