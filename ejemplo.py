from PyIFSview.PyIFSview import PyIFSview
from astropy.io import fits
import numpy as np
#from PyIFSview.der_snr import der_snr


#test(cube_fits_file = 'NGC2253.fits.gz')
#test(cube_fits_file = 'J0210_cube.fits')
#test(cube_fits_file = 'NGC2906.fits.fz', origin_cube='MUSE')
PyIFSview(cube_fits_file = 'manga-7443-12703-LOGCUBE.fits.gz', origin_cube='MANGA')
#test(cube_fits_file = '37050_blue_8_Y13SAR1_P005_15T018.fits.gz', origin_cube='SAMI')


#image_data = hdu_list['FLUX'].data

# hdu_list = fits.open(cube_fits_file) 

# image_data = hdu_list[1].data
# image_data[image_data == 1] = np.nan # convert 0 to nan
# naxis1 = hdu_list[1].header['NAXIS1']
# naxis2 = hdu_list[1].header['NAXIS2']
# naxis3 = hdu_list[1].header['NAXIS3']
# crval3 = hdu_list[1].header['CRVAL3']
# cdelt3 = hdu_list[1].header['CD3_3']

# #print(hdu_list[0].header)

# # Read the error spectra if available. Otherwise estimate the errors with the der_snr algorithm
# if len(hdu_list) == 3:
#     image_error  = hdu_list[2].data
# elif len(hdu_list) == 2:
#     print("No error extension found. Estimating the error spectra with the der_snr algorithm")
#     image_error = np.zeros( image_data.shape )
#     for i in range( 0, image_data.shape[1] ):
#         image_error[:,i] = der_snr( image_data[:,i] )
