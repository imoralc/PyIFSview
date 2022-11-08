################################################################################

"""
Copyright (C) 2022, Ignacio del Moral Castro
E-mail: ignaciodelmoralcastro@gmail.com

################################################################################

NAME:

PyIFSview()

PURPOSE:

PyIFSview is a tool developed to interactively visualize integral field spectroscopy (IFS) data, 
such as CALIFA, MaNGA, SAMI or MUSE data.

CALLING SEQUENCE:

from PyIFSview import PyIFSview

PyIFSview(cube_fits_file, origin_cube=None, slide=None, f_min =None, f_max =None, l_min =None, 
          l_max =None, c_min =None, c_max =None, x0=None, y0=None, snr_lmax=None, snr_lmin=None)


INPUT PARAMETERS

- `cube_fits_file`: Path and name of the IFS datacube to read

KEYWORDS:

- origin_cube: str optional. This keyword allow us to indicate the origin of the data.
  Possible options: CALIFA, MUSE, MaNGA, SAMI, KOALA or CAVITY

- slide: float optional.

- f_min and f_max: float optional. Flux limits of the spectrum's plot.

- l_min` and l_max: float optional. Wavelenght limits to show the spectrum

- c_min` and c_max: float optional. Flux limits of the imshow's plot.

- x0 and y0: int optional. Coordinates of the spectrum to show.

- snr_lmax and snr_lmin: float optional. Wavelenght limits to compute the S/N


REQUIRED ROUTINES:

- Numpy: Python library used for working with arrays
- Matplotlib: Python 2D plotting library
- Astropy: astronomy library

NODIFICATION HISTORY:

V1.0.0 -- Created by Ignacio del Moral-Castro

"""

########################################################################
#     importing libraries                                              #
########################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox
from astropy.io import fits
import matplotlib.gridspec as gridspec
from PyIFSview.der_snr import der_snr

# Update fuego color map
import matplotlib.colors as colors
fuego_color_map = colors.LinearSegmentedColormap.from_list("fuego", ((0.25, 0, 0),  (0.5,0,0),    (1, 0, 0), (1, 0.5, 0), (1, 0.75, 0), (1, 1, 0), (1, 1, 1)), N=256, gamma=1.0)
fuego_color_map.set_bad('lightgray')
plt.register_cmap(cmap=fuego_color_map)


class PyIFSview(object):
    
    def __init__(self, cube_fits_file, origin_cube=None, slide=None, f_min =None, f_max =None, l_min =None, l_max =None, c_min =None, c_max =None, x0=None, y0=None, snr_lmax=None, snr_lmin=None):

        ### Checking initial parameters

        self.verbose= False

        if slide != None:
            self.slide = slide
        else:
            self.slide = 0

        if c_min != None:
            self.c_min = c_min
        else:
            self.c_min = -1.

        if c_max != None:
            self.c_max = c_max
        else:
            self.c_max = 0.

        if f_min != None:
            self.f_min = f_min
        else:
            self.f_min = 0.

        if f_max != None:
            self.f_max = f_max
        else:
            self.f_max = 0.

        if l_min != None:
            self.l_min = l_min
        else:
            self.l_min = 0.

        if l_max != None:
            self.l_max = l_max
        else:
            self.l_max = 0.

        if x0 != None:

            self.x0 = x0
        else:
            self.x0 = 0.

        if y0 != None:
            self.y0 = y0
        else:
            self.y0 = 0.

        self.cmap0 =""
        self.color_list = ['jet', 'gist_gray','viridis', 'gnuplot', 'gnuplot2', 'cubehelix', 'nipy_spectral', 'RdBu', 'fuego']
        
        ###################### 

        ### Open cube

        if origin_cube == 'CALIFA':
            self.open_CALIFA(cube_fits_file, snr_lmin, snr_lmax)
        elif origin_cube == 'MUSE':
            self.open_MUSE(cube_fits_file)
        elif origin_cube == 'MANGA':
            self.open_MANGA(cube_fits_file)
        elif origin_cube == 'SAMI':
            self.open_SAMI(cube_fits_file)
        elif origin_cube == 'CAVITY':
            self.open_CAVITY(cube_fits_file)
        elif origin_cube == None:
            self.hdu_list = fits.open(cube_fits_file) 
            
            try:
                self.image_data = self.hdu_list[0].data
                self.image_data[self.image_data == 0] = np.nan # convert 0 to nan
                self.naxis1 = self.hdu_list[0].header['NAXIS1']
                self.naxis2 = self.hdu_list[0].header['NAXIS2']
                self.naxis3 = self.hdu_list[0].header['NAXIS3']
                self.crval3 = self.hdu_list[0].header['CRVAL3']
                try:
                    self.cdelt3 = self.hdu_list[0].header['CDELT3']
                except:
                    self.cdelt3 = self.hdu_list[0].header['CD3_3'] 


                # Sum datacube
                self.sum_cube = np.nansum(self.image_data, axis=0)
                self.sum_cube[self.sum_cube == 0] = np.nan

                self.c_max_sum = np.nanmax(np.nansum(self.image_data, axis=0))
                #self.c_min = np.nanmin(np.nansum(self.image_data, axis=0))
                self.c_min_sum = np.nanpercentile(self.sum_cube, 50)

                # Wavelength:
                self.wave2 = np.arange(self.crval3, self.crval3+self.cdelt3*self.naxis3, self.cdelt3) 

                # Read the error spectra if available. Otherwise estimate the errors with the der_snr algorithm
                if len(self.hdu_list) > 2:
                    self.image_error = self.hdu_list[1].data
                elif len(self.hdu_list) <= 2:
                    print("No error extension found. Estimating the error spectra with the der_snr algorithm")
                    self.image_error = np.zeros( self.image_data.shape )
                    for i in range( 0, self.image_data.shape[1] ):
                        self.image_error[:,i] = der_snr( self.image_data[:,i] )

                if snr_lmin != None:
                    self.snr_lmin = snr_lmin
                else:
                    self.snr_lmin = np.nanmin(self.wave2)

                if snr_lmax != None:
                    self.snr_lmax = snr_lmax
                else:
                    self.snr_lmax = np.nanmax(self.wave2)

                self.idx_snr = np.where(np.logical_and( self.wave2 >= self.snr_lmin, self.wave2 <= self.snr_lmax ) )[0]
                self.signal = np.nanmedian(self.image_data[self.idx_snr,:],axis=0)
                self.noise  = np.abs(np.nanmedian(np.sqrt(self.image_error[self.idx_snr,:]),axis=0))
                self.snr    = self.signal / self.noise

                self.c_max_snr = np.nanmax(self.snr)
                #self.c_min = np.nanmin(np.nansum(self.image_data, axis=0))
                self.c_min_snr = np.nanpercentile(self.snr, 60)

            except:
                self.image_data = self.hdu_list[1].data
                self.image_data[self.image_data == 0] = np.nan # convert 0 to nan
                self.naxis1 = self.hdu_list[1].header['NAXIS1']
                self.naxis2 = self.hdu_list[1].header['NAXIS2']
                self.naxis3 = self.hdu_list[1].header['NAXIS3']
                self.crval3 = self.hdu_list[1].header['CRVAL3']
                self.cdelt3 = self.hdu_list[1].header['CD3_3'] 

                # Sum datacube
                self.sum_cube = np.nansum(self.image_data, axis=0)
                self.sum_cube[self.sum_cube == 0] = np.nan

                self.c_max_sum = np.nanmax(np.nansum(self.image_data, axis=0))
                #self.c_min = np.nanmin(np.nansum(self.image_data, axis=0))
                self.c_min_sum = np.nanpercentile(self.sum_cube, 50)

                # Wavelength:
                self.wave2 = np.arange(self.crval3, self.crval3+self.cdelt3*self.naxis3, self.cdelt3) 

                # Read the error spectra if available. Otherwise estimate the errors with the der_snr algorithm
                if len(self.hdu_list) > 2:
                    self.image_error  = self.hdu_list[2].data
                elif len(self.hdu_list) <= 2:
                    print("No error extension found. Estimating the error spectra with the der_snr algorithm")
                    self.image_error = np.zeros( self.image_data.shape )
                    for i in range( 0, self.image_data.shape[1] ):
                        self.image_error[:,i] = der_snr( self.image_data[:,i] )

                if snr_lmin != None:
                    self.snr_lmin = snr_lmin
                else:
                    self.snr_lmin = np.nanmin(self.wave2)

                if snr_lmax != None:
                    self.snr_lmax = snr_lmax
                else:
                    self.snr_lmax = np.nanmax(self.wave2)

                self.idx_snr = np.where(np.logical_and( self.wave2 >= self.snr_lmin, self.wave2 <= self.snr_lmax ) )[0]
                self.signal = np.nanmedian(self.image_data[self.idx_snr,:],axis=0)
                self.noise  = np.abs(np.nanmedian(np.sqrt(self.image_error[self.idx_snr,:]),axis=0))
                self.snr    = self.signal / self.noise

                self.c_max_snr = np.nanmax(self.snr)
                #self.c_min = np.nanmin(np.nansum(self.image_data, axis=0))
                self.c_min_snr = np.nanpercentile(self.snr, 60)

        # Check we don't have negative values (for plotting log)
        minimo = np.nanmin(self.image_data)
        if self.verbose: print("Min value found in this cube:", minimo)
        if minimo < 0: 
            self.image_data= self.image_data + np.abs(minimo) +1E-20
            if self.verbose: 
                print("----> As this value is negative, adopting as min value ",np.nanmin(self.image_data))
        
        


        self.fig = plt.figure(figsize=(16,9))
        spec2 = gridspec.GridSpec(ncols=3, nrows=3)
        self.ax1 = self.fig.add_subplot(spec2[1:3, 0:2])   # Image bottom left
        self.ax2 = self.fig.add_subplot(spec2[0, 0:3])     # Spectrum top row
        spec2.update(hspace=0.45)

        #### figure 2: spectrum

        self.wave = self.wave2

        if snr_lmin != None:
            self.snr_lmin = snr_lmin
        else:
            self.snr_lmin = np.nanmin(self.wave)

        if snr_lmax != None:
            self.snr_lmax = snr_lmax
        else:
            self.snr_lmax = np.nanmax(self.wave)

        if self.slide == 0: 
            self.slide = int(self.naxis3/2) 

        if self.x0 ==0 : 
            self.x0 = int(self.naxis1/2)
        if self.y0 ==0 : 
            self.y0 = int(self.naxis2/2)

        
        if self.cmap0 == "" : 
            self.cmap0 = 'jet'
        elif self.cmap0 not in self.color_list : 
            self.cmap0 = 'jet'

        if self.c_min == -1: 
            self.c_min = np.nanpercentile(self.image_data, 1)
        if self.c_max ==0 : 
            self.c_max = np.nanmax(self.image_data) 

        if self.l_min == 0 :
            self.l_min = min(self.wave) 
        if self.l_max == 0: 
            self.l_max = max(self.wave)
        
        ### Defining first plot with imshow

        self.norm=colors.LogNorm(vmin=self.c_min, vmax=self.c_max)

        self.l = self.ax1.imshow(self.image_data[self.slide], origin='lower', norm=self.norm)
        self.cbar = self.fig.colorbar(self.l, orientation="vertical", pad=0.02, ax=self.ax1)
        self.cbar.set_label('Flux [ erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')

        self.l4, = self.ax1.plot(self.x0, self.y0, '+', color='black')
        
        ### Defining second plot with the spectrum

        self.ax2.set_title('spaxel'+' '+str(self.x0)+','+str(self.y0))
        spectrum = self.image_data[:,self.y0, self.x0]   ### data reads lambda, y0, x0 
        self.l2, = self.ax2.plot(self.wave[:], spectrum, label='emision')
        self.ax2.set_xlim(self.l_min, self.l_max)
        self.ax2.set_xlabel('Wavelength [$\AA$]')
        self.ax2.minorticks_on()
        self.ax2.set_ylabel('Flux [ erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')

        ##### Defining interactive Slide for wavelenght slide

        axcolor = 'lightgoldenrodyellow'
        axfreq = plt.axes([0.125, 0.92, 0.776, 0.05], facecolor=axcolor)

        self.sfreq = Slider(axfreq, '', 0, self.naxis3-2, valinit=self.slide)
        self.sfreq.valtext.set_visible(False)
        
        self.sfreq.on_changed(self.update)

        ##### Defining interactive button for reset
        
        resetax = plt.axes([0.85, 0.2, 0.05, 0.05])
        self.button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
        self.button.on_clicked(self.reset)

        ##### Defining interactive button for cmaps
   
        rax = plt.axes([0.65, 0.20, 0.15, 0.22])
        
        color_list = ['jet', 'gist_gray','viridis', 'gnuplot', 'gnuplot2', 'cubehelix', 'nipy_spectral', 'RdBu', 'fuego']
    
        self.radio = RadioButtons(rax, (color_list), active=color_list.index(self.cmap0))
        self.radio.on_clicked(self.colorfunc)

        ##### Choosing linear / log map / power
        
        scaleax = plt.axes([0.85, 0.25, 0.05, 0.1])
        self.scale = RadioButtons(scaleax, ('linear', 'log', 'power 2'), active=1)
        self.scale.on_clicked(self.scalefunc)

        ###### Choosing cube / integrated / S/N plot

        sscaleax = plt.axes([0.15, 0.25, 0.05, 0.1])
        self.sscale = RadioButtons(sscaleax, ('slide', 'sum', 'snr'), active=0)
        self.sscale.on_clicked(self.select_map)

        ##### Defining interactive limits for cmaps
    
        ax_cmin = plt.axes([0.65, 0.1, 0.19, 0.03])
        ax_cmax  = plt.axes([0.65, 0.15, 0.19, 0.03])
        
        self.color_min =  np.nanpercentile(self.image_data, 1)    # These two must be the absolute min and max, otherwise this points couldn't be reached by slider!!!
        self.color_max =  np.nanmax(self.image_data)
        
        self.s_cmin = Slider(ax_cmin, 'vmin', self.color_min, self.color_max, valinit=self.c_min, valfmt='%5.2E')
        self.s_cmin.valtext.set_visible(False)
        self.s_cmax = Slider(ax_cmax, 'vmax', self.color_min, self.color_max, valinit=self.c_max, valfmt='%5.2E')
        self.s_cmax.valtext.set_visible(False)
        
        self.s_cmin.on_changed(self.update_cmin)
        self.s_cmax.on_changed(self.update_cmax)

        ########  Text boxes with values
        
        self.pcmin = "{:.3e}".format(self.c_min)    
        axbox1 = plt.axes([0.85, 0.1, 0.06, 0.03])
        self.text_box1 = TextBox(axbox1, '', initial=self.pcmin)
        self.text_box1.on_submit(self.input_cmin)

        self.pcmax = "{:.3e}".format(self.c_max)    
        axbox2 = plt.axes([0.85, 0.15, 0.06, 0.03])
        self.text_box2 = TextBox(axbox2, '', initial=self.pcmax)
        self.text_box2.on_submit(self.input_cmax)


         ##### Defining interactive limits for flux in spectrum
    
        if self.f_max == 0: 
            self.f_max = np.nanmax(self.image_data)
        
        ax_fmin = plt.axes([0.65, 0.45, 0.19, 0.03])
        ax_fmax  = plt.axes([0.65, 0.5, 0.19, 0.03])

        
        self.f_min_real = np.nanmin([-np.nanmin(self.image_data),np.nanmin(self.image_data)])   #(np.abs(f_min) -f_max)/10. 
        self.f_max_real = np.nanmax(self.image_data)
            
        self.s_fmin = Slider(ax_fmin, 'Flux$_{min}$', self.f_min_real, self.f_max_real, valinit=self.f_min, valfmt='%5.2E')
        #s_fmin.valtext.set_visible(False)
        self.s_fmax = Slider(ax_fmax, 'Flux$_{max}$', self.f_min_real, self.f_max_real, valinit=self.f_max, valfmt='%5.2E')
        #s_fmax.valtext.set_visible(False)
        
        self.s_fmin.on_changed(self.update_fmin)
        self.s_fmax.on_changed(self.update_fmax)


        ##### Defining interactive limits for wavelength in spectrum
    
        ax_lmin = plt.axes([0.65, 0.55, 0.19, 0.03])
        ax_lmax  = plt.axes([0.65, 0.6, 0.19, 0.03])
        
        self.s_lmin = Slider(ax_lmin, '$\lambda_{min}$ ($\AA$)', np.min(self.wave), np.max(self.wave), valinit=self.l_min, valfmt='%7.2f')
        #s_lmin.valtext.set_visible(False)
        self.s_lmax = Slider(ax_lmax, '$\lambda_{max}$ ($\AA$)', np.min(self.wave), np.max(self.wave), valinit=self.l_max, valfmt='%7.2f')
        #s_lmax.valtext.set_visible(False)
        
        self.s_lmin.on_changed(self.update_lmin)
        self.s_lmax.on_changed(self.update_lmax)

        ##### Defining moving plane by plane 
    
        planedown = plt.axes([0.108, 0.92, 0.017, 0.05])
        planeup = plt.axes([0.901, 0.92, 0.017, 0.05])
        
        self.button_down = Button(planedown, '-', color='lightgoldenrodyellow', hovercolor='0.975')
        self.button_up = Button(planeup, '+', color='lightgoldenrodyellow', hovercolor='0.975')
            
        self.button_up.on_clicked(self.planeup_event)
        self.button_down.on_clicked(self.planedown_event)

        #### Onclick
    
        cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)  

        plt.show()

    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------

    ##### Functions for setting the wavelenght slide

    def update(self,val):
        self.freq = self.sfreq.val
        self.l.set_data(self.image_data[int(self.freq)])
        self.ax1.set_title(str(self.wave[int(self.freq)])+' $\AA$')
        self.fig.canvas.draw_idle()

    def planeup_event(self,event):
        self.sfreq.val = self.sfreq.val+1
        if self.sfreq.val > self.naxis3:
            self.sfreq.val = self.naxis3
        #self.sfreq.set_val(self.slide)
        self.freq = self.sfreq.val
        #self.ptitle = np.str(np.round(self.wave[self.slide],2))+' $\AA$'  
        self.l.set_data(self.image_data[int(self.freq)])
        #self.l3.set_data([self.wave[self.slide], self.wave[self.slide]],[-10,10])
        #ptitle = np.str(np.round(self.wave[self.slide],2))+' $\AA$'    
        self.ax1.set_title(str(self.wave[int(self.freq)])+' $\AA$')
        self.fig.canvas.draw_idle()
            #print("CAN'T GO HIGHER!!")

    def planedown_event(self,event):
        self.sfreq.val = self.sfreq.val-1
        if self.sfreq.val < 1: 
            self.sfreq.val = 1
        #self.sfreq.set_val(self.slide)
        self.freq = self.sfreq.val
        #self.ptitle = np.str(np.round(self.wave[self.slide],2))+' $\AA$'  
        self.l.set_data(self.image_data[int(self.freq)])
        #self.l3.set_data([self.wave[self.slide], self.wave[self.slide]],[-10,10])
        #ptitle = np.str(np.round(self.wave[self.slide],2))+' $\AA$'    
        self.ax1.set_title(str(self.wave[int(self.freq)])+' $\AA$')
        self.fig.canvas.draw_idle()
            #print("CAN'T GO LOWER!!")

    ##### Functions for input boxes

    def input_cmin(self, val):
        self.c_min = np.float(val)
        #_cmax = _cmin *100.
        self.l.set_clim([self.c_min, self.c_max])
        #plt.draw()

    def input_cmax(self, val):
        self.c_max = np.float(val)
        #_cmax = _cmin *100.
        self.l.set_clim([self.c_min, self.c_max])
        #plt.draw()

    ##### Functions for setting the colour map and the colour range

    """" 
    These functions allow to modify: 

    -) the colour map choosing among several options (more o different options can be add)
    -) the colour range choosing interactively the vmin and vmax values

    """
    def colorfunc(self, label):
        #self.cmap0=label
        self.l.set_cmap(label)
        #plt.draw()
        self.fig.canvas.draw_idle()

    def update_cmin(self, val, s=None):
        #f_min = s_fmin.val
        #c_cmin = self.s_cmin.val
        self.l.set_clim(self.s_cmin.val, self.s_cmax.val)
        self.text_box1.set_val("{:.3e}".format(self.s_cmin.val))
        plt.draw()

    def update_cmax(self, val, s=None):
        #f_max = s_fmax.val
        self.l.set_clim(self.s_cmin.val, self.s_cmax.val)
        self.text_box2.set_val("{:.3e}".format(self.s_cmax.val))
        #plt.draw()

    ##### Functions for setting the colour scale

    def scalefunc(self, label):
        if label =="log":
            self.norm=colors.LogNorm(vmin=self.c_min, vmax=self.c_max)
        elif label =="linear":
            self.norm =colors.Normalize(vmin=self.c_min, vmax=self.c_max)
        elif label =='power 2':   
            self.norm =colors.PowerNorm(2, vmin=self.c_min, vmax=self.c_max)       
        self.ax1.cla()
        self.cbar.remove()
        self.l = self.ax1.imshow(self.image_data[int(self.slide)], origin='lower',cmap=self.cmap0, norm=self.norm)
        self.ax1.set_xlabel('spaxel_x ')
        self.ax1.set_ylabel('spaxel_y')
        self.ax1.set_title(str(self.wave[int(self.slide)])+' $\AA$')
        #l5, = ax1.plot(x0, y0, '+', color='white', ms=10.5, mew=1.1)
        self.l4 = self.ax1.plot(self.x0, self.y0, '+', color='black')
        self.cbar = self.fig.colorbar(self.l, orientation="vertical", pad=0.02, ax=self.ax1)
        self.cbar.set_label('Flux [ erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')
        #radio.set_active(self.color_list.index(self.cmap0))
        #plt.draw()
        self.fig.canvas.draw_idle()


    ##### Functions for choosing among cube / integrated / S/N maps 

    def select_map(self, label):
        self.ax1.cla()
        self.cbar.remove()
        if label =="slide":
            self.l = self.ax1.imshow(self.image_data[int(self.slide)], origin='lower',cmap=self.cmap0, norm=self.norm)
        elif label =="sum":
            self.l = self.ax1.imshow(self.sum_cube, vmin=self.c_min_sum, vmax= self.c_max_sum, origin='lower',cmap=self.cmap0)
        elif label =="snr":
            self.l = self.ax1.imshow(self.snr, origin='lower', vmin=self.c_min_snr, vmax= self.c_max_snr, cmap=self.cmap0)

        self.ax1.set_xlabel('spaxel_x ')
        self.ax1.set_ylabel('spaxel_y')
        self.ax1.set_title(str(self.wave[int(self.slide)])+' $\AA$')
        #l5, = ax1.plot(x0, y0, '+', color='white', ms=10.5, mew=1.1)
        self.l4 = self.ax1.plot(self.x0, self.y0, '+', color='black')
        self.cbar = self.fig.colorbar(self.l, orientation="vertical", pad=0.02, ax=self.ax1)
        self.cbar.set_label('Flux [ erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')
        #radio.set_active(self.color_list.index(self.cmap0))
        #plt.draw()
        self.fig.canvas.draw_idle()
        
    ##### Function for interactive changing spaxel cliking in image or changing slide clicking in spectrum

    def onclick(self,event):
        # Check if the click was in ax1 and update spectrum
        if event.inaxes in [self.ax1]: 
            if int(event.xdata) < self.naxis1 and int(event.ydata) < self.naxis2:
                self.indexx = int(event.xdata)
                self.indexy = int(event.ydata)
                self.spectrum = self.image_data[:,self.indexy, self.indexx]
                self.l2.set_data(self.wave[:], self.spectrum)
                self.ax2.set_title('spaxel'+' '+str(self.indexx)+','+str(self.indexy))
                self.w_low_index=np.abs(np.array(self.wave)-self.l_min).argmin()
                self.w_high_index=np.abs(np.array(self.wave)-self.l_max).argmin()
                #self.f_max = np.nanmax(self.image_data[self.w_low_index: self.w_high_index, self.indexy, self.indexx])
                #self.ax2.set_ylim(self.f_min, self.f_max)
                self.fig.canvas.draw_idle()
                self.l4.set_data(self.indexx, self.indexy)
                self.l2.set_data(self.wave[:], self.image_data[:, self.indexy,self.indexx])
                #plt.draw()
                self.fig.canvas.draw_idle()
        # Check if the click was in ax2 and update slide (image)
        elif  event.inaxes in [self.ax2]:
            #print("\n\nClick event")
            self.val=np.abs(np.array(self.wave)-event.xdata).argmin()
            self.slide = int(self.val)
            self.sfreq.set_val(self.slide)

    ##### Functions for setting the flux range

    def update_fmin(self, val, s=None):
        #f_min = s_fmin.val
        self.ax2.set_ylim(self.s_fmin.val, self.s_fmax.val)
        #plt.draw()

    def update_fmax(self, val, s=None):
        #f_max = s_fmax.val
        self.ax2.set_ylim(self.s_fmin.val, self.s_fmax.val)
        #plt.draw()

    ##### Functions for setting the wavelength range

    """" 
    These functions allow to modify the wavelength range cliking on the bars to select any value

    """

    def update_lmin(self, val, s=None):
        #l_min = s_lmin.val
        self.ax2.set_xlim(self.s_lmin.val,self.s_lmax.val)
        #plt.draw()

    def update_lmax(self, val, s=None):
        #l_max = s_lmax.val
        self.ax2.set_xlim(self.l_min.val,self.s_lmax.val)
        #plt.draw()

    ##### Functions for reset parameters

    def reset(self, event):
        self.sfreq.reset()
        #self.cmap0.reset()
        #self.samp.reset()
        self.s_lmin.reset() 
        self.s_lmax.reset()
        self.s_cmin.reset() 
        self.s_cmax.reset()
        self.s_fmin.reset()
        self.s_fmax.reset()

    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------

    ##### Functions for opening specific datacubes

    """" 
    These functions allow to open specific datacubes for different large surveys:

    +) CALIFA
    +) MaNGA
    +) SAMI
    +) MUSE DATA
    +) CAVITY (ongoing survey)
    +) KOALA (ongoing survey)
    """

    def open_CALIFA(self, cube_fits_file, snr_lmin=None, snr_lmax=None):

        self.hdu_list = fits.open(cube_fits_file) 

        self.image_data = self.hdu_list[0].data
        self.image_data[self.image_data == 0] = np.nan # convert 0 to nan
        self.naxis1 = self.hdu_list[0].header['NAXIS1']
        self.naxis2 = self.hdu_list[0].header['NAXIS2']
        self.naxis3 = self.hdu_list[0].header['NAXIS3']
        self.crval3 = self.hdu_list[0].header['CRVAL3']
        self.cdelt3 = self.hdu_list[0].header['CDELT3']

        # Sum datacube
        self.sum_cube = np.nansum(self.image_data, axis=0)
        self.sum_cube[self.sum_cube == 0] = np.nan

        self.c_max_sum = np.nanmax(np.nansum(self.image_data, axis=0))
        #self.c_min = np.nanmin(np.nansum(self.image_data, axis=0))
        self.c_min_sum = np.nanpercentile(self.sum_cube, 50)

        # Wavelength:
        self.wave2 = np.arange(self.crval3, self.crval3+self.cdelt3*self.naxis3, self.cdelt3) 

        # Read the error spectra if available. Otherwise estimate the errors with the der_snr algorithm
        if len(self.hdu_list) > 2:
            self.image_error  = self.hdu_list[1].data
        elif len(self.hdu_list) <= 2:
            print("No error extension found. Estimating the error spectra with the der_snr algorithm")
            self.image_error = np.zeros( self.image_data.shape )
            for i in range( 0, self.image_data.shape[0] ):
                self.image_error[:,i] = der_snr( self.image_data[:,i] )

        if snr_lmin != None:
            self.snr_lmin = snr_lmin
        else:
            self.snr_lmin = np.nanmin(self.wave2)

        if snr_lmax != None:
            self.snr_lmax = snr_lmax
        else:
            self.snr_lmax = np.nanmax(self.wave2)

        self.idx_snr = np.where(np.logical_and( self.wave2 >= self.snr_lmin, self.wave2 <= self.snr_lmax ) )[0]
        self.signal = np.nanmedian(self.image_data[self.idx_snr,:],axis=0)
        self.noise  = np.abs(np.nanmedian(np.sqrt(self.image_error[self.idx_snr,:]),axis=0))
        self.snr    = self.signal / self.noise
        
        self.c_max_snr = np.nanmax(self.snr)
        #self.c_min = np.nanmin(np.nansum(self.image_data, axis=0))
        self.c_min_snr = np.nanpercentile(self.snr, 60)
        #print(self.c_min_snr, self.c_max_snr)

    
    def open_MUSE(self, cube_fits_file, snr_lmin=None, snr_lmax=None):

        self.hdu_list = fits.open(cube_fits_file) 

        self.image_data = self.hdu_list[1].data
        self.image_data[self.image_data == 1] = np.nan # convert 0 to nan
        self.naxis1 = self.hdu_list[1].header['NAXIS1']
        self.naxis2 = self.hdu_list[1].header['NAXIS2']
        self.naxis3 = self.hdu_list[1].header['NAXIS3']
        self.crval3 = self.hdu_list[1].header['CRVAL3']
        self.cdelt3 = self.hdu_list[1].header['CD3_3']

        # Sum datacube
        self.sum_cube = np.nansum(self.image_data, axis=0)
        self.sum_cube[self.sum_cube == 0] = np.nan

        self.c_max_sum = np.nanmax(np.nansum(self.image_data, axis=0))
        #self.c_min = np.nanmin(np.nansum(self.image_data, axis=0))
        self.c_min_sum = np.nanpercentile(self.sum_cube, 50)

        # Read the error spectra if available. Otherwise estimate the errors with the der_snr algorithm
        if len(self.hdu_list) == 3:
            self.image_error  = self.hdu_list[2].data
        elif len(self.hdu_list) == 2:
            print("No error extension found. Estimating the error spectra with the der_snr algorithm")
            self.image_error = np.zeros( self.image_data.shape )
            for i in range( 0, self.image_data.shape[1] ):
                self.image_error[:,i] = der_snr( self.image_data[:,i] )

        # Wavelength:
        self.wave2 = np.arange(self.crval3, self.crval3+self.cdelt3*self.naxis3, self.cdelt3) 

        if snr_lmin != None:
            self.snr_lmin = snr_lmin
        else:
            self.snr_lmin = np.nanmin(self.wave2)

        if snr_lmax != None:
            self.snr_lmax = snr_lmax
        else:
            self.snr_lmax = np.nanmax(self.wave2)

        self.idx_snr = np.where(np.logical_and( self.wave2 >= self.snr_lmin, self.wave2 <= self.snr_lmax ) )[0]
            #print(self.idx_snr)
            #self.idx_snr = np.nan_to_num(self.idx_snr)
        self.signal = np.nanmedian(self.image_data[self.idx_snr,:],axis=0)
        self.noise  = np.abs(np.nanmedian(np.sqrt(self.image_error[self.idx_snr,:]),axis=0))
        self.snr    = self.signal / self.noise
        
        self.c_max_snr = np.nanmax(self.snr)
            #self.c_min = np.nanmin(np.nansum(self.image_data, axis=0))
        self.c_min_snr = np.nanpercentile(self.snr, 60)
        print(self.c_min_snr, self.c_max_snr)

    def open_MANGA(self, cube_fits_file, snr_lmin=None, snr_lmax=None):

        self.hdu_list = fits.open(cube_fits_file) 

        self.image_data = self.hdu_list['FLUX'].data
        self.image_data[self.image_data == 1] = np.nan # convert 0 to nan
        self.naxis1 = self.hdu_list[1].header['NAXIS1']
        self.naxis2 = self.hdu_list[1].header['NAXIS2']
        self.naxis3 = self.hdu_list[1].header['NAXIS3']
        self.crval3 = self.hdu_list[1].header['CRVAL3']
        self.cdelt3 = self.hdu_list[1].header['CD3_3']

        # Sum datacube
        self.sum_cube = np.nansum(self.image_data, axis=0)
        self.sum_cube[self.sum_cube == 0] = np.nan

        self.c_max_sum = np.nanmax(np.nansum(self.image_data, axis=0))
        #self.c_min = np.nanmin(np.nansum(self.image_data, axis=0))
        self.c_min_sum = np.nanpercentile(self.sum_cube, 50)

        #self.image_error  = self.hdu_list['ERR'].data

        self.image_error = np.zeros( self.image_data.shape )
        for i in range(0, self.image_data.shape[1]):
            self.image_error[:,i] = der_snr( self.image_data[:,i] )


        # Wavelength:
        self.wave2 = self.hdu_list['WAVE'].data 

        if snr_lmin != None:
            self.snr_lmin = snr_lmin
        else:
            self.snr_lmin = np.nanmin(self.wave2)

        if snr_lmax != None:
            self.snr_lmax = snr_lmax
        else:
            self.snr_lmax = np.nanmax(self.wave2)

        self.idx_snr = np.where(np.logical_and( self.wave2 >= self.snr_lmin, self.wave2 <= self.snr_lmax ) )[0]
            #print(self.idx_snr)
            #self.idx_snr = np.nan_to_num(self.idx_snr)
        self.signal = np.nanmedian(self.image_data[self.idx_snr,:],axis=0)
        self.noise  = np.abs(np.nanmedian(np.sqrt(self.image_error[self.idx_snr,:]),axis=0))
        self.snr    = self.signal / self.noise
        
        self.c_max_snr = np.nanmax(self.snr)
            #self.c_min = np.nanmin(np.nansum(self.image_data, axis=0))
        self.c_min_snr = np.nanpercentile(self.snr, 60)
        #print(self.c_min_snr, self.c_max_snr)

    def open_SAMI(self, cube_fits_file, snr_lmin=None, snr_lmax=None):

        self.hdu_list = fits.open(cube_fits_file) 

        self.image_data = self.hdu_list[0].data
        self.image_data[self.image_data == 1] = np.nan # convert 0 to nan
        self.naxis1 = self.hdu_list[0].header['NAXIS1']
        self.naxis2 = self.hdu_list[0].header['NAXIS2']
        self.naxis3 = self.hdu_list[0].header['NAXIS3']
        self.crval3 = self.hdu_list[0].header['CRVAL3']
        self.cdelt3 = self.hdu_list[0].header['CDELT3']

        # Sum datacube
        self.sum_cube = np.nansum(self.image_data, axis=0)
        self.sum_cube[self.sum_cube == 0] = np.nan

        self.c_max_sum = np.nanmax(np.nansum(self.image_data, axis=0))
        #self.c_min = np.nanmin(np.nansum(self.image_data, axis=0))
        self.c_min_sum = np.nanpercentile(self.sum_cube, 50)

        # Read the error spectra if available. Otherwise estimate the errors with the der_snr algorithm
        if len(self.hdu_list) > 2:
            self.image_error  = self.hdu_list[1].data
        elif len(self.hdu_list) <= 2:
            print("No error extension found. Estimating the error spectra with the der_snr algorithm")
            self.image_error = np.zeros( self.image_data.shape )
            for i in range( 0, self.image_data.shape[0] ):
                self.image_error[:,i] = der_snr( self.image_data[:,i] )

        # Wavelength:
        self.wave2 = np.arange(self.crval3, self.crval3+self.cdelt3*self.naxis3, self.cdelt3)

        if snr_lmin != None:
            self.snr_lmin = snr_lmin
        else:
            self.snr_lmin = np.nanmin(self.wave2)

        if snr_lmax != None:
            self.snr_lmax = snr_lmax
        else:
            self.snr_lmax = np.nanmax(self.wave2)

        self.idx_snr = np.where(np.logical_and( self.wave2 >= self.snr_lmin, self.wave2 <= self.snr_lmax ) )[0]
            #print(self.idx_snr)
            #self.idx_snr = np.nan_to_num(self.idx_snr)
        self.signal = np.nanmedian(self.image_data[self.idx_snr,:],axis=0)
        self.noise  = np.abs(np.nanmedian(np.sqrt(self.image_error[self.idx_snr,:]),axis=0))
        self.snr    = self.signal / self.noise
        
        self.c_max_snr = np.nanmax(self.snr)
            #self.c_min = np.nanmin(np.nansum(self.image_data, axis=0))
        self.c_min_snr = np.nanpercentile(self.snr, 60)
        #print(self.c_min_snr, self.c_max_snr)

    def open_CAVITY(self, cube_fits_file, snr_lmin=None, snr_lmax=None):

        self.hdu_list = fits.open(cube_fits_file) 

        self.image_data = self.hdu_list[0].data
        self.image_data[self.image_data == 1] = np.nan # convert 0 to nan
        self.naxis1 = self.hdu_list[0].header['NAXIS1']
        self.naxis2 = self.hdu_list[0].header['NAXIS2']
        self.naxis3 = self.hdu_list[0].header['NAXIS3']
        self.crval3 = self.hdu_list[0].header['CRVAL3']
        self.cdelt3 = self.hdu_list[0].header['CD3_3']

        # Sum datacube
        self.sum_cube = np.nansum(self.image_data, axis=0)
        self.sum_cube[self.sum_cube == 0] = np.nan

        self.c_max_sum = np.nanmax(np.nansum(self.image_data, axis=0))
        #self.c_min = np.nanmin(np.nansum(self.image_data, axis=0))
        self.c_min_sum = np.nanpercentile(self.sum_cube, 50)
        
        # Read the error spectra if available. Otherwise estimate the errors with the der_snr algorithm
        if len(self.hdu_list) > 2:
            self.image_error  = self.hdu_list[1].data
        elif len(self.hdu_list) <= 2:
            print("No error extension found. Estimating the error spectra with the der_snr algorithm")
            self.image_error = np.zeros( self.image_data.shape )
            for i in range( 0, self.image_data.shape[0] ):
                self.image_error[:,i] = der_snr( self.image_data[:,i] )

        # Wavelength:
        self.wave2 = np.arange(self.crval3, self.crval3+self.cdelt3*self.naxis3, self.cdelt3)

        if snr_lmin != None:
            self.snr_lmin = snr_lmin
        else:
            self.snr_lmin = np.nanmin(self.wave2)

        if snr_lmax != None:
            self.snr_lmax = snr_lmax
        else:
            self.snr_lmax = np.nanmax(self.wave2)

        self.idx_snr = np.where(np.logical_and( self.wave2 >= self.snr_lmin, self.wave2 <= self.snr_lmax ) )[0]
            #print(self.idx_snr)
            #self.idx_snr = np.nan_to_num(self.idx_snr)
        self.signal = np.nanmedian(self.image_data[self.idx_snr,:],axis=0)
        self.noise  = np.abs(np.nanmedian(np.sqrt(self.image_error[self.idx_snr,:]),axis=0))
        self.snr    = self.signal / self.noise
        
        self.c_max_snr = np.nanmax(self.snr)
            #self.c_min = np.nanmin(np.nansum(self.image_data, axis=0))
        self.c_min_snr = np.nanpercentile(self.snr, 60)
        #print(self.c_min_snr, self.c_max_snr)

    def open_KOALA(self, cube_fits_file, snr_lmin=None, snr_lmax=None):

        self.hdu_list = fits.open(cube_fits_file) 

        self.image_data = self.hdu_list[1].data
        self.image_data[self.image_data == 1] = np.nan # convert 0 to nan
        self.naxis1 = self.hdu_list[1].header['NAXIS1']
        self.naxis2 = self.hdu_list[1].header['NAXIS2']
        self.naxis3 = self.hdu_list[1].header['NAXIS3']
        self.crval3 = self.hdu_list[1].header['CRVAL3']
        self.cdelt3 = self.hdu_list[1].header['CD3_3']

        # Sum datacube
        self.sum_cube = np.nansum(self.image_data, axis=0)
        self.sum_cube[self.sum_cube == 0] = np.nan

        self.c_max_sum = np.nanmax(np.nansum(self.image_data, axis=0))
        #self.c_min = np.nanmin(np.nansum(self.image_data, axis=0))
        self.c_min_sum = np.nanpercentile(self.sum_cube, 50)

        # Original:
        self.wave2 = np.arange(self.crval3, self.crval3+self.cdelt3*self.naxis3, self.cdelt3) #wavelength

        if snr_lmin != None:
            self.snr_lmin = snr_lmin
        else:
            self.snr_lmin = np.nanmin(self.wave2)

        if snr_lmax != None:
            self.snr_lmax = snr_lmax
        else:
            self.snr_lmax = np.nanmax(self.wave2)

        self.idx_snr = np.where(np.logical_and( self.wave2 >= self.snr_lmin, self.wave2 <= self.snr_lmax ) )[0]
            #print(self.idx_snr)
            #self.idx_snr = np.nan_to_num(self.idx_snr)
        self.signal = np.nanmedian(self.image_data[self.idx_snr,:],axis=0)
        self.noise  = np.abs(np.nanmedian(np.sqrt(self.image_error[self.idx_snr,:]),axis=0))
        self.snr    = self.signal / self.noise
        
        self.c_max_snr = np.nanmax(self.snr)
            #self.c_min = np.nanmin(np.nansum(self.image_data, axis=0))
        self.c_min_snr = np.nanpercentile(self.snr, 60)
        #print(self.c_min_snr, self.c_max_snr)

    # =====================================================================================
    # =====================================================================================
    # =====================================================================================
    