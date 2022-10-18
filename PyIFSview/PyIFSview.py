import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from astropy.io import fits
import matplotlib.gridspec as gridspec

# Angel update fuego color map
import matplotlib.colors as colors
fuego_color_map = colors.LinearSegmentedColormap.from_list("fuego", ((0.25, 0, 0),  (0.5,0,0),    (1, 0, 0), (1, 0.5, 0), (1, 0.75, 0), (1, 1, 0), (1, 1, 1)), N=256, gamma=1.0)
fuego_color_map.set_bad('lightgray')
plt.register_cmap(cmap=fuego_color_map)


class test(object):
    
    def __init__(self, cube_fits_file, slide=None, f_min =None, f_max =None, l_min =None, l_max =None, c_min =None, c_max =None, x0=None, y0=None):

        self.verbose= False


        if slide != None:
            self.slide = slide
        else:
            self.slide = 0
        

        if c_min != None:
        #self.f_min = 0.
            self.c_min = c_min
        else:
            self.c_min = -1.

        if c_max != None:
        #self.f_max = 0.
            self.c_max = c_max
        else:
            self.c_max = 0.

        if f_min != None:
            self.f_min = f_min
            print(f_min)
        else:
            self.f_min = 0.

        if f_max != None:
            self.f_max = f_max
            print(f_max)
        else:
            self.f_max = 0.
            print(self.f_max)


        if l_min != None:
        #self.f_min = 0.
            self.l_min = l_min
        else:
            self.l_min = 0.

        if l_max != None:
        #self.f_max = 0.
            self.l_max = l_max
        else:
            self.l_max = 0.

        if x0 != None:
        #self.f_max = 0.
            self.x0 = x0
        else:
            self.x0 = 0.

        if y0 != None:
        #self.f_max = 0.
            self.y0 = y0
        else:
            self.y0 = 0.


        self.cmap0 =""
        self.color_list = ['jet', 'gist_gray','viridis', 'gnuplot', 'gnuplot2', 'cubehelix', 'nipy_spectral', 'RdBu', 'fuego']
        

        #cube_fits_file = 'NGC2253.fits.gz'

        hdu_list = fits.open(cube_fits_file) 
        
        try:
            self.image_data = hdu_list[0].data
            self.image_data[self.image_data == 0] = np.nan # convert 0 to nan
            self.naxis1 = hdu_list[0].header['NAXIS1']
            self.naxis2 = hdu_list[0].header['NAXIS2']
            self.naxis3 = hdu_list[0].header['NAXIS3']
            self.crval3 = hdu_list[0].header['CRVAL3']
            self.cdelt3 = hdu_list[0].header['CDELT3']
            
        except:
            self.image_data = hdu_list[1].data
            self.image_data[self.image_data == 0] = np.nan # convert 0 to nan
            self.naxis1 = hdu_list[1].header['NAXIS1']
            self.naxis2 = hdu_list[1].header['NAXIS2']
            self.naxis3 = hdu_list[1].header['NAXIS3']
            self.crval3 = hdu_list[1].header['CRVAL3']
            self.cdelt3 = hdu_list[1].header['CD3_3'] 


        # Check we don't have negative values (for plotting log)
        minimo = np.nanmin(self.image_data)
        if self.verbose: print("Min value found in this cube:", minimo)
        if minimo < 0: 
            self.image_data= self.image_data + np.abs(minimo) +1E-20
            if self.verbose: 
                print("----> As this value is negative, adopting as min value ",np.nanmin(self.image_data))
        
        # Original:
        self.wave = np.arange(self.crval3, self.crval3+self.cdelt3*self.naxis3, self.cdelt3) #wavelength


        #self.fig, ax = plt.subplots()
        self.fig = plt.figure(figsize=(16,9))
        #self.fig.canvas.mpl_disconnect(fig.canvas.manager.key_press_handler_id)  # This is for avoiding an annoying warning
        spec2 = gridspec.GridSpec(ncols=3, nrows=3)
        self.ax1 = self.fig.add_subplot(spec2[1:3, 0:2])   # Image bottom left
        self.ax2 = self.fig.add_subplot(spec2[0, 0:3])     # Spectrum top row
        spec2.update(hspace=0.45)



        #### figure 2: spectrum

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
        
        ### define first plot with imshow

        self.norm=colors.LogNorm(vmin=self.c_min, vmax=self.c_max)

        self.l = self.ax1.imshow(self.image_data[self.slide], origin='lower', norm=self.norm)
        self.cbar = self.fig.colorbar(self.l, orientation="vertical", pad=0.02, ax=self.ax1)
        self.cbar.set_label('Flux [ erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')

        self.l4, = self.ax1.plot(self.x0, self.y0, '+', color='black')
        
        ### define second plot with the spectrum


        self.ax2.set_title('spaxel'+' '+str(self.x0)+','+str(self.y0))
        spectrum = self.image_data[:,self.y0, self.x0]   ############################## CAREFUL: data reads y0, x0 !!!!!!!!!!!!!!!!!!!!!!
        self.l2, = self.ax2.plot(self.wave[:], spectrum, label='emision')
        self.ax2.set_xlim(self.l_min, self.l_max)
        self.ax2.set_xlabel('Wavelength [$\AA$]')
        self.ax2.minorticks_on()
        self.ax2.set_ylabel('Flux [ erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')

        ##### defining interactive Slide for wavelenght slide

        axcolor = 'lightgoldenrodyellow'
        axfreq = plt.axes([0.125, 0.92, 0.776, 0.05], facecolor=axcolor)

        self.sfreq = Slider(axfreq, '', 0, self.naxis3-2, valinit=self.slide)
        self.sfreq.valtext.set_visible(False)
        
        self.sfreq.on_changed(self.update)

        ##### defining interactive button for reset
        
        resetax = plt.axes([0.85, 0.2, 0.05, 0.05])
        self.button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
        self.button.on_clicked(self.reset)

        ##### defining interactive button for cmaps
   
        rax = plt.axes([0.65, 0.20, 0.15, 0.22])
        
        color_list = ['jet', 'gist_gray','viridis', 'gnuplot', 'gnuplot2', 'cubehelix', 'nipy_spectral', 'RdBu', 'fuego']
    
        self.radio = RadioButtons(rax, (color_list), active=color_list.index(self.cmap0))
        self.radio.on_clicked(self.colorfunc)


        ##### choosing linear / log map / power
        
        scaleax = plt.axes([0.85, 0.25, 0.05, 0.1])
        self.scale = RadioButtons(scaleax, ('linear', 'log', 'power 2'), active=1)
        self.scale.on_clicked(self.scalefunc)


        ##### defining interactive limits for cmaps
    
        ax_cmin = plt.axes([0.65, 0.1, 0.19, 0.03])
        ax_cmax  = plt.axes([0.65, 0.15, 0.19, 0.03])
        
        #s_cmin = Slider(ax_cmin, 'min', min_value, c_max, valinit=c_min, valfmt='%5.2E')
        #s_cmin = Sliderlog(ax_cmin, 'min', min_value, c_max, valinit=np.log10(c_min), valfmt='%5.2E')
        self.color_min =  np.nanpercentile(self.image_data, 1)    # These two must be the absolute min and max, otherwise this points couldn't be reached by slider!!!
        self.color_max =  np.nanmax(self.image_data)
        
        self.s_cmin = Slider(ax_cmin, 'vmin', self.color_min, self.color_max, valinit=self.c_min, valfmt='%5.2E')
        #s_cmin.valtext.set_visible(False)
        self.s_cmax = Slider(ax_cmax, 'vmax', self.color_min, self.color_max, valinit=self.c_max, valfmt='%5.2E')
        #s_cmax.valtext.set_visible(False)
        
        self.s_cmin.on_changed(self.update_cmin)
        self.s_cmax.on_changed(self.update_cmax)

         ##### defining interactive limits for flux in spectrum
    
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


        ##### defining interactive limits for wavelength in spectrum
    
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

        #### onclick
    
        cid = self.fig.canvas.mpl_connect('button_press_event', self.onclick)  

        plt.show()

    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------

    def update(self,val):
        #amp = self.samp.val
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
        self.l.set_clim(self.s_cmin.val, self.s_cmax.val)
        #plt.draw()

    def update_cmax(self, val, s=None):
        #f_max = s_fmax.val
        self.l.set_clim(self.s_cmin.val, self.s_cmax.val)
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
        
    ##### Function for interactive changing spaxel cliking in image or changing slide clicking in spectrum

    def onclick(self,event):
        #global f_min,f_max,slide, ax1
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
        self.ax2.set_xlim(self.l_min,self.s_lmax.val)
        #plt.draw()

    ##### Functions for reset parameters

    def reset(self, event):
        self.sfreq.reset()
        #self.cmap0.reset()
        #self.samp.reset()
        self.s_cmin.reset() 
        self.s_cmax.reset()
        self.s_fmin.reset()
        self.s_fmax.reset()

    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
