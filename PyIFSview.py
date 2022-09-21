'''
PyIFSview - View & interact with IFS datacubes

This code contains the visualisation tool developed by Ignacio del Moral-Castro (IAC, ignaciodelmoralcastro@gmail.com). 
The basic functionality of this interactive tool is to visualise IFU data (e.g. CALIFA, MANGA, SAMI ...)

Updated for KOALA by Angel López-Sánchez (AAO-MQ, Angel.Lopez-Sanchez@mq.edu.au)

Examples:

    - From command line:
    
> python PyIFSview.py /DATA/KOALA/Jamila/20180227/385R/He2_10Ar_385R_JAMILA.fits s=318 l_min=6500.2 x=43 y=35 l_max=6800. c_min=1E-18 c_max=1E-14 cmap="gist_gray" f_min=0. f_max=2.E-14

    - From Python:
        
> import os     (only needed once!)
> os.system("python PyIFSview.py /DATA/KOALA/Jamila/20180227/385R/He2_10Ar_385R_JAMILA.fits  s=318 l_min=6500. x=43 y=35 l_max=6800.")

'''
version="Version 1.0 - 1 Aug 2021"

########################################################################
#     importing libraries                                              #
########################################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, TextBox
from astropy.io import fits
import matplotlib.gridspec as gridspec
#from astropy.visualization import simple_norm
import sys

#import warnings
#warnings.simplefilter('error', UserWarning)

# Angel update fuego color map
import matplotlib.colors as colors
fuego_color_map = colors.LinearSegmentedColormap.from_list("fuego", ((0.25, 0, 0),  (0.5,0,0),    (1, 0, 0), (1, 0.5, 0), (1, 0.75, 0), (1, 1, 0), (1, 1, 1)), N=256, gamma=1.0)
fuego_color_map.set_bad('lightgray')
plt.register_cmap(cmap=fuego_color_map)

# Set global parameters

global l,l3,ax1
global cbar
global norm
global slide
global c_max, c_min
global l_max,l_min
global f_min,f_max
global x0,y0
global verbose

verbose= False
slide = 0
c_max = 0.
c_min = -1.
f_min = 0.
f_max =0.
l_max =0.
l_min =0.
x0=0
y0=0
cmap0 =""
color_list = ['jet', 'gist_gray','viridis', 'gnuplot', 'gnuplot2', 'cubehelix', 'nipy_spectral', 'RdBu', 'fuego']


# Log scale by default
#norm=colors.Normalize()
#norm=colors.LogNorm()
#norm=colors.PowerNorm(gamma=1/2.)


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# class Sliderlog(Slider):

#     """
#     Logarithmic slider.

#     Takes in every method and function of the matplotlib's slider.

#     Set slider to *val* visually so the slider still is lineat but display 10**val next to the slider.

#     Return 10**val to the update function (func)
#     """

#     def set_val(self, val):

#         xy = self.poly.xy
#         if self.orientation == 'vertical':
#             xy[1] = 0, val
#             xy[2] = 1, val
#         else:
#             xy[2] = val, 1
#             xy[3] = val, 0
#         self.poly.xy = xy
#         self.valtext.set_text(self.valfmt % 10**val)   # Modified to display 10**val instead of val
#         if self.drawon:
#             self.ax.figure.canvas.draw_idle()
#         self.val = val
#         if not self.eventson:
#             return
#         for cid, func in self.observers.items():
#                 func(10**val)
                
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

##### Function for changing the colour scale linear/log/power


def scalefunc(label):
    global l
    global cbar
    global radio
    if label =="log":
        norm=colors.LogNorm(vmin=c_min, vmax=c_max)
    elif label =="linear":
        norm =colors.Normalize(vmin=c_min, vmax=c_max)
    elif label =='power 2':   
        norm =colors.PowerNorm(2, vmin=c_min, vmax=c_max)       
    ax1.cla()
    cbar.remove()
    l = ax1.imshow(image_data[slide], origin='lower',cmap=cmap0, norm=norm)
    l5, = ax1.plot(x0, y0, '+', color='white', ms=10.5, mew=1.1)
    l4, = ax1.plot(x0, y0, '+', color='black')
    cbar = fig.colorbar(l, orientation="vertical", pad=0.02, ax=ax1)
    cbar.set_label('Flux [ erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')
    radio.set_active(color_list.index(cmap0))
    #plt.draw()

##### Functions for update value of the slide

def update_slide(val):
    slide = np.int(np.round(val))
    text_box10.set_val(slide)       # Updating the text box will do the rest!!

def input_slide(val):
    global slide
    slide = np.int(np.float(val))
    sslide.set_val(slide)
    ptitle = np.str(np.round(wave[slide],2))+' $\AA$'  
    l.set_data(image_data[slide])
    l3.set_data([wave[slide], wave[slide]],[-10,10])
    ptitle = np.str(np.round(wave[slide],2))+' $\AA$'    
    ax1.set_title(ptitle)
    #fig.canvas.draw_idle()   
    if verbose: print("Update slide:  ",slide,"   at   ",ptitle)
    
def planeup_event(event):
    global slide
    slide = slide+1
    if slide > naxis3:
        slide = naxis3
        #print("CAN'T GO HIGHER!!")
    text_box10.set_val(slide)       # Updating the the text box will do the rest!!

def planedown_event(event):
    global slide
    slide = slide-1
    if slide < 1: 
        slide = 1
        #print("CAN'T GO LOWER!!")
    text_box10.set_val(slide)       # Updating the the text box will do the rest!!


##### Functions for setting the colour map and the colour range

def colorfunc(label):
    global cmap0,x0,y0
    cmap0=label
    l.set_cmap(cmap0)
    #plt.draw()
    
def update_cmin(val, s=None):
    global c_min,c_max
    c_min = s_cmin.val
    l.set_clim([c_min, c_max])
    text_box1.set_val("{:.3e}".format(c_min))
    #plt.draw()

def update_cmax(val, s=None):
    global c_min,c_max
    c_max = s_cmax.val
    l.set_clim([c_min, c_max])
    text_box2.set_val("{:.3e}".format(c_max))
    #plt.draw()

##### Functions for setting the wavelength range

def update_lmin(val, s=None):
    global l_min,l_max
    l_min = s_lmin.val
    ax2.set_xlim(l_min,l_max)
    text_box4.set_val(np.round(l_min,2))
    #plt.draw()

def update_lmax(val, s=None):
    global l_min,l_max
    l_max = s_lmax.val
    ax2.set_xlim(l_min,l_max)
    text_box3.set_val(np.round(l_max,2))
    #plt.draw()

##### Functions for setting the flux range

def update_fmin(val, s=None):
    global f_min,f_max
    f_min = s_fmin.val
    ax2.set_ylim(f_min,f_max)
    text_box5.set_val("{:5.2e}".format(f_min))
    #plt.draw()

def update_fmax(val, s=None):
    global f_min,f_max
    f_max = s_fmax.val
    ax2.set_ylim(f_min,f_max)
    text_box6.set_val("{:5.2e}".format(f_max))
    #plt.draw()

##### Functions for input boxes

def input_cmin(val):
    global c_min
    c_min = np.float(val)
    #_cmax = _cmin *100.
    l.set_clim([c_min, c_max])
    #plt.draw()

def input_cmax(val):
    global c_max
    c_max = np.float(val)
    #_cmin = _cmax /100.
    l.set_clim([c_min, c_max])
    #plt.draw()    
    
def input_lmin(val):
    global l_max,l_min
    l_min = np.float(val)
    s_lmin.set_val(l_min)
    ax2.set_xlim(l_min,l_max)
    #plt.draw()    

def input_lmax(val):
    global l_max,l_min
    l_max = np.float(val)
    s_lmax.set_val(l_max)
    ax2.set_xlim(l_min,l_max)
    #plt.draw()

def input_fmin(val):
    global f_max,f_min
    f_min = np.float(val)
    s_fmin.set_val(f_min)
    ax2.set_ylim(f_min,f_max)
    s_fmin.set_val(f_min)
    #plt.draw()

def input_fmax(val):
    global f_max,f_min
    f_max = np.float(val)
    s_fmax.set_val(f_max)
    ax2.set_ylim(f_min,f_max)
    #plt.draw()

def input_spaxel(x,y):
    global f_min, f_max
    l4.set_data(x,y)
    l5.set_data(x,y)
    spectrum = image_data[:,y, x]
    l2.set_data(wave[:], spectrum)
    ax2.set_title('spaxel'+' '+str(x)+','+str(y))
    w_low_index=np.abs(np.array(wave)-l_min).argmin()
    w_high_index=np.abs(np.array(wave)-l_max).argmin()
    f_max = np.nanmax(image_data[w_low_index:w_high_index,y,x])
    ax2.set_ylim(f_min,f_max)
    text_box6.set_val("{:5.2e}".format(f_max))   
    text_box7.set_val(x)
    text_box8.set_val(y)
    #plt.draw()
    #if verbose: print("Update to spaxel: ",x,y)
    
def input_x0(val):
    global y0,x0
    x0=np.int(val)
    input_spaxel(x0,y0)
    #if verbose: print("input_x0: ",x0)
def input_y0(val):
    global x0,y0
    y0=np.int(val)
    input_spaxel(x0,y0)        
    #if verbose: print("input_y0: ",y0)


##### Function for interactive changing spaxel cliking in image or changing slide clicking in spectrum

def onclick(event):
    global f_min,f_max,slide
    # Check if the click was in ax1 and update spectrum
    if event.inaxes in [ax1]: 
        if int(event.xdata) < naxis1 and int(event.ydata) < naxis2:
            indexx = int(np.round(event.xdata))
            indexy = int(np.round(event.ydata))
            input_spaxel(indexx,indexy)
            #text_box7.set_val(indexx)
            #text_box8.set_val(indexy)
    # Check if the click was in ax2 and update slide (image)
    elif  event.inaxes in [ax2]:
        #print("\n\nClick event")
        val=np.abs(np.array(wave)-event.xdata).argmin()
        slide = np.int(val)
        sslide.set_val(slide)

##### Function for reset

def reset(event):
    """Reset the interactive plots to inital values: slide and vmin and vmax."""
    sslide.reset()
    s_cmin.reset() 
    s_cmax.reset()
    s_fmin.reset()
    s_fmax.reset()
    s_lmin.reset()
    s_lmax.reset()
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------


if __name__ == "__main__": 

    # Read cube_fits_file
    
    #cube_fits_file = '/DATA/KOALA/Python/PACE/385R/POX4_385R_20201029.fits' 
    #cube_fits_file = "/DATA/KOALA/Jamila/20180227/385R/He2_10Ar_385R_JAMILA.fits"
    cube_fits_file = sys.argv[1]
    
    
    # Read arguments, 
    # ANGEL: I know this is probably not the optimal way of doing it but it works...
    
    argv = sys.argv[2:]
    for value in argv:
        if value[:6] == "l_min=":
            l_min=np.float(value[6:])
        elif value[:6] == "l_max=":
            l_max=np.float(value[6:])
        elif value[:6] == "slide=":
            slide=np.int(value[6:])
        elif value[:2] == "s=":
            slide=np.int(value[2:])
        elif value[:6] == "c_min=":
            c_min=np.float(value[6:])
        elif value[:6] == "c_max=":
            c_max=np.float(value[6:])    
        elif value[:6] == "f_min=":
            f_min=np.float(value[6:])
        elif value[:6] == "f_max=":
            f_max=np.float(value[6:])  
        elif value[:2] == "x=":
            x0=np.int(value[2:])
        elif value[:2] == "y=":
            y0=np.int(value[2:])
        elif value[:5] == "cmap=":
            cmap0=value[5:]
    

#   IF TRYING THIS EVERYTHING BREAKS...

#     PyIFSview(cube_fits_file, slide = slide, c_max = c_max, c_min = c_min, f_min = f_min, f_max =f_max, 
#                   l_max =l_max, l_min =l_min, x0=x0, y0=y0, cmap0 =cmap0) 

# def PyIFSview(cube_fits_file, slide = 0, c_max = 0., c_min = -1. ,f_min = -1., f_max =0.,
#               l_max =0., l_min =0., x0=0, y0=0, cmap0 =""):
    
    ########################################################################
    #     opening data and defining wavelength range                       #
    ########################################################################
    
    hdu_list = fits.open(cube_fits_file) 
    
    try:
      image_data = hdu_list[0].data
      image_data[image_data == 0] = np.nan # convert 0 to nan
      naxis1 = hdu_list[0].header['NAXIS1']
      naxis2 = hdu_list[0].header['NAXIS2']
      naxis3 = hdu_list[0].header['NAXIS3']
      crval3 = hdu_list[0].header['CRVAL3']
      cdelt3 = hdu_list[0].header['CDELT3']
    
    except:
      image_data = hdu_list[1].data
      image_data[image_data == 0] = np.nan # convert 0 to nan
      naxis1 = hdu_list[1].header['NAXIS1']
      naxis2 = hdu_list[1].header['NAXIS2']
      naxis3 = hdu_list[1].header['NAXIS3']
      crval3 = hdu_list[1].header['CRVAL3']
      cdelt3 = hdu_list[1].header['CD3_3'] 
    
    
    # Check we don't have negative values (for plotting log)
    minimo = np.nanmin(image_data)
    if verbose: print("Min value found in this cube:", minimo)
    if minimo < 0: 
        image_data= image_data + np.abs(minimo) +1E-20
        if verbose: print("----> As this value is negative, adopting as min value ",np.nanmin(image_data))
    
    # Original:
    #wave = np.arange(crval3, crval3+cdelt3*naxis3, cdelt3) #wavelength
    # PyKOALA:
    #wave = np.arange(crval3-cdelt3*naxis3/2., crval3+cdelt3*naxis3/2., cdelt3) #wavelength
    wave = np.arange(crval3-cdelt3*naxis3/2., crval3+cdelt3*(naxis3-1)/2. , cdelt3) #wavelength
      
    
    ########################################################################
    #                plotting figures                                      #
    ########################################################################
    
    fig = plt.figure(figsize=(16,9))
    fig.canvas.mpl_disconnect(fig.canvas.manager.key_press_handler_id)  # This is for avoiding an annoying warning
    spec2 = gridspec.GridSpec(ncols=3, nrows=3)
    ax1 = fig.add_subplot(spec2[1:3, 0:2])   # Image bottom left
    ax2 = fig.add_subplot(spec2[0, 0:3])     # Spectrum top row
    spec2.update(hspace=0.45)
    
    #### initial values
    
    if cmap0 == "" : cmap0 = 'jet'
    elif cmap0 not in color_list : cmap0 = 'jet'
    
    if slide == 0: slide = int(naxis3/2) 
    if c_min == -1: c_min = np.nanpercentile(image_data, 1)
    if c_max ==0 : c_max = np.nanmax(image_data)
    
    #### figure 1: MAP at slide with imshow
    
    if x0 ==0 : x0 = int(naxis1/2)
    if y0 ==0 : y0 = int(naxis2/2)
    
    ax1.set_title(str(np.round(wave[slide],2))+' $\AA$')
    
    norm=colors.LogNorm(vmin=c_min, vmax=c_max)
    
    l = ax1.imshow(image_data[slide], origin='lower', cmap=cmap0, norm=norm)
    
    l5, = ax1.plot(x0, y0, '+', color='white', ms=10.5, mew=1.1)
    l4, = ax1.plot(x0, y0, '+', color='black')
    
    ax1.set_xlim(0,naxis1)
    ax1.set_ylim(0,naxis2)
    ax1.axes.tick_params(axis='both',which='minor',direction='out',length=10)
    
    ax1.set_xlabel('spaxel_x ')
    ax1.set_ylabel('spaxel_y')
    
    cbar = fig.colorbar(l, orientation="vertical", pad=0.02, ax=ax1)
    cbar.set_label('Flux [ erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')
    
    
    #### figure 2: spectrum
    
    ax2.set_title('spaxel'+' '+str(x0)+','+str(y0))
    spectrum = image_data[:,y0, x0]   ############################## CAREFUL: KOALA data reads y0, x0 !!!!!!!!!!!!!!!!!!!!!!
    l2, = ax2.plot(wave[:], spectrum, label='emision')
    ax2.set_xlabel('Wavelength [$\AA$]')
    ax2.minorticks_on()
    ax2.set_ylabel('Flux [ erg cm$^{-2}$ s$^{-1}$ A$^{-1}$]')
    
    
    l3, = ax2.plot([wave[slide],wave[slide]],[-10,10]) # interactive vertical line
                      
    if l_min == 0 :l_min = min(wave) 
    if l_max == 0: l_max = max(wave)
    
    ax2.set_xlim(l_min, l_max)
            
    w_low_index=np.abs(np.array(wave)-l_min).argmin()
    w_high_index=np.abs(np.array(wave)-l_max).argmin()
    if f_max == 0:  f_max =np.nanmax(image_data[w_low_index:w_high_index,y0, x0])
    
    ax2.set_ylim(f_min,f_max)
    
    
    #### defining interactive Slide 
    
    axfreq = plt.axes([0.125, 0.92, 0.776, 0.05])
    sslide = Slider(axfreq, '', 0, naxis3-2, valinit=slide, valfmt="%i")
    sslide.valtext.set_visible(False)
      
    slide_values_ = np.arange(0, naxis3-2, len(wave)/8) 
    slide_values=[]
    for i in range(len(slide_values)):
    	slide_values.extend(np.int(slide_values[i]))
    wave_values = np.round(wave[slide_values],2 )
    
    axfreq.set_xticks(slide_values, minor = False)
    axfreq.set_xticklabels(wave_values)
    
    sslide.on_changed(update_slide)
      
    pslide = " {}".format(slide)
    axbox10 = plt.axes([0.925, 0.92, 0.03, 0.05])
    text_box10 = TextBox(axbox10, '', initial=pslide)
    text_box10.on_submit(input_slide)
    
    
    ##### Defining moving plane by plane 
    
    planedown = plt.axes([0.108, 0.92, 0.017, 0.05])
    planeup = plt.axes([0.901, 0.92, 0.017, 0.05])
    
    button_down = Button(planedown, '-', color='lightgoldenrodyellow', hovercolor='0.975')
    button_up = Button(planeup, '+', color='lightgoldenrodyellow', hovercolor='0.975')
        
    button_up.on_clicked(planeup_event)
    button_down.on_clicked(planedown_event)
    
    
    ##### defining interactive button for cmaps
    
    rax = plt.axes([0.65, 0.20, 0.15, 0.22])
    
    #radio = RadioButtons(rax, ('jet', 'gist_gray','viridis', 'gnuplot', 'gnuplot2', 'cubehelix', 'nipy_spectral', 'RdBu', 'fuego'), active=0)
    
    radio = RadioButtons(rax, (color_list), active=color_list.index(cmap0))
    radio.on_clicked(colorfunc)
    
    ##### defining interactive limits for wavelength in spectrum
    
    ax_lmin = plt.axes([0.65, 0.55, 0.19, 0.03])
    ax_lmax  = plt.axes([0.65, 0.6, 0.19, 0.03])
    
    s_lmin = Slider(ax_lmin, 'l_min', np.min(wave), np.max(wave), valinit=l_min, valfmt='%7.2f')
    s_lmin.valtext.set_visible(False)
    s_lmax = Slider(ax_lmax, 'l_max', np.min(wave), np.max(wave), valinit=l_max, valfmt='%7.2f')
    s_lmax.valtext.set_visible(False)
    
    s_lmin.on_changed(update_lmin)
    s_lmax.on_changed(update_lmax)
    
    ##### defining interactive limits for flux in spectrum
    
    if f_max == 0: f_max = np.nanmax(image_data)
    
    ax_fmin = plt.axes([0.65, 0.45, 0.19, 0.03])
    ax_fmax  = plt.axes([0.65, 0.5, 0.19, 0.03])

    
    f_min_real = np.nanmin([-np.nanmin(image_data),np.nanmin(image_data)])   #(np.abs(f_min) -f_max)/10. 
    f_max_real = np.nanmax(image_data)
        
    s_fmin = Slider(ax_fmin, 'f_min', f_min_real, f_max_real, valinit=f_min, valfmt='%5.2E')
    s_fmin.valtext.set_visible(False)
    s_fmax = Slider(ax_fmax, 'f_max', f_min_real, f_max_real, valinit=f_max, valfmt='%5.2E')
    s_fmax.valtext.set_visible(False)
    
    s_fmin.on_changed(update_fmin)
    s_fmax.on_changed(update_fmax)
    
    ##### defining interactive limits for cmaps
    
    ax_cmin = plt.axes([0.65, 0.1, 0.19, 0.03])
    ax_cmax  = plt.axes([0.65, 0.15, 0.19, 0.03])
    
    #s_cmin = Slider(ax_cmin, 'min', min_value, c_max, valinit=c_min, valfmt='%5.2E')
    #s_cmin = Sliderlog(ax_cmin, 'min', min_value, c_max, valinit=np.log10(c_min), valfmt='%5.2E')
    color_min =  np.nanpercentile(image_data, 1)    # These two must be the absolute min and max, otherwise this points couldn't be reached by slider!!!
    color_max =  np.nanmax(image_data)
    
    s_cmin = Slider(ax_cmin, 'c_min', color_min, color_max, valinit=c_min, valfmt='%5.2E')
    s_cmin.valtext.set_visible(False)
    s_cmax = Slider(ax_cmax, 'c_max', color_min, color_max, valinit=c_max, valfmt='%5.2E')
    s_cmax.valtext.set_visible(False)
    
    s_cmin.on_changed(update_cmin)
    s_cmax.on_changed(update_cmax)
    
    
    ########  Text boxes with values
        
    pcmin = "{:.3e}".format(c_min)    
    axbox1 = plt.axes([0.85, 0.1, 0.06, 0.03])
    text_box1 = TextBox(axbox1, '', initial=pcmin)
    text_box1.on_submit(input_cmin)
        
    pcmax = "{:.3e}".format(c_max)
    axbox2 = plt.axes([0.85, 0.15, 0.06, 0.03])
    text_box2 = TextBox(axbox2, '', initial=pcmax)
    text_box2.on_submit(input_cmax)
     
    plmin = "{:.2f}".format(l_min)
    axbox4 = plt.axes([0.85, 0.55, 0.06, 0.03])
    text_box4 = TextBox(axbox4, '', initial=plmin)
    text_box4.on_submit(input_lmin)
     
    plmax = "{:.2f}".format(l_max)
    axbox3 = plt.axes([0.85, 0.6, 0.06, 0.03])
    text_box3 = TextBox(axbox3, '', initial=plmax)
    text_box3.on_submit(input_lmax)
        
    pfmin = "{:.3e}".format(f_min)
    axbox5 = plt.axes([0.85, 0.45, 0.06, 0.03])
    text_box5 = TextBox(axbox5, '', initial=pfmin)
    text_box5.on_submit(input_fmin)
    
    pfmax = "{:.3e}".format(f_max)
    axbox6 = plt.axes([0.85, 0.5, 0.06, 0.03])
    text_box6 = TextBox(axbox6, '', initial=pfmax)
    text_box6.on_submit(input_fmax)
    
    
    px = "{:}".format(x0)
    axbox7 = plt.axes([0.85, 0.4, 0.04, 0.03])
    text_box7 = TextBox(axbox7, 'x =', initial=px)
    text_box7.on_submit(input_x0)
    
    py = "{:}".format(y0)
    axbox8 = plt.axes([0.85, 0.36, 0.04, 0.03])
    text_box8 = TextBox(axbox8, 'y =', initial=py)
    text_box8.on_submit(input_y0)
    
    ##### choosing linear / log map
    
    scaleax = plt.axes([0.85, 0.25, 0.05, 0.1])
    scale = RadioButtons(scaleax, ('linear', 'log', 'power 2'), active=1)
    scale.on_clicked(scalefunc)
    
    ####
    
    cid = fig.canvas.mpl_connect('button_press_event', onclick)    
    
    ### reset button
    
    resetax = plt.axes([0.85, 0.2, 0.05, 0.05])
    button = Button(resetax, 'Reset', color='lightgoldenrodyellow', hovercolor='0.975')
        
    button.on_clicked(reset)
    
    plt.gcf().canvas.set_window_title('PyIFSview in file: '+cube_fits_file)
    plt.show()