import spectrograph as spec
import telescope as tel
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits

import math
import sys

def observe_spectrum(instrument, texp, input_wv, input_spec, skyspec):#skyfile="eso_newmoon_radiance.txt"):

    #if (os.path.isfile(skyfile)):
    #    full_path = skyfile
    #else:
    #    full_path = os.environ['COATINGS_PATH']+skyfile
    
    # ESO dark night sky spectrum is in weird units: photons/m2/s/micron/arcsec2
    # interpolate immediately onto the instrument's wavelength grid
    #skyspec= np.genfromtxt(full_path,usecols=[0,1],names=['waves_nm','skyflux'])
    sky   = np.interp(instrument.waves,skyspec['waves_nm'],skyspec['skyflux'])

    # Object containing operational parameters of the telescope
    magellan = tel.Telescope()

    # Poisson noise is counting statistics, and the thing we count are e- in the CCD
    # that originated from sky photons.
    # "sky" is in photons/m2/s/micron/arcsec^2, we need to turn this into electrons (e-)
    # which requires multiplying by all factors in the demonimator 
    skyphotons = sky * \
        magellan.Atel/(100**2) * \
        texp * \
        (instrument.waves/1.0e3)/instrument.R * \
        instrument.fiber.Afib * \
        instrument.throughput * \
        magellan.throughput(instrument.waves) #e-

    # Convert from energy units (ergs) to counted photons
    h = 6.6e-27
    c = 3.0e17   # nm/sec
    input_photons = input_spec / (h * c/input_wv)

    objspec = np.interp(instrument.waves,input_wv,input_photons)

    # Same calculation as above, except assume all of the light goes down one fiber
    # So, you don't need to integrate over the fiber area, but do need everything else
    objphotons = objspec * \
        magellan.Atel * \
        texp * \
        (instrument.waves)/instrument.R * \
        instrument.throughput * \
        magellan.throughput(instrument.waves) #e-

    readnoise = instrument.sensor.rn #e-
    dark      = instrument.sensor.dark * texp #e-/s * s

    # This calculates the number of pixels that a fiber subtends. When we simulate the extracted the spectrum, we need
    # to quadrature sum read noise from ALL pixels in the profile, not just one pixel.
    pix_resel = (instrument.fiber.dFib) * (instrument.f_cam / instrument.f_col) / instrument.sensor.pixelsize
    
    skynoise = np.sqrt(skyphotons)
    totnoise = np.sqrt(skyphotons + dark*np.ceil(pix_resel) + readnoise**2*np.ceil(pix_resel))
    
    if any([math.isnan(n) for n in totnoise]):
        breakpoint()

    return objphotons, totnoise

#    plt.plot(instrument.waves, skynoise/totnoise)
#    plt.ylim([0,1])
#    plt.show()



def observe(instrument, texp, input_spectrum='lbg_shapley.fits', telescope='magellan', 
            skyfile="eso_newmoon_radiance.txt"):

    if (os.path.isfile(skyfile)):
        full_path = skyfile
    else:
        full_path = os.environ['COATINGS_PATH']+skyfile


    # Total noise is composed of
    # (1) Poisson noise from the sky, which typically dominates,
    # (2) Poisson noise from the object, which matters when the object is brighter than the sky
    # (3) Read noise from the sensor
    #
    # Since (1) is usually the main noise source, we get a night sky spectrum from the literature

    
    # ESO dark night sky spectrum is in weird units: photons/m2/s/micron/arcsec2
    # interpolate immediately onto the instrument's wavelength grid

    skyspec= np.genfromtxt(full_path,usecols=[0,1],names=['waves_nm','skyflux'])
    sky   = np.interp(instrument.waves,skyspec['waves_nm'],skyspec['skyflux'])

    # Object containing operational parameters of the telescope
    magellan = tel.Telescope()

    # Poisson noise is counting statistics, and the thing we count are e- in the CCD
    # that originated from sky photons.

    # "sky" is in photons/m2/s/micron/arcsec^2, we need to turn this into electrons (e-)
    # which requires multiplying by all factors in the demonimator
    # 
    # (a) magellan.Atel/(100**2) is the telescope area in m2
    # (b) texp is the number of seconds in the exposure
    # (c) instrument.fiber.Afib is the area of one fiber in arcsec^2 (bigger fibers would take in more sky)
    # (d) The tricky one is that "per micron," which is the bandwidth of one pixel. It's not the same as the
    #       dlambda of one pixel, becuase you get light spread from adjacent wavelengths from the wings of the
    #       spectrograph line spread function (LSF).  The spectral resolution R=lambda/dlambda where
    #       dlambda is one "resolution element" and is the effective bandpass of a pixel after convolution with
    #       the LSF. This means that dlambda = lambda / R, and the factoe of 1.0e3 converts from nm to microns.
    #       That's the meaning of the (instrument.waves/1.0e3)/instrument.R line.
    #
    # Finally, to get from sky photons per pixel interval to photoelectrons in the CCD, you need to account
    # for throughput losses in the telescope and instrument.
    
    skyphotons = sky * \
        magellan.Atel/(100**2) * \
        texp * \
        (instrument.waves/1.0e3)/instrument.R * \
        instrument.fiber.Afib * \
        instrument.throughput * \
        magellan.throughput(instrument.waves)

    # The input wavelengths are in the lab reference frame, the (1+2.5) shifts the galaxy to a redshift z=2.5.
    with fits.open(input_spectrum) as hdul:
        input_wv   = hdul[2].data * (1+2.5) / 10.0
        input_spec = hdul[0].data


    # This is a noramlization constant taken from the Shapley (2003) paper that presented this spectrum.
    fnu_shapley = 3.4e-30

    # Requirement is r=25 magnitudes
    m = 25.0
    h = 6.6e-27
    c = 3.0e17   # nm/sec

    # The flux of an object with AB magnitude of 25 (see definition of AB magnitudes)
    fnu = 10**(-0.4*(m+48.6))

    # Renormalize the spectrum to the correct magnitude
    input_spec *= fnu / fnu_shapley

    # Convert from fnu (erg/cm2/s/Hz) to flambda (erg/cm2/s/AA)
    input_spec *= c / (input_wv)**2

    # Convert from energy units (ergs) to counted photons
    input_photons = input_spec / (h * c/input_wv)

    objspec = np.interp(instrument.waves,input_wv,input_photons)

    # Same calculation as above, except assume all of the light goes down one fiber
    # So, you don't need to integrate over the fiber area, but do need everything else
    objphotons = objspec * \
        magellan.Atel * \
        texp * \
        (instrument.waves)/instrument.R * \
        instrument.throughput * \
        magellan.throughput(instrument.waves)

    readnoise = instrument.sensor.rn
    dark      = instrument.sensor.dark * texp

    # This calculates the number of pixels that a fiber subtends. When we simulate the extracted the spectrum, we need
    # to quadrature sum read noise from ALL pixels in the profile, not just one pixel.
    pix_resel = (instrument.fiber.dFib/1e3) * (instrument.f_cam / instrument.f_col) / instrument.sensor.pixelsize

    print ("Profile (pix) = ", pix_resel)
    print ("Dark noise/profile  = ", np.sqrt(dark * np.ceil(pix_resel)))
    print ("RN/profile    = ", readnoise * np.sqrt(np.ceil(pix_resel)))
    
    skynoise = np.sqrt(skyphotons)
    totnoise = np.sqrt(skyphotons + dark*np.ceil(pix_resel) + readnoise**2*np.ceil(pix_resel))

    # Not sure what that last return value is, must have been work in progress.
    return objphotons, totnoise, skynoise, np.sqrt(2*dark * np.ceil(pix_resel))

    plt.plot(instrument.waves, skynoise/totnoise)
    plt.ylim([0,1])
    plt.show()

def skynoise(texp=600.0):

    llamas_blue = spec.Spectrograph('LLAMAS_BLUE')
    llamas_green = spec.Spectrograph('LLAMAS_GREEN')
    llamas_red = spec.Spectrograph('LLAMAS_RED')

    llamas_blue.build_model('llamas_blue.def')
    llamas_green.build_model('llamas_green.def')
    llamas_red.build_model('llamas_red.def')

    objphotons_blue,  totnoise_blue,  skynoise_blue = observe(llamas_blue, texp)
    objphotons_green, totnoise_green, skynoise_green = observe(llamas_green, texp)
    objphotons_red,   totnoise_red,   skynoise_red = observe(llamas_red, texp)

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(llamas_blue.waves, objphotons_blue)
    ax.plot(llamas_blue.waves, totnoise_blue, alpha=0.5, color='b')
#    ax.plot(llamas_blue.waves, skynoise_blue**2, alpha=0.5, color='c')
    ax.plot(llamas_blue.waves, 0.01*skynoise_blue**2, alpha=0.5, color='c')

    ax.plot(llamas_green.waves, objphotons_green, color='g')
    ax.plot(llamas_green.waves, totnoise_green, alpha=0.5, color='g')
#    ax.plot(llamas_green.waves, skynoise_green**2, alpha=0.5, color='c')
    ax.plot(llamas_green.waves, 0.01*skynoise_green**2, alpha=0.5, color='c')

    ax.plot(llamas_red.waves, objphotons_red, color='r')
    ax.plot(llamas_red.waves, totnoise_red, alpha=0.2, color='r')
#    ax.plot(llamas_red.waves, skynoise_red**2, alpha=0.5, color='c')
    ax.plot(llamas_red.waves, 0.01*skynoise_red**2, alpha=0.2, color='c')

    ax.set_ylim([0,200*texp/3000.0])

    fig.savefig('lbg_sky_10hr.pdf')

    plt.show()

    stop


    ax.set_xlim(350,970)
    ax.set_ylim(0,10)

    ax.plot(llamas_blue.waves, skynoise_blue/llamas_blue.sensor.rn)
    ax.plot(llamas_green.waves, skynoise_green/llamas_green.sensor.rn)
    ax.plot(llamas_red.waves, skynoise_red/llamas_red.sensor.rn)
    ax.plot([350,970],[1.0,1.0],linestyle='-')
    ax.set_xlabel("Wavelength (nm)")
    ax.set_ylabel("Sky Poisson Noise / Read Noise")
    label = "texp (s) = ", texp
    ax.text(380,9,label)

    ax.plot([373.,373.],[0,2],color='k',alpha=0.5)
    ax.plot([486.,486.],[0,2],color='k',alpha=0.5)
    ax.plot([501.,501.],[0,2],color='k',alpha=0.5)
    ax.plot([656.,656.],[0,2],color='k',alpha=0.5)
    ax.plot([656.,656.],[0,2],color='k',alpha=0.5)
    ax.plot([850.,850.],[0,2],color='k',alpha=0.5)

#    fig.savefig("llamas_sky_readnoise.pdf")

#    plt.show()

    fig2 = plt.figure()
    ax = fig2.add_subplot(111)

    ax.set_xlim(350,970)
    ax.set_ylim(0,200*texp/600.0)

    ax.plot(llamas_blue.waves, skynoise_blue**2)
    ax.plot(llamas_green.waves, skynoise_green**2)
    ax.plot(llamas_red.waves, skynoise_red*2)
    ax.plot(llamas_blue.waves, totnoise_blue,alpha=0.3,color='b')
    ax.plot(llamas_green.waves, totnoise_green,alpha=0.3,color='g')
    ax.plot(llamas_red.waves, totnoise_red,alpha=0.3,color='r')
    
    plt.show()
