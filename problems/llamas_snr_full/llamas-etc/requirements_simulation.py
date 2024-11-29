import spectrograph as spec
import numpy as np
import matplotlib.pyplot as plt
import observe
from astropy.io import fits
import scipy.signal as ss

import argparse
import os
import sys

def snr(signal, noise):
    return signal/noise


#################################
# Continuum Sensitivity
# R = 25 Lyman Break Galaxy SNR=3 in 5 hours
#
def continuum_req(do_plot, llamas_red, llamas_green, llamas_blue, llamas_waves, texp=5*3600):
    with fits.open('lbg_shapley.fits') as hdul:
        input_wv   = hdul[2].data * (1+2.5) / 10.0
        input_spec = hdul[0].data
    
    fnu_shapley = 3.4e-30
    Rmag = 25.0
    c = 3.0e17   # nm/sec
    fnu = 10**(-0.4*(Rmag+48.6))
    input_spec *= fnu / fnu_shapley
    input_flam = input_spec * c / (input_wv)**2 

    if do_plot: #plot spectra
        plt.plot(input_wv,input_flam, color = 'k')
        plt.title("Continuum Sensitivity (Shapley spectrum)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("flux [erg/cm2/s/Å]")
        plt.show()

    red_1_photons, red_1_noise = observe.observe_spectrum(llamas_red, texp, input_wv, input_flam)
    gre_1_photons, gre_1_noise = observe.observe_spectrum(llamas_green, texp, input_wv, input_flam)
    blu_1_photons, blu_1_noise = observe.observe_spectrum(llamas_blue, texp, input_wv, input_flam)
    spectral_npix = 4
    signal_1 = np.array(spectral_npix * np.concatenate([blu_1_photons,gre_1_photons,red_1_photons]))
    noise_1 = np.array(np.sqrt(spectral_npix) * np.concatenate([blu_1_noise,gre_1_noise,red_1_noise]))
    snr_1 = [signal_1[i]/noise_1[i] for i in range(0, len(signal_1))]
    
    if do_plot: #plot signal
        plt.plot(llamas_waves,signal_1)
        plt.title("Continuum Sensitivity (Shapley spectrum)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("counts")
        plt.show()

    print("CONTINUUM SENSITIVITY")
    print("Target SNR: 3")
    print("Median SNR: " + str(np.median(snr_1)))

    #plt.plot(llamas_waves, snr_1)
    if do_plot:
        plt.plot(llamas_red.waves, snr(red_1_photons,red_1_noise), color='r')
        plt.plot(llamas_green.waves, snr(gre_1_photons,gre_1_noise), color='g')
        plt.plot(llamas_blue.waves, snr(blu_1_photons,blu_1_noise), color='b')
        plt.plot(llamas_waves, llamas_waves*0 + 3, color='orange')
        plt.title("Continuum Sensitivity (SNR vs wavelength)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("SNR")
        plt.show()
        
    return np.mean(snr_1)


###########################################
# Transient Continuum Sensitivity
# R = 21 transient in 600s SNR10 per resel
#
# Template is in erg/cm2/sec/A x 1e15
def transient_continuum_req(do_plot, llamas_red, llamas_green, llamas_blue, llamas_waves, texp=600):
    sn1a_spec = np.genfromtxt('MBSE_SPECTRA/sn2000fa/sn2000fa-20001216.flm',dtype=[('wave','f8'),('flux','f8')])
    sn1a_spec['flux'] *= 1e-15
    Rmag_reference     = 16.4
    Rmag_requirement   = 21.0
    scalefac = 10**(0.4*(Rmag_reference-Rmag_requirement))
    sn1a_spec['flux'] *= scalefac
    sn1a_spec['wave'] /= 10.0 # convert to nm

    if do_plot: #plot spectra
        plt.plot(sn1a_spec['wave'],sn1a_spec['flux'], color = 'k')
        plt.title("Transient Continuum Sensitivity (Supernova spectrum)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("flux [erg/cm2/s/Å]")
        plt.show()

    red_2_photons, red_2_noise = observe.observe_spectrum(llamas_red, texp, sn1a_spec['wave'], sn1a_spec['flux'])
    gre_2_photons, gre_2_noise = observe.observe_spectrum(llamas_green, texp, sn1a_spec['wave'], sn1a_spec['flux'])
    blu_2_photons, blu_2_noise = observe.observe_spectrum(llamas_blue, texp, sn1a_spec['wave'], sn1a_spec['flux'])
    spectral_npix = 4
    signal_2 = np.array(spectral_npix * np.concatenate([blu_2_photons,gre_2_photons,red_2_photons]))
    noise_2 = np.array(np.sqrt(spectral_npix) * np.concatenate([blu_2_noise,gre_2_noise,red_2_noise]))
    snr_2 = [signal_2[i]/noise_2[i] for i in range(0, len(signal_2))]

    print("TRANSIENT CONTINUUM SENSITIVITY")
    print("Target SNR: 10")
    print("Median SNR: " + str(np.median(snr_2)))
    
    if do_plot: #plot signal
        plt.plot(llamas_waves,signal_2)
        plt.title("Transient Continuum Sensitivity (Supernova spectrum)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("counts")
        plt.show()

    #plt.plot(llamas_waves, snr_2)
    if do_plot:
        plt.plot(llamas_red.waves, snr(red_2_photons,red_2_noise), color='r')
        plt.plot(llamas_green.waves, snr(gre_2_photons,gre_2_noise), color='g')
        plt.plot(llamas_blue.waves, snr(blu_2_photons,blu_2_noise), color='b')
        plt.plot(llamas_waves, llamas_waves*0 + 10, color='orange')
        plt.title("Transient Continuum Sensitivity (SNR vs wavelength)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("SNR")
        plt.show()
        
    return np.mean(snr_2)


##########################################################
# Emission Sensitivity
# LSB emission 3e-18 erg/cm2/s/A/arcsec SNR=5 in 1 night or less
#
def emission_req(do_plot, llamas_red, llamas_green, llamas_blue, llamas_waves, texp=8*3600):
    sampling = 0.05 
    wave_oversampled = 300 + sampling * np.arange((1000-300)/sampling)
    spec_oversampled = wave_oversampled * 0 + 4.5e-18

    if do_plot: #plot spectra
        plt.plot(wave_oversampled, spec_oversampled, color = 'k')
        plt.title("Emission Sensitivity (Oversampled spectrum)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("flux [erg/cm2/s/Å]")
        plt.show()

    red_3_photons, red_3_noise = observe.observe_spectrum(llamas_red, texp, wave_oversampled, spec_oversampled)
    gre_3_photons, gre_3_noise = observe.observe_spectrum(llamas_green, texp, wave_oversampled, spec_oversampled)
    blu_3_photons, blu_3_noise = observe.observe_spectrum(llamas_blue, texp, wave_oversampled, spec_oversampled)
    spectral_npix = 4
    signal_3 = np.array(spectral_npix * np.concatenate([blu_3_photons,gre_3_photons,red_3_photons]))
    noise_3 = np.array(np.sqrt(spectral_npix) * np.concatenate([blu_3_noise,gre_3_noise,red_3_noise]))
    snr_3 = [signal_3[i]/noise_3[i] for i in range(0, len(signal_3))]

    print("EMISSION SENSITIVITY")
    print("Target SNR: 5")
    print("Median SNR: " + str(np.median(snr_3)))
    
    if do_plot: #plot signal
        plt.plot(llamas_waves,signal_3)
        plt.title("Emission Sensitivity (Oversampled spectrum)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("counts")
        plt.show()

    #plt.plot(llamas_waves, snr_3)
    if do_plot:
        plt.plot(llamas_red.waves, snr(red_3_photons,red_3_noise), color='r')
        plt.plot(llamas_green.waves, snr(gre_3_photons,gre_3_noise), color='g')
        plt.plot(llamas_blue.waves, snr(blu_3_photons,blu_3_noise), color='b')
        plt.plot(llamas_waves, llamas_waves*0 + 5, color='orange')
        plt.title("Emission Sensitivity (SNR vs wavelength)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("SNR")
        plt.show()

    return np.mean(snr_3)

    
################################################################
# Line Sensitivity
# SNR=5 in 3 hours for emission line with flux > 1e-18 erg/cm2/s for lambda < 700
#
def line_req(do_plot, llamas_red, llamas_green, llamas_blue, llamas_waves, texp=3*3600):
    sampling = 0.05
    wave_oversampled = 300 + sampling * np.arange((1000-300)/sampling)
    spec_oversampled = wave_oversampled * 0
    peakpos = [400,500,600,700,800,900]

    sigmas = []
    for pk in peakpos:
        d_lambda_doppler = pk * (200 / 3e8)
        d_lambda_instr = pk / 2200
        sigma = np.sqrt(d_lambda_doppler**2 + d_lambda_instr**2) #peak width
        sigmas.append(sigma)
        lineprof = 1e-18 * 1/(np.sqrt(2*np.pi)*sigma) * np.exp(-(wave_oversampled-pk)**2/(2*sigma**2))
        spec_oversampled += lineprof

    if do_plot: #plot spectra
        plt.plot(wave_oversampled, spec_oversampled, color = 'k')
        plt.title("Line Sensitivity (Comb spectrum)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("flux [erg/cm2/s/Å]")
        plt.show()

    red_4_photons, red_4_noise = observe.observe_spectrum(llamas_red, texp, wave_oversampled, spec_oversampled)
    gre_4_photons, gre_4_noise = observe.observe_spectrum(llamas_green, texp, wave_oversampled, spec_oversampled)
    blu_4_photons, blu_4_noise = observe.observe_spectrum(llamas_blue, texp, wave_oversampled, spec_oversampled)
    signal_4 = np.array(np.concatenate([blu_4_photons,gre_4_photons,red_4_photons]))
    noise_4 = np.array(np.concatenate([blu_4_noise,gre_4_noise,red_4_noise]))
    snr_4 = [signal_4[i]/noise_4[i] for i in range(0, len(signal_4))]

    print("LINE SENSITIVITY")
    print("Target Peak SNR: 5")
    peaks, props = ss.find_peaks(snr_4, distance=50)

    #how far out in nm we want to count the peak
    sig_nm = np.array([(sig*2) for sig in sigmas]) 
    #conversion from nm to px at each peak position
    nm_per_px = np.array([np.diff(llamas_waves)[x] for x in peaks]) 
    #how far out in px we want to count the peak
    num_px = [int(np.ceil(x)) for x in sig_nm / nm_per_px] 

    #adding up the snr (signal linearly, noise in quadrature) of each pixel in the peak
    peak_snrs = []
    for x,pk in enumerate(peaks):
        peak_signal, peak_noise = 0,0
        for n in range(0,num_px[x]): #count pixels to the left and middle
           peak_signal += signal_4[pk-n]
           peak_noise += (noise_4[pk-n])**2
        for n in range(1,num_px[x]): #count pixels to the right
           peak_signal += signal_4[pk+n]
           peak_noise += (noise_4[pk+n])**2
        peak_snrs.append(peak_signal / np.sqrt(peak_noise))
        
    print("Peak SNR's: ",peak_snrs)
    print("Median SNR: " + str(np.median(peak_snrs)))
    
    if do_plot: #plot signal
        plt.plot(llamas_waves,signal_4)
        plt.title("Line Sensitivity (Comb spectrum)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("counts")
        plt.show()

    #plt.plot(llamas_waves, snr_4)
    if do_plot:
        plt.plot(llamas_red.waves, snr(red_4_photons,red_4_noise), color='r')
        plt.plot(llamas_green.waves, snr(gre_4_photons,gre_4_noise), color='g')
        plt.plot(llamas_blue.waves, snr(blu_4_photons,blu_4_noise), color='b')
        plt.plot(llamas_waves, llamas_waves*0 + 5, color='orange')
        plt.scatter(peakpos, peak_snrs, color='k')
        plt.title("Line Sensitivity (SNR vs wavelength)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("SNR")
        plt.show()
        
    return np.mean([snr_4[x] for x in peaks])

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('req', type=str, help='Options are: ContinuumSensitivity, TransientContinuumSensitivity, EmissionSensitivity, LineSensitivity', nargs='?', default="")
    parser.add_argument('coatings_path', type=str, help='System global path of "llamas-etc/COATINGS"', nargs='?', default="C:/LLAMAS/SIMULATIONS/llamas-etc/COATINGS/")
    parser.add_argument('test_data_path', type=str, help='System global path of "LLAMAS/TEST DATA/"', nargs='?', default="C:/LLAMAS/TEST DATA/")
    parser.add_argument('index', type=str, help='Number of the spectrograph; options are 1..8; default is blank', nargs='?', default="")
    args = parser.parse_args()
    
    #set path variables
    os.environ["COATINGS_PATH"]=args.coatings_path
    os.environ["TEST_DATA_PATH"]=args.test_data_path
    
    #build model
    llamas_red = spec.Spectrograph('LLAMAS_RED')
    llamas_blue = spec.Spectrograph('LLAMAS_BLUE')
    llamas_green = spec.Spectrograph('LLAMAS_GREEN') 

    llamas_red.build_model('llamas_red'+args.index+'.def')
    llamas_blue.build_model('llamas_blue'+args.index+'.def')
    llamas_green.build_model('llamas_green'+args.index+'.def')
    llamas_waves = np.array(np.concatenate([llamas_blue.waves,llamas_green.waves,llamas_red.waves]))
    
    #call function
    if args.req == "":
        continuum_req(True, llamas_red, llamas_green, llamas_blue, llamas_waves)
        transient_continuum_req(True, llamas_red, llamas_green, llamas_blue, llamas_waves)
        emission_req(True, llamas_red, llamas_green, llamas_blue, llamas_waves)
        line_req(True, llamas_red, llamas_green, llamas_blue, llamas_waves)
    
    elif args.req == "ContinuumSensitivity":
        continuum_req(False, llamas_red, llamas_green, llamas_blue, llamas_waves)
    
    elif args.req == "TransientContinuumSensitivity":
        transient_continuum_req(False, llamas_red, llamas_green, llamas_blue, llamas_waves)
        
    elif args.req == "EmissionSensitivity":
        emission_req(False, llamas_red, llamas_green, llamas_blue, llamas_waves)
    
    elif args.req == "LineSensitivity":
        line_req(False, llamas_red, llamas_green, llamas_blue, llamas_waves)
        