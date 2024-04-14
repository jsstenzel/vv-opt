import spectrograph as spec
import numpy as np
import matplotlib.pyplot as plt
import observe
from astropy.io import fits
import scipy.signal as ss

import argparse
import os
import sys

os.environ['COATINGS_PATH']="C:/Users/June/Desktop/Repositories/LLAMAS/SIMULATIONS/llamas-etc-uq/COATINGS/"
os.environ['TEST_DATA_PATH']="C:/Users/June/Desktop/Repositories/LLAMAS/TEST DATA/"

def get_snr(signal,noise):
    return np.array([signal[i]/noise[i] for i in range(0, len(signal))])
    

#################################
# Continuum Sensitivity
# R = 25 Lyman Break Galaxy SNR=3 in 5 hours
#
def continuum_req(do_plot, llamas_red, llamas_green, llamas_blue, llamas_waves, texp = 5*3600):
    with fits.open('lbg_shapley.fits') as hdul:
        input_wv   = hdul[2].data * (1+2.5) / 10.0
        input_spec = hdul[0].data
    
    fnu_shapley = 3.4e-30
    Rmag = 25.0
    c = 3.0e17   # nm/sec
    fnu = 10**(-0.4*(Rmag+48.6))
    input_spec *= fnu / fnu_shapley
    input_flam = input_spec * c / (input_wv)**2 
    print("continuum min: "+str(np.min(input_flam[3:])),flush=True)
    print("continuum avg: "+str(np.mean(input_flam)))

    if do_plot: #plot spectra
        plt.plot(input_wv,input_flam, color = 'k')
        plt.title("Continuum Sensitivity (Shapley spectrum)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("flux [erg/cm2/s/Å]")
        plt.show()

    red_photons, red_noise = observe.observe_spectrum(llamas_red, texp, input_wv, input_flam)
    gre_photons, gre_noise = observe.observe_spectrum(llamas_green, texp, input_wv, input_flam)
    blu_photons, blu_noise = observe.observe_spectrum(llamas_blue, texp, input_wv, input_flam)
    signal = red_photons + gre_photons + blu_photons
    red_snr = get_snr(red_photons, red_noise)
    gre_snr = get_snr(gre_photons, gre_noise)
    blu_snr = get_snr(blu_photons, blu_noise)
    snr = red_snr + gre_snr + blu_snr
    
    if do_plot: #plot signal
        plt.plot(llamas_waves,signal)
        plt.title("Continuum Sensitivity (Shapley spectrum)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("counts")
        plt.show()

    print("CONTINUUM SENSITIVITY")
    print("Target SNR: 3")
    print("Average SNR: " + str(np.mean(snr)))

    #plt.plot(llamas_waves, snr)
    if do_plot:
        plt.plot(llamas_waves, get_snr(red_photons,red_noise), color='r')
        plt.plot(llamas_waves, get_snr(gre_photons,gre_noise), color='g')
        plt.plot(llamas_waves, get_snr(blu_photons,blu_noise), color='b')
        plt.plot(llamas_waves, snr, color='k')
        plt.plot(llamas_waves, llamas_waves*0 + 3, color='orange')
        plt.title("Continuum Sensitivity (SNR vs wavelength)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("SNR")
        plt.show()
        
    return np.mean(snr)


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
    print("continuum min: "+str(np.min(sn1a_spec['flux'])))
    print("continuum avg: "+str(np.mean(sn1a_spec['flux'])))

    if do_plot: #plot spectra
        plt.plot(sn1a_spec['wave'],sn1a_spec['flux'], color = 'k')
        plt.title("Transient Continuum Sensitivity (Supernova spectrum)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("flux [erg/cm2/s/Å]")
        plt.show()

    red_photons, red_noise = observe.observe_spectrum(llamas_red, texp, sn1a_spec['wave'], sn1a_spec['flux'])
    gre_photons, gre_noise = observe.observe_spectrum(llamas_green, texp, sn1a_spec['wave'], sn1a_spec['flux'])
    blu_photons, blu_noise = observe.observe_spectrum(llamas_blue, texp, sn1a_spec['wave'], sn1a_spec['flux'])
    signal = blu_photons + gre_photons + red_photons
    noise = blu_noise + gre_noise + red_noise
    snr = [signal[i]/noise[i] for i in range(0, len(signal))]

    print("TRANSIENT CONTINUUM SENSITIVITY")
    print("Target SNR: 10")
    print("Average SNR: " + str(np.mean(snr)))
    
    if do_plot: #plot signal
        plt.plot(llamas_waves,signal)
        plt.title("Transient Continuum Sensitivity (Supernova spectrum)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("counts")
        plt.show()

    #plt.plot(llamas_waves, snr)
    if do_plot:
        plt.plot(llamas_red.waves, get_snr(red_photons,red_noise), color='r')
        plt.plot(llamas_green.waves, get_snr(gre_photons,gre_noise), color='g')
        plt.plot(llamas_blue.waves, get_snr(blu_photons,blu_noise), color='b')
        plt.plot(llamas_waves, snr, color='k')
        plt.plot(llamas_waves, llamas_waves*0 + 10, color='orange')
        plt.title("Transient Continuum Sensitivity (SNR vs wavelength)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("SNR")
        plt.show()
        
    return np.mean(snr)


##########################################################
# Emission Sensitivity
# LSB emission 3e-19 erg/cm2/s/A/fiber SNR=5 in 1 night or less
#
def emission_req(do_plot, llamas_red, llamas_green, llamas_blue, llamas_waves, texp=8*3600):
    sampling = 0.05 
    wave_oversampled = 300 + sampling * np.arange((1000-300)/sampling)
    spec_oversampled = wave_oversampled * 0 + 3e-19

    if do_plot: #plot spectra
        plt.plot(wave_oversampled, spec_oversampled, color = 'k')
        plt.title("Emission Sensitivity (Oversampled spectrum)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("flux [erg/cm2/s/Å]")
        plt.show()

    red_photons, red_noise = observe.observe_spectrum(llamas_red, texp, wave_oversampled, spec_oversampled)
    gre_photons, gre_noise = observe.observe_spectrum(llamas_green, texp, wave_oversampled, spec_oversampled)
    blu_photons, blu_noise = observe.observe_spectrum(llamas_blue, texp, wave_oversampled, spec_oversampled)
    signal = blu_photons + gre_photons + red_photons
    noise = blu_noise + gre_noise + red_noise
    snr = [signal[i]/noise[i] for i in range(0, len(signal))]

    print("EMISSION SENSITIVITY")
    print("Target SNR: 5")
    print("Average SNR: " + str(np.mean(snr)))
    
    if do_plot: #plot signal
        plt.plot(llamas_waves,signal)
        plt.title("Emission Sensitivity (Oversampled spectrum)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("counts")
        plt.show()

    #plt.plot(llamas_waves, snr)
    if do_plot:
        plt.plot(llamas_red.waves, get_snr(red_photons,red_noise), color='r')
        plt.plot(llamas_green.waves, get_snr(gre_photons,gre_noise), color='g')
        plt.plot(llamas_blue.waves, get_snr(blu_photons,blu_noise), color='b')
        plt.plot(llamas_waves, snr, color='k')
        plt.plot(llamas_waves, llamas_waves*0 + 5, color='orange')
        plt.title("Emission Sensitivity (SNR vs wavelength)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("SNR")
        plt.show()

    return np.mean(snr)

    
################################################################
# Line Sensitivity
# SNR=5 in 3 hours for emission line with flux > 1e-18 erg/cm2/s
#
def line_req(do_plot, llamas_red, llamas_green, llamas_blue, llamas_waves, texp=3*3600):
    sampling = 0.05
    wave_oversampled = 300 + sampling * np.arange((1000-300)/sampling)
    spec_oversampled = wave_oversampled * 0
    peakpos = [400,500,600,700,800,900]

    for pk in peakpos:
        sigma = 0.5
        lineprof = 1e-18 * 1/(np.sqrt(2*np.pi)*sigma) * np.exp(-(wave_oversampled-pk)**2/(2*sigma**2))
        spec_oversampled += lineprof

    if do_plot: #plot spectra
        plt.plot(wave_oversampled, spec_oversampled, color = 'k')
        plt.title("Line Sensitivity (Comb spectrum)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("flux [erg/cm2/s/Å]")
        plt.show()

    red_photons, red_noise = observe.observe_spectrum(llamas_red, texp, wave_oversampled, spec_oversampled)
    gre_photons, gre_noise = observe.observe_spectrum(llamas_green, texp, wave_oversampled, spec_oversampled)
    blu_photons, blu_noise = observe.observe_spectrum(llamas_blue, texp, wave_oversampled, spec_oversampled)
    signal = blu_photons + gre_photons + red_photons
    noise = blu_noise + gre_noise + red_noise
    snr = [signal[i]/noise[i] for i in range(0, len(signal))]

    print("LINE SENSITIVITY")
    print("Target Peak SNR: 5")
    peaks, props = ss.find_peaks(snr)
    print("Average SNR: " + str(np.mean([snr[x] for x in peaks])))
    
    if do_plot: #plot signal
        plt.plot(llamas_waves,signal)
        plt.title("Line Sensitivity (Comb spectrum)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("counts")
        plt.show()

    #plt.plot(llamas_waves, snr)
    if do_plot:
        plt.plot(llamas_red.waves, get_snr(red_photons,red_noise), color='r')
        plt.plot(llamas_green.waves, get_snr(gre_photons,gre_noise), color='g')
        plt.plot(llamas_blue.waves, get_snr(blu_photons,blu_noise), color='b')
        plt.plot(llamas_waves, snr, color='k')
        plt.plot(llamas_waves, llamas_waves*0 + 5, color='orange')
        plt.title("Line Sensitivity (SNR vs wavelength)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("SNR")
        plt.show()
        
    return np.mean(snr)
    
    
def run_snr(do_plot, red_def, green_def, blue_def, texp=8*3600):
    sampling = 0.05 
    wave_oversampled = 300 + sampling * np.arange((1000-300)/sampling)
    spec_oversampled = wave_oversampled * 0 + 3e-19

    llamas_red = spec.Spectrograph('LLAMAS_RED')
    llamas_blue = spec.Spectrograph('LLAMAS_BLUE')
    llamas_green = spec.Spectrograph('LLAMAS_GREEN')
    llamas_red.build_model(red_def)
    llamas_blue.build_model(blue_def)
    llamas_green.build_model(green_def)
    llamas_waves = llamas_green.waves
    
    red_photons, red_noise = observe.observe_spectrum(llamas_red, texp, wave_oversampled, spec_oversampled)
    gre_photons, gre_noise = observe.observe_spectrum(llamas_green, texp, wave_oversampled, spec_oversampled)
    blu_photons, blu_noise = observe.observe_spectrum(llamas_blue, texp, wave_oversampled, spec_oversampled)
    signal = blu_photons + gre_photons + red_photons
    noise = blu_noise + gre_noise + red_noise
    snr = [signal[i]/noise[i] for i in range(0, len(signal))]

    #print("EMISSION SENSITIVITY")
    #print("Target SNR: 5")
    #print("Average SNR: " + str(np.mean(snr)))

    #plt.plot(llamas_waves, snr)
    """
    if do_plot:
        plt.plot(llamas_red.waves, get_snr(red_photons,red_noise), color='r')
        plt.plot(llamas_green.waves, get_snr(gre_photons,gre_noise), color='g')
        plt.plot(llamas_blue.waves, get_snr(blu_photons,blu_noise), color='b')
        plt.plot(llamas_waves, snr, color='k')
        plt.plot(llamas_waves, llamas_waves*0 + 5, color='orange')
        plt.title("Emission Sensitivity (SNR vs wavelength)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("SNR")
        plt.show()
    """

    return np.mean(snr) #this needs to be done differently

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('req', type=str, help='Options are: ContinuumSensitivity, TransientContinuumSensitivity, EmissionSensitivity, LineSensitivity', nargs='?', default="")
    parser.add_argument('coatings_path', type=str, help='System global path of "llamas-etc/COATINGS"', nargs='?', default="C:/LLAMAS/SIMULATIONS/llamas-etc/COATINGS/")
    parser.add_argument('test_data_path', type=str, help='System global path of "LLAMAS/TEST DATA/"', nargs='?', default="C:/LLAMAS/TEST DATA/")
    parser.add_argument('index', type=str, help='Number of the spectrograph; options are 1..8; default is blank', nargs='?', default="1")
    args = parser.parse_args()
    
    #set path variables
    #os.environ["COATINGS_PATH"]=args.coatings_path
    #os.environ["TEST_DATA_PATH"]=args.test_data_path
    #os.environ['COATINGS_PATH']="C:/Users/June/Desktop/Repositories/LLAMAS/SIMULATIONS/llamas-etc/COATINGS/"
    #os.environ['TEST_DATA_PATH']="C:/Users/June/Desktop/Repositories/LLAMAS/TEST DATA/"

    #build model
    llamas_red = spec.Spectrograph('LLAMAS_RED')
    llamas_blue = spec.Spectrograph('LLAMAS_BLUE')
    llamas_green = spec.Spectrograph('LLAMAS_GREEN') 

    llamas_red.build_model('llamas_red'+args.index+'.def')
    llamas_blue.build_model('llamas_blue'+args.index+'.def')
    llamas_green.build_model('llamas_green'+args.index+'.def')
    llamas_waves = llamas_green.waves
    
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
        