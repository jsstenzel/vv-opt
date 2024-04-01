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
    return [signal[i]/noise[i] for i in range(0, len(signal))]


#################################
# Continuum Sensitivity
# R = 25 Lyman Break Galaxy SNR=3 in 5 hours
#
def continuum_req(do_plot, llamas_red, llamas_green, llamas_blue, llamas_waves):
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

    red_1_photons, red_1_noise = observe.observe_spectrum(llamas_red, 5*3600, input_wv, input_flam)
    gre_1_photons, gre_1_noise = observe.observe_spectrum(llamas_green, 5*3600, input_wv, input_flam)
    blu_1_photons, blu_1_noise = observe.observe_spectrum(llamas_blue, 5*3600, input_wv, input_flam)
    signal_1 = np.array(np.concatenate([blu_1_photons,gre_1_photons,red_1_photons]))
    noise_1 = np.array(np.concatenate([blu_1_noise,gre_1_noise,red_1_noise]))
    snr_1 = [signal_1[i]/noise_1[i] for i in range(0, len(signal_1))]
    
    if do_plot: #plot signal
        plt.plot(llamas_waves,signal_1)
        plt.title("Continuum Sensitivity (Shapley spectrum)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("counts")
        plt.show()

    print("CONTINUUM SENSITIVITY")
    print("Target SNR: 3")
    print("Average SNR: " + str(np.mean(snr_1)))

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
def transient_continuum_req(do_plot, llamas_red, llamas_green, llamas_blue, llamas_waves):
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

    red_2_photons, red_2_noise = observe.observe_spectrum(llamas_red, 600, sn1a_spec['wave'], sn1a_spec['flux'])
    gre_2_photons, gre_2_noise = observe.observe_spectrum(llamas_green, 600, sn1a_spec['wave'], sn1a_spec['flux'])
    blu_2_photons, blu_2_noise = observe.observe_spectrum(llamas_blue, 600, sn1a_spec['wave'], sn1a_spec['flux'])
    signal_2 = np.array(np.concatenate([blu_2_photons,gre_2_photons,red_2_photons]))
    noise_2 = np.array(np.concatenate([blu_2_noise,gre_2_noise,red_2_noise]))
    snr_2 = [signal_2[i]/noise_2[i] for i in range(0, len(signal_2))]

    print("TRANSIENT CONTINUUM SENSITIVITY")
    print("Target SNR: 10")
    print("Average SNR: " + str(np.mean(snr_2)))
    
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
# LSB emission 3e-19 erg/cm2/s/A/fiber SNR=5 in 1 night or less
#
def emission_req(do_plot, llamas_red, llamas_green, llamas_blue, llamas_waves):
    sampling = 0.05 
    wave_oversampled = 300 + sampling * np.arange((1000-300)/sampling)
    spec_oversampled = wave_oversampled * 0 + 3e-19

    if do_plot: #plot spectra
        plt.plot(wave_oversampled, spec_oversampled, color = 'k')
        plt.title("Emission Sensitivity (Oversampled spectrum)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("flux [erg/cm2/s/Å]")
        plt.show()

    red_3_photons, red_3_noise = observe.observe_spectrum(llamas_red, 8*3600, wave_oversampled, spec_oversampled)
    gre_3_photons, gre_3_noise = observe.observe_spectrum(llamas_green, 8*3600, wave_oversampled, spec_oversampled)
    blu_3_photons, blu_3_noise = observe.observe_spectrum(llamas_blue, 8*3600, wave_oversampled, spec_oversampled)
    signal_3 = np.array(np.concatenate([blu_3_photons,gre_3_photons,red_3_photons]))
    noise_3 = np.array(np.concatenate([blu_3_noise,gre_3_noise,red_3_noise]))
    snr_3 = [signal_3[i]/noise_3[i] for i in range(0, len(signal_3))]

    print("EMISSION SENSITIVITY")
    print("Target SNR: 5")
    print("Average SNR: " + str(np.mean(snr_3)))
    
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
# SNR=5 in 3 hours for emission line with flux > 1e-18 erg/cm2/s
#
def line_req(do_plot, llamas_red, llamas_green, llamas_blue, llamas_waves):
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

    red_4_photons, red_4_noise = observe.observe_spectrum(llamas_red, 3*3600, wave_oversampled, spec_oversampled)
    gre_4_photons, gre_4_noise = observe.observe_spectrum(llamas_green, 3*3600, wave_oversampled, spec_oversampled)
    blu_4_photons, blu_4_noise = observe.observe_spectrum(llamas_blue, 3*3600, wave_oversampled, spec_oversampled)
    signal_4 = np.array(np.concatenate([blu_4_photons,gre_4_photons,red_4_photons]))
    noise_4 = np.array(np.concatenate([blu_4_noise,gre_4_noise,red_4_noise]))
    snr_4 = [signal_4[i]/noise_4[i] for i in range(0, len(signal_4))]

    print("LINE SENSITIVITY")
    print("Target Peak SNR: 5")
    peaks, props = ss.find_peaks(snr_4)
    print("Average SNR: " + str(np.mean([snr_4[x] for x in peaks])))
    
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
        plt.title("Line Sensitivity (SNR vs wavelength)")
        plt.xlabel("wavelength [nm]")
        plt.ylabel("SNR")
        plt.show()
        
    return np.mean(snr_4)

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('req', type=str, help='Options are: ContinuumSensitivity, TransientContinuumSensitivity, EmissionSensitivity, LineSensitivity', nargs='?', default="")
    parser.add_argument('coatings_path', type=str, help='System global path of "llamas-etc/COATINGS"', nargs='?', default="Q:/Repositories/LLAMAS/SIMULATIONS/llamas-etc/COATINGS/")
    parser.add_argument('test_data_path', type=str, help='System global path of "LLAMAS/TEST DATA/"', nargs='?', default="Q:/Repositories/LLAMAS/TEST DATA/")
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
        