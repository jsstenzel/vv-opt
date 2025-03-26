import numpy as np
import os
import math
import matplotlib.pyplot as plt
import surfaces
#import importlib

class Spectrograph:
    def __init__(self,name):
        self.name       =    name
        self.R          =    0.0
        self.wv_min     =    0.0
        self.wv_max     =    0.0
        self.f_col      =    0.0
        self.f_cam      =    0.0
        self.Fratio_col =    0.0
        self.Fratio_cam =    0.0
        self.Dbeam      =    0.0
        self.blade_obscuration = 1.0
        self.nelements  =    0
        self.elements   =    []           # Array of optical element objects
        self.sensor     =    []
        self.grating    =    []
        self.fiber      =    FiberFeed()
        self.waves      =    0.0
        self.throughput =    0.0

    def build_model(self,config_file):
        # Read in the instrument model parameters from file
        with open(config_file) as input:
            for line in input:
                arr = line.split()
                if (arr[0] == '#'):
                    continue
                elif (arr[0] == 'RESOLUTION'):
                    self.R = float(arr[1])
                elif (arr[0] == 'F_COL'):
                    self.f_col = float(arr[1])
                elif (arr[0] == 'F_CAM'):
                    self.f_cam = float(arr[1])
                elif (arr[0] == 'FRATIO_COL'):
                    self.Fratio_col = float(arr[1])
                elif (arr[0] == 'FRATIO_CAM'):
                    self.FRatio_cam = float(arr[1])
                elif (arr[0] == 'D_BEAM'):
                    self.D_beam = float(arr[1])
                elif (arr[0] == 'WAVE_MIN'):
                    self.wv_min = float(arr[1])
                elif (arr[0] == 'WAVE_MAX'):
                    self.wv_max = float(arr[1])
                elif (arr[0] == 'BLADE_OBSCURE'):
                    self.blade_obscuration = float(arr[1])
                elif (arr[0] == "SPECTROGRAPH"):
                    print("\n",arr[1])
                elif (arr[0] == 'SENSOR'):
                    self.sensor = Sensor(arr[1],arr[2],float(arr[3]),
                                         float(arr[4]),float(arr[5]))
                elif (arr[0] == 'GRATING'):
                    self.grating = surfaces.DiffractionGrating \
                        (arr[1],float(arr[2]),float(arr[3]),float(arr[4]),arr[5],arr[6],arr[7])
                elif (arr[0] == 'ELEMENT'):
                    self.nelements += 1
                    #importlib.reload(surfaces)
                    tmp = surfaces.OpticalElement(arr[1],int(arr[2]),arr[3],
                                                  arr[4],arr[5],arr[6])
                    self.elements.append(tmp)

                elif(arr[0] == 'FIBER_THETA'):
                    self.fiber.theta_fib = float(arr[1])
                    self.fiber.Afib = np.pi * (self.fiber.theta_fib/2.0)**2
                elif(arr[0] == 'FIBER_DCORE'):
                    self.fiber.dFib = float(arr[1])
                elif(arr[0] == 'FIBER_FRD'):
                    self.fiber.frd_loss = float(arr[1])
                elif(arr[0] == 'FIBER_ELEM'):
                    self.fiber.nelements += 1
                    #importlib.reload(surfaces)
                    tmp = surfaces.OpticalElement(arr[1],int(arr[2]),arr[3],
                                                  arr[4],arr[5],arr[6])
                    self.fiber.elements.append(tmp)


            # Now the parameters are read in, so calculate the central wavelength
            # of each pixel and store along with the total throguhput.
            #disp = (self.wv_max-self.wv_min)/float(self.sensor.naxis1)
            #self.waves  = self.wv_min+np.arange(self.sensor.naxis1)*disp
            self.waves = self.calculate_composite_waves()
            self.throughput = self.calc_throughput(self.waves)
            
    #A hack I'm implementing to fix the throughput interpolations
    #This is because I changed how spectrograph works to use the instrument bandpass instead of limb bandpass
    #now, the throughput wavelengths will be correctly interpolated into the 2048*3 pixels
    def calculate_composite_waves(self):
        blue_low = 350.0
        blue_high = 480.0
        green_low = 480.0
        green_high = 690.0
        red_low = 690.0
        red_high = 975.0
        
        blue_disp = (blue_high-blue_low)/float(self.sensor.naxis1)
        green_disp = (green_high-green_low)/float(self.sensor.naxis1)
        red_disp = (red_high-red_low)/float(self.sensor.naxis1)
        
        blue_waves = blue_low+np.arange(self.sensor.naxis1)*blue_disp
        green_waves = green_low+np.arange(self.sensor.naxis1)*green_disp
        red_waves = red_low+np.arange(self.sensor.naxis1)*red_disp
        composite_waves = np.concatenate((blue_waves,green_waves,red_waves))
        return composite_waves

    def calc_throughput(self,input_wave,nofront=False):
        composite_throughput = np.ones(len(input_wave))
        # Glass Coating transmission
        for i in range(self.nelements):
            composite_throughput *= self.elements[i].throughput(input_wave)
        # Factor in grating blaze function
        composite_throughput *= self.grating.blaze.throughput(input_wave)
        # Factor in CCD QE
        composite_throughput *= self.sensor.qe.throughput(input_wave)
        # Factor in the Fiber Feed
        if (nofront == False):
            composite_throughput *= self.fiber.throughput(input_wave)
        # Vignetting from the fiber blade
        composite_throughput *= self.blade_obscuration
        
        for i,t in enumerate(composite_throughput):
            if t <= 0.0:
                composite_throughput[i] = 0.0
            elif t >= 1.0:
                composite_throughput[i] = 1.0
        
        return composite_throughput

###############################

class FiberFeed:
    def __init__(self):
        self.theta_fib = 0     # Projected diameter on the sky 
        self.Afib      = np.pi * (self.theta_fib/2.0)**2 # Solid angle
        self.dFib      = 0     # Physical diameter
        self.frd_loss  = 0.00  # Throguhput loss from focal ratio deg.
        self.nelements = 0
        self.elements  = []

    def throughput(self, input_wave):
        composite_throughput = np.ones(len(input_wave))
        composite_throughput *= (1-self.frd_loss)
        for i in range(self.nelements):
            composite_throughput *= self.elements[i].throughput(input_wave)
        return composite_throughput

#############################

class Sensor:
    def __init__(self,name,qefile,rn,dark,pixelsize,vendor='e2v'):
        print("   ", name, qefile, rn, dark)
        self.rn           = rn    # e- RMS
        self.dark         = dark  # e-/sec/pix
        self.pixelsize    = pixelsize # mm
        self.bin_spatial  = 1.0   
        self.bin_spectral = 1.0
        self.gain         = 1.0
        self.naxis1       = 2048.0
        self.naxis2       = 2048.0
        self.qe = surfaces.OpticalSurface(name,'CCD',qefile,vendor,status='model')
        self.qe.loadThroughputTab()
