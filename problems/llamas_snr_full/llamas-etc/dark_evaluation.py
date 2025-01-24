import os
import sys
import spectrograph as spec
import matplotlib.pyplot as plt  
import observe
import numpy as np

if 'COATINGS_PATH' not in os.environ:
   print("Path variables not set, please check and run 'llamas_env.bat' and try again.")

llamas_red = spec.Spectrograph('LLAMAS_RED')
llamas_blue = spec.Spectrograph('LLAMAS_BLUE')
llamas_green = spec.Spectrograph('LLAMAS_GREEN') 

llamas_red.build_model('llamas_red.def')
llamas_blue.build_model('llamas_blue.def')
llamas_green.build_model('llamas_green.def')

objphotons_red,totnoise_red,skynoise_red,dark_red = observe.observe(llamas_red, 1800)
objphotons_green,totnoise_green,skynoise_green,dark_green = observe.observe(llamas_green, 1800)
objphotons_blue,totnoise_blue,skynoise_blue,dark_blue = observe.observe(llamas_blue, 1800)

#Requirements
llamas_blue.sensor.dark = 0.001
llamas_green.sensor.dark = 0.001
llamas_red.sensor.dark = 0.001
objphotons_red,totnoise_red_reqt,skynoise_red,dark_red_reqt = observe.observe(llamas_red, 1800)
objphotons_blue,totnoise_blue_reqt,skynoise_blue,dark_blue_reqt = observe.observe(llamas_blue, 1800)
objphotons_green,totnoise_green_reqt,skynoise_green,dark_green_reqt = observe.observe(llamas_green, 1800)

rtmp = np.sqrt(skynoise_red**2+dark_red**2)/skynoise_red
gtmp = np.sqrt(skynoise_green**2+dark_green**2)/skynoise_green
btmp = np.sqrt(skynoise_blue**2+dark_blue**2)/skynoise_blue

'''
plt.plot(rtmp,c='r')
plt.show()

plt.plot(gtmp,c='g')
plt.show()

plt.plot(btmp,c='b')
plt.show()
'''

#In-house measurements
llamas_blue.sensor.dark = 0.00267
llamas_green.sensor.dark = 0.00267
llamas_red.sensor.dark = 0.00238
llamas_blue.sensor.rn = 2.35
llamas_green.sensor.rn = 2.35
llamas_red.sensor.rn = 2.35
objphotons_red,totnoise_red_asbuilt,skynoise_red,dark_red_asbuilt = observe.observe(llamas_red, 1800)
objphotons_blue,totnoise_blue_asbuilt,skynoise_blue,dark_blue_asbuilt = observe.observe(llamas_blue, 1800)
objphotons_green,totnoise_green_asbuilt,skynoise_green,dark_green_asbuilt = observe.observe(llamas_green, 1800)

plt.plot(llamas_red.waves, totnoise_red_asbuilt/totnoise_red_reqt, color='r')
plt.plot(llamas_blue.waves, totnoise_blue_asbuilt/totnoise_blue_reqt, color='b')
plt.plot(llamas_green.waves, totnoise_green_asbuilt/totnoise_green_reqt, color='g')

plt.suptitle('As-measured at -90C: 0.002 (measured) vs 0.001 (reqt) e-/sec/pix')
plt.ylim(0.975,1.1)
plt.xlim(350,1000)
plt.xlabel("Wavelength (nm)")
plt.ylabel("Total Noise Ratio [(Dark=0.00267)/(Dark=0.001)]")

plt.show()

