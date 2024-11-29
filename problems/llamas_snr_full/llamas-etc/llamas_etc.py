import os
import sys
from llamas_plot_throughput import plot_throughput

llamas_path = "Q:\\Repositories\\LLAMAS\\"
sys.path.append(llamas_path + "SIMULATIONS\\llamas-etc-copy\\")
os.environ["COATINGS_PATH"]= llamas_path + "SIMULATIONS\\llamas-etc-copy\\COATINGS\\"
os.environ["TEST_DATA_PATH"]=llamas_path + "\\TEST DATA\\"
os.environ["REQUIREMENTS_PATH"]="Q:\\Dropbox (MIT)\\grad school\\LLAMAS\\LLAMAS\\"



def_list = [["llamas_red1.def","llamas_green1.def","llamas_blue1.def"]]


pmfile = "performance_measures.def"

##Choose whether or not to update performance measures file
if '--parse' in sys.argv:
   #Populate the temporary data files from MagicDraw requirements file
   exec(open("parse_requirements.py").read())


##Run calculations
import spectrograph as spec
import numpy as np
import matplotlib.pyplot as plt

'''
total_peak_med_thru
total_mid_min_thru
total_edge_min_thru
red_camera_avg_thru
green_camera_avg_thru
blue_camera_avg_thru
red_camera_min_thru
green_camera_min_thru
blue_camera_min_thru
fiberrun_med_thru
fiberrun_min_thru
spect_peak_med_thru
spect_mid_min_thru
spect_edge_min_thru
collimator_avg_thru
rg_dichroic_avg_thru
rg_dichroic_min_thru
bg_dichroic_avg_thru
bg_dichroic_min_thru
mla_min_thru
'''

print("Writing to report.txt...")
reportfile = open("report.txt", 'w+')
sys.stdout = reportfile

for num,spectro in enumerate(def_list):
   print("===SPECTROGRAPH " + str(num+1) + "===")

   llamas_red = spec.Spectrograph('LLAMAS_RED')
   llamas_red.build_model(spectro[0])
   llamas_green = spec.Spectrograph('LLAMAS_GREEN')
   llamas_green.build_model(spectro[1])
   llamas_blue = spec.Spectrograph('LLAMAS_BLUE')
   llamas_blue.build_model(spectro[2])

   waves = llamas_red.waves #they should all be the same!
   thru_all = np.add(llamas_blue.throughput, llamas_green.throughput, llamas_red.throughput) #all covering the same range - VPH's will enforce te spectral cutoffs
   thru_all_peak, thru_all_mid, thru_all_edge = [], [], []
   for wave,thru in zip(waves, thru_all):
      if wave >= 400.0 and wave <= 900.0:
         thru_all_peak.append(thru)
      if wave >= 350.0 and wave <= 970.0:
         thru_all_mid.append(thru)
      if wave < 400.0 or wave > 900.0:
         thru_all_edge.append(thru)
   total_peak_med_thru = np.median(thru_all_peak)
   total_mid_min_thru = np.amin(thru_all_mid)
   total_edge_min_thru = np.amin(thru_all_edge)
   
   thru_spect = np.append([],(llamas_blue.calc_throughput(llamas_blue.waves,nofront=True), \
                llamas_green.calc_throughput(llamas_green.waves,nofront=True), \
                llamas_red.calc_throughput(llamas_red.waves,nofront=True) ))
   thru_spect_peak, thru_spect_mid, thru_spect_edge = [], [], []
   for wave,thru in zip(waves, thru_spect):
      if wave >= 400.0 and wave <= 900.0:
         thru_spect_peak.append(thru)
      if wave >= 350.0 and wave <= 970.0:
         thru_spect_mid.append(thru)
      if wave < 400.0 or wave > 900.0:
         thru_spect_edge.append(thru)
   spect_peak_med_thru = np.median(thru_spect_peak)
   spect_mid_min_thru = np.amin(thru_spect_mid)
   spect_edge_min_thru = np.amin(thru_spect_edge)
   
   #"camera" = L1-L7 in this context
   blue_camera_thru = llamas_blue.elements[4].throughput(llamas_blue.waves) * \
                      llamas_blue.elements[5].throughput(llamas_blue.waves) * \
                      llamas_blue.elements[6].throughput(llamas_blue.waves) * \
                      llamas_blue.elements[7].throughput(llamas_blue.waves) * \
                      llamas_blue.elements[8].throughput(llamas_blue.waves) * \
                      llamas_blue.elements[9].throughput(llamas_blue.waves) * \
                      llamas_blue.elements[10].throughput(llamas_blue.waves)
   green_camera_thru = llamas_green.elements[4].throughput(llamas_green.waves) * \
                       llamas_green.elements[5].throughput(llamas_green.waves) *\
                       llamas_green.elements[6].throughput(llamas_green.waves) *\
                       llamas_green.elements[7].throughput(llamas_green.waves) *\
                       llamas_green.elements[8].throughput(llamas_green.waves) *\
                       llamas_green.elements[9].throughput(llamas_green.waves) *\
                       llamas_green.elements[10].throughput(llamas_green.waves)
   red_camera_thru = llamas_red.elements[3].throughput(llamas_red.waves) * \
                     llamas_red.elements[4].throughput(llamas_red.waves) *\
                     llamas_red.elements[5].throughput(llamas_red.waves) *\
                     llamas_red.elements[6].throughput(llamas_red.waves) *\
                     llamas_red.elements[7].throughput(llamas_red.waves) *\
                     llamas_red.elements[8].throughput(llamas_red.waves) *\
                     llamas_red.elements[9].throughput(llamas_red.waves)
   red_camera_avg_thru = np.average(red_camera_thru)
   green_camera_avg_thru = np.average(green_camera_thru)
   blue_camera_avg_thru = np.average(blue_camera_thru)
   red_camera_min_thru = np.amin(red_camera_thru)
   green_camera_min_thru = np.amin(green_camera_thru)
   blue_camera_min_thru = np.amin(blue_camera_thru)
   
   collimator_avg_thru = np.average(np.append([],( \
                         llamas_blue.elements[0].throughput(llamas_blue.waves),
                         llamas_green.elements[0].throughput(llamas_green.waves),
                         llamas_red.elements[0].throughput(llamas_red.waves) )))
   rg_dichroic_avg_thru = np.average(np.append([],(
                          llamas_blue.elements[1].throughput(llamas_blue.waves), \
                          llamas_green.elements[1].throughput(llamas_green.waves), \
                          llamas_red.elements[1].throughput(llamas_red.waves) )))
   rg_dichroic_min_thru = np.amin(np.append([],( \
                          llamas_blue.elements[1].throughput(llamas_blue.waves), \
                          llamas_green.elements[1].throughput(llamas_green.waves), \
                          llamas_red.elements[1].throughput(llamas_red.waves) )))
   bg_dichroic_avg_thru = np.average(np.append([],( \
                          llamas_blue.elements[2].throughput(llamas_blue.waves), \
                          llamas_green.elements[2].throughput(llamas_green.waves) )))
   bg_dichroic_min_thru = np.amin(np.append([],( \
                          llamas_blue.elements[2].throughput(llamas_blue.waves), \
                          llamas_green.elements[2].throughput(llamas_green.waves) )))
   
   fiberrun_med_thru = np.median(np.append([],( \
                       llamas_blue.fiber.throughput(llamas_blue.waves), \
                       llamas_green.fiber.throughput(llamas_green.waves), \
                       llamas_red.fiber.throughput(llamas_red.waves) )))
   fiberrun_min_thru = np.amin(np.append([],( \
                       llamas_blue.fiber.throughput(llamas_blue.waves), \
                       llamas_green.fiber.throughput(llamas_green.waves), \
                       llamas_red.fiber.throughput(llamas_red.waves) )))
   
   mla_min_thru = np.amin(np.append([],(
                  llamas_blue.fiber.elements[0].throughput(llamas_blue.waves), \
                  llamas_green.fiber.elements[0].throughput(llamas_green.waves), \
                  llamas_red.fiber.elements[0].throughput(llamas_red.waves) )))
   
   
   
   ##Check performance against requirements
   print("Checking against requirements...")
   targets = {}
   with open(pmfile, 'r') as f:
      lines = f.readlines()
      for line in lines:
         pair = line.split("\t")
         targets[pair[0]] = pair[1].rstrip()
   #print(targets)
   
   def check_req(name, value, operator, target_key):
   
      if operator=='<':
         if value < float(targets[target_key]):
            result = "PASS"
         else:
            result = "FAIL"
      elif operator=='<=':
         if value <= float(targets[target_key]):
            result = "PASS"
         else:
            result = "FAIL"
      elif operator=='>':
         if value > float(targets[target_key]):
            result = "PASS"
         else:
            result = "FAIL"
      elif operator=='>=':
         if value >= float(targets[target_key]):
            result = "PASS"
         else:
            result = "FAIL"
      elif operator=='==':
         if value == float(targets[target_key]):
            result = "PASS"
         else:
            result = "FAIL"
      elif operator=='!=':
         if value == float(targets[target_key]):
            result = "PASS"
         else:
            result = "FAIL"
      else:
         print("ERROR")
      print('{:<44s}{:<20s}{:<18s}{:<6s}{:<6s}'.format(name+":",str(value),"Required Value: ",targets[target_key], result))
   
   
   print("\nRequirement: LLAMAS Total Thoughput")
   check_req("Median total throughput 350-970nm", total_peak_med_thru*100, ">=", "TOTAL_PEAK_MED_THROUGH")
   check_req("Min. total throughput 400-900nm", total_mid_min_thru*100, ">=", "TOTAL_MID_MIN_THROUGH")
   check_req("Min. total throughput <350 >970nm", total_edge_min_thru*100, ">=", "TOTAL_EDGE_MIN_THROUGH")
   print("\n> Requirement: Spectrograph Thoughput")
   check_req("Median spectrograph throughput 350-970nm", spect_peak_med_thru*100, ">=", "SPECT_PEAK_MED_THROUGH")
   check_req("Min. spectrograph throughput 400-900nm", spect_mid_min_thru*100, ">=", "SPECT_MID_MIN_THROUGH")
   check_req("Min. spectrograph throughput <350 >970nm", spect_edge_min_thru*100, ">=", "SPECT_EDGE_MIN_THROUGH")
   print("\n>> Requirement: Collimator Thoughput")
   check_req("Avg. Collimator Throughput", collimator_avg_thru*100, ">=", "COLLIMATOR_AVG_THROUGH")
   print("\n>> Requirement: Dichroic Thoughput")
   check_req("Avg. RG Dichroic Throughput", rg_dichroic_avg_thru*100, ">=", "DICHROIC_AVG_THROUGH")
   check_req("Min. RG Dichroic Throughput", rg_dichroic_min_thru*100, ">=", "DICHROIC_MIN_THROUGH")
   check_req("Avg. BG Dichroic Throughput", bg_dichroic_avg_thru*100, ">=", "DICHROIC_AVG_THROUGH")
   check_req("Min. BG Dichroic Throughput", bg_dichroic_min_thru*100, ">=", "DICHROIC_MIN_THROUGH")
   print("\n>> Requirement: Camera Thoughput")
   check_req("Min. Red Camera Throughput", red_camera_min_thru*100, ">=", "CAMERA_MIN_THROUGH")
   check_req("Min. Green Camera Throughput", green_camera_min_thru*100, ">=", "CAMERA_MIN_THROUGH")
   check_req("Min. Blue Camera Throughput", blue_camera_min_thru*100, ">=", "CAMERA_MIN_THROUGH")
   check_req("Avg. Red Camera Throughput", red_camera_avg_thru*100, ">=", "CAMERA_AVG_THROUGH")
   check_req("Avg. Green Camera Throughput", green_camera_avg_thru*100, ">=", "CAMERA_AVG_THROUGH")
   check_req("Avg. Blue Camera Throughput", blue_camera_avg_thru*100, ">=", "CAMERA_AVG_THROUGH")
   print("\n> Requirement: Fiber-Run Thoughput")
   check_req("Median Fiber-Run Throughput", fiberrun_med_thru*100, ">=", "FIBER_RUN_MED_THROUGH")
   check_req("Min. Fiber-Run Throughput", fiberrun_min_thru*100, ">=", "FIBER_RUN_MIN_THROUGH")
   print("\n>> Requirement: MLA Thoughput")
   check_req("Min. MLA Throughput", mla_min_thru*100, ">=", "MLA_MIN_THROUGH")
   print("\n\n")

   if '--graph'  in sys.argv:
      plot_throughput(spectro[0], spectro[1], spectro[2])


reportfile.close()