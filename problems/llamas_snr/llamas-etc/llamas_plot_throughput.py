import spectrograph as spec
import numpy as np
import matplotlib.pyplot as plt

def plot_throughput(rfile, gfile, bfile):
   llamas_blue = spec.Spectrograph('LLAMAS_BLUE')
   llamas_blue.build_model(bfile)
   plt.plot(llamas_blue.waves, llamas_blue.throughput, c='b')

   llamas_green = spec.Spectrograph('LLAMAS_GREEN')
   llamas_green.build_model(gfile)
   plt.plot(llamas_green.waves, llamas_green.throughput, c='g')

   llamas_red = spec.Spectrograph('LLAMAS_RED')
   llamas_red.build_model(rfile)
   plt.plot(llamas_red.waves, llamas_red.throughput, c='r')

   waves_all = [llamas_blue.waves,llamas_green.waves,llamas_red.waves]
   thru_all = [llamas_blue.throughput,llamas_green.throughput,llamas_red.throughput]
   print(np.median(thru_all))
   
   plt.hlines(y=0.38, xmin=350, xmax=970, linewidth=2, color='orange')
   plt.hlines(y=0.3, xmin=400, xmax=900, linewidth=2, color='black')
   plt.hlines(y=0.1, xmin=350, xmax=400, linewidth=2, color='black')
   plt.hlines(y=0.1, xmin=900, xmax=970, linewidth=2, color='black')

   plt.plot(llamas_blue.waves, llamas_blue.throughput + llamas_green.throughput + llamas_red.throughput, linestyle=(0, (4,4)), c='k')

   plt.title('Total Throughput')
   plt.ylim(-0.05,1.05)
   plt.xlabel("Wavelength [nm]")
   plt.ylabel("Throughput ratio")
   plt.show()

   print("FIBER RUN")
   plt.plot(llamas_blue.waves, llamas_blue.fiber.throughput(llamas_blue.waves), c='b')
   plt.plot(llamas_green.waves, llamas_green.fiber.throughput(llamas_green.waves), c='g')
   plt.plot(llamas_red.waves, llamas_red.fiber.throughput(llamas_red.waves), c='r')

   fiber_run_med = np.median([llamas_blue.fiber.throughput(llamas_blue.waves), \
                    llamas_green.fiber.throughput(llamas_green.waves), \
                    llamas_red.fiber.throughput(llamas_red.waves)])
   print(fiber_run_med)

   plt.title('Fiber Run Throughput')
   plt.ylim(-0.05,1.05)
   plt.xlabel("Wavelength [nm]")
   plt.ylabel("Throughput ratio")
   plt.show()
   
   plt.title('Fiber Run Throughput')
   plt.ylim(0,1)
   plt.xlabel("Wavelength [nm]")
   plt.ylabel("Throughput ratio")
   plt.show()

   print("SPECTROGRAPH")
   plt.plot(llamas_blue.waves, llamas_blue.calc_throughput(llamas_blue.waves,nofront=True), c='b')
   plt.plot(llamas_green.waves, llamas_green.calc_throughput(llamas_green.waves,nofront=True), c='g')
   plt.plot(llamas_red.waves, llamas_red.calc_throughput(llamas_red.waves,nofront=True), c='r')
   
   plt.plot(llamas_blue.waves, llamas_blue.calc_throughput(llamas_blue.waves,nofront=True) + llamas_green.calc_throughput(llamas_green.waves,nofront=True) + llamas_red.calc_throughput(llamas_red.waves,nofront=True), linestyle=(0, (4,4)), c='k')


   med = np.median([llamas_blue.calc_throughput(llamas_blue.waves,nofront=True), \
                    llamas_green.calc_throughput(llamas_green.waves,nofront=True), \
                    llamas_red.calc_throughput(llamas_red.waves,nofront=True)])
   print(med)
   red_med = np.median(llamas_red.calc_throughput(llamas_red.waves,nofront=True))
   green_med = np.median(llamas_red.calc_throughput(llamas_red.waves,nofront=True))
   blue_med = np.median(llamas_red.calc_throughput(llamas_red.waves,nofront=True))
   
   plt.hlines(y=0.45, xmin=350, xmax=970, linewidth=2, color='orange')
   plt.hlines(y=0.2, xmin=400, xmax=900, linewidth=2, color='black')
   plt.hlines(y=0.15, xmin=350, xmax=400, linewidth=2, color='black')
   plt.hlines(y=0.15, xmin=900, xmax=970, linewidth=2, color='black')

   plt.title('Spectrograph Throughput')
   plt.ylim(-0.05,1.05)
   plt.xlabel("Wavelength [nm]")
   plt.ylabel("Throughput ratio")
   plt.show()

   print("PURE FIBER")

   # This is the AR coating + internal transmission
   b = llamas_blue.fiber.elements[1].throughput(llamas_blue.waves) * \
       llamas_blue.fiber.elements[2].throughput(llamas_blue.waves)

   g = llamas_green.fiber.elements[1].throughput(llamas_green.waves) * \
       llamas_green.fiber.elements[2].throughput(llamas_green.waves)

   r = llamas_red.fiber.elements[1].throughput(llamas_red.waves) * \
       llamas_red.fiber.elements[2].throughput(llamas_red.waves)

   w = [llamas_blue.waves,llamas_green.waves,llamas_red.waves]
   ft = [b,g,r]
   plt.plot(w,ft)
   plt.title('Fiber Throughput')
   plt.ylim(-0.05,1.05)
   plt.xlabel("Wavelength [nm]")
   plt.ylabel("Throughput ratio")
   plt.show()

   print(np.median(ft))

   print("DICHROICS")
   rg_dichroic_b = llamas_blue.elements[1].throughput(llamas_blue.waves)
   rg_dichroic_g = llamas_green.elements[1].throughput(llamas_green.waves)
   rg_dichroic_r = llamas_red.elements[1].throughput(llamas_red.waves)

   bg_dichroic_b = llamas_blue.elements[2].throughput(llamas_blue.waves)
   bg_dichroic_g = llamas_green.elements[2].throughput(llamas_green.waves)

   plt.plot(llamas_blue.waves, rg_dichroic_b, c='b')
   plt.plot(llamas_green.waves, rg_dichroic_g, c='g')
   plt.plot(llamas_red.waves, rg_dichroic_r, c='r')

   plt.title('RG Dichroic Throughput')
   plt.ylim(-0.05,1.05)
   plt.xlabel("Wavelength [nm]")
   plt.ylabel("Throughput ratio")
   plt.show()

   plt.plot(llamas_blue.waves, bg_dichroic_b, c='b')
   plt.plot(llamas_green.waves, bg_dichroic_g, c='g')

   plt.title('BG Dichroic Throughput')
   plt.ylim(-0.05,1.05)
   plt.xlabel("Wavelength [nm]")
   plt.ylabel("Throughput ratio")
   plt.show()

   return red_med, green_med, blue_med, fiber_run_med


if __name__ == '__main__':
   plot_throughput("llamas_red.def","llamas_green.def","llamas_blue.def")