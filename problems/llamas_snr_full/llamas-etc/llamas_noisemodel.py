import spectrograph as spec
import matplotlib.pyplot as plt

llamas_blue = spec.Spectrograph('LLAMAS_BLUE')
llamas_blue.build_model('llamas_blue.def')

plt.plot(llamas_blue.waves,llamas_blue.throughput)
plt.show()

