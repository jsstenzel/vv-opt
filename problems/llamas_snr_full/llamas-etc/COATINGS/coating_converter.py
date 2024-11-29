import re

def convert(filename, doFlip, changePercent):
	wavelengths = []
	throughputs = []
	
	with open(filename, "r") as f:
		for line in f:
			words = line.split()
			print(words)
			if len(words)==0:
				continue
			if re.search('[a-zA-Z]', words[0]):
				continue
			if not re.search('[0-9]', words[0]):
				continue
			wavelengths.append(float(words[0].replace('#','')))
			throughputs.append(float(words[1]))
			
	if changePercent:
		for i,t in enumerate(throughputs):
			throughputs[i] = t/100.0
			
	if doFlip:
		for i,t in enumerate(throughputs):
			throughputs[i] = 1.0 - t
			
	with open(filename, "w") as f:
		for w,t in zip(wavelengths, throughputs):
			f.write(str(w)+'\t'+str(t)+'\n')

def dichroic(filename, split, flip=False):
	wavelengths = []
	throughputs = []
	with open(filename, "r") as f:
		for line in f:
			words = line.split()
			if len(words)==0:
				continue
			if "#" not in words[0]:
				wavelengths.append(float(words[0]))
				throughputs.append(float(words[1]))
				
	for i,pt in enumerate(wavelengths):
		if float(pt) < split or (float(pt) >= split and flip):
			throughputs[i] = 1.0 - throughputs[i]

	with open(filename+"split"+str(split), "w") as f:
		for w,t in zip(wavelengths, throughputs):
			f.write(str(w)+'\t'+str(t)+'\n')

	return wavelengths, throughputs
	

if __name__ == '__main__':  
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument('-f', metavar='string', required=True, help='File to convert')
	parser.add_argument('--flip', action='store_true', help='Change from reflectance to transmittance')
	parser.add_argument('--percent', action='store_true', help='Change total 100 to total 1')
	parser.add_argument('-s', metavar='string', required=False, help='Split dichroic at this wavelength')
	args = parser.parse_args()
	
	#We assume that we want all coating files to be throughput, out of a total of 1
	#use --flip if its a reflectance file
	#use --percent if its out of 100
	if args.s:
		dichroic(args.f, float(args.s))
	if args.flip or args.percent:
		convert(args.f, args.flip, args.percent)