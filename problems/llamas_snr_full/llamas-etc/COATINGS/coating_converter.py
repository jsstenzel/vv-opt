import sys
import re
import pandas as pd
import os
import csv

def convert(filename, doFlip, centToFrac, fracToCent, doSort):
	wavelengths = []
	throughputs = []
	
	#check file extension type
	if filename.endswith(".xlsx"):
		writefile = filename[:-5]+'.txt'
		df = pd.read_excel(filename)
		wavelengths = [float(i) for i in list(df.iloc[:, 0])]
		throughputs = [float(i) for i in list(df.iloc[:, 1])]
	elif filename.endswith(".csv"):
		writefile = filename[:-4]+'.txt'
		with open(filename) as csvfile:
			csvreader = csv.reader(csvfile, delimiter=',') #could be trouble someday
			for row in csvreader:
				if len(row)==0:
					continue
				if re.search('[a-zA-Z]', row[0]):
					continue
				if not re.search('[0-9]', row[0]):
					continue
				wavelengths.append(float(row[0].replace('#','')))
				throughputs.append(float(row[1]))
	elif filename.endswith(".txt"):
		writefile = filename
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
	else:
		print("Sorry, don't recognize that file extension.")
		sys.exit()
			
	if centToFrac:
		for i,t in enumerate(throughputs):
			throughputs[i] = t/100.0
			
	if fracToCent: #I'm probably only changing fractions to percents if i messed something up
		for i,t in enumerate(throughputs):
			throughputs[i] = t*100.0
			
	if doFlip:
		for i,t in enumerate(throughputs):
			throughputs[i] = 1.0 - t
			
	if doSort:
		wavelengths, throughputs = (list(t) for t in zip(*sorted(zip(wavelengths, throughputs))))
			
	#if extension is txt, edit it
	with open(writefile, "w+") as f:
		for w,t in zip(wavelengths, throughputs):
			f.write(str(w)+'\t'+str(t)+'\n')
	#else, make a new file with the same name and the txt extension

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
	parser.add_argument('--cent_to_frac', action='store_true', help='Change total 100 to total 1')
	parser.add_argument('--frac_to_cent', action='store_true', help='Change total 1 to total 100')
	parser.add_argument('--sort', action='store_true', help='Sort the data from shortest to longest wavelength')
	parser.add_argument('-s', metavar='string', required=False, help='Split dichroic at this wavelength')
	args = parser.parse_args()
	
	#We assume that we want all coating files to be throughput, out of a total of 1
	#use --flip if its a reflectance file
	#use --percent if its out of 100 and we want to turn it to a fraction
	if args.s:
		dichroic(args.f, float(args.s))
	if args.flip or args.cent_to_frac or args.frac_to_cent or args.sort:
		convert(args.f, args.flip, args.cent_to_frac, args.frac_to_cent, args.sort)
	#else:
	#	print("What's that?")