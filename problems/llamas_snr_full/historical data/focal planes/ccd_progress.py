"""

"""

import sys
import os
import matplotlib
import matplotlib.pyplot as plt

req_rn = 2.5
req_dark = 0.001
req_gain = 1.0
req_gainmin = 1.0 - 0.2
req_gainmax = 1.0 + 0.2
compare = True
squish = True
uncertainty = False
oneplot = True

#camera, read noise, dark current, gain
ccd_data = [
	["SN20001", "Red Qual", "nodither", [2.32,.0237], [0.00238,.00015], [0.999,.0102]],
	["SN20003", "Blue Qual", "dither", [2.35,.0243], [0.00267,5.27e-05], [1.088,.0112]],
	#
	["SN20004", "Blue 1", "dither", [2.31,0.0236], [0.000865,3.16e-05], [0.998,.0102]],
	["SN20005", "Green 1", "dither", [2.4,.0246], [0.00196,.000338], [1.037,.0106]],
	#
	["SN20006", "Red 1", "nodither", [2.38,.0245], [0.000976,.000394], [1.059,.0109]],
	["SN20007", "Blue 2", "dither", [2.26,.0235], [0.00176,2.91e-05], [1.094,.0114]],
	["SN20008", "Green 2", "dither", [2.77,.0281], [0.00198,4.47e-05], [1.010,.0103]],
	#
	["SN20009", "Red 2", "nodither", [2.56,.0261], [0.00113,.000272], [0.903,.0092]],
	["SN20010", "Blue 3", "dither", [2.44,.0253], [0.00105,2.91e-05], [1.033,.0107]],
	["SN20011", "Green 3", "dither", [2.69,.0223], [0.00147,5.17e-05], [1.061,.0088]],
	#
	["SN20012", "Red 3", "nodither", [2.44,.025], [0.00113,.000248], [1.031,.0106]],
	["SN20013", "Blue 4", "dither", [2.44,.0211], [0.00107,4.9e-05], [1.073,.0086]],
	["SN20014", "Green 4", "dither", [2.54,0.0262], [0.00147,0.000131], [1.061,0.0109],
	#
	["SN20015", "Red 4", "nodither", [2.35,.0241], [0.00095,.000202], [1.007,.0103]],
	["SN20016", "Blue 5", "dither", [2.65,.0271], [0.00107,4.87e-05], [1.053,.0108]],
	["SN20017", "Green 5", "dither", [2.41,.0247], [0.00117,.000496], [1.029,.0105]],
	#
	["SN20018", "Red 5", "nodither", [2.43,.0248], [0.00191,.000131], [1.051,.0107]],
	#
	["SN20019", "Blue 6", "dither", [2.26,.023], [0.00054,2.34e-05], [1.025,.0104]],
	["SN20020", "Green 6", "dither", [2.66,.0273], [0.00121,.00024], [1.063,.0109]],
	["SN20021", "Red 6", "nodither", [2.4,.0246], [0.00117,.000159], [1.069,.0109]],
	#
	["SN20022", "Blue 7", "dither", [2.59,.0267], [0.000689,2.8e-05], [1.071,.011]],
	["SN20023", "Green 7", "dither", [2.69,.0281], [0.000643,8.05e-05], [1.007,.0105]],
	["SN20024", "Red 7", "nodither", [2.63,.0268], [0.00171,.000112], [1.044,.0106]],
	#
	["SN20025", "Green 8", "dither", [2.47,.0251], [0.00137,4.12e-05], [1.049,.0107]],
	["SN20026", "Red 8", "dither", [2.48,.0253], [0.00076,.000212], [1.032,.0105]],
	["SN20002", "Blue 8", "dither", [2.48,0.0256], [0.00137,8.18E-05], [1.073,0.0111]]
	]
	
perf_data = [
	["SN20001", "Red Qual", "nodither", 2.62, 0.0018, 1.12],
	["SN20003", "Blue Qual", "dither", 2.41, 0.0046, 1.13],
	#
	["SN20004", "Blue 1", "dither", 2.4, 0.0009, 1.05],
	["SN20005", "Green 1", "dither", 2.5, 0.00057, 1.08],
	#
	["SN20006", "Red 1", "nodither", 2.7, 0.0010, 1.16],
	["SN20007", "Blue 2", "dither", 2.3, 0.0009, 1.12],
	["SN20008", "Green 2", "dither", 2.9, 0.0029, 1.05],
	#
	["SN20009", "Red 2", "nodither", 2.9, 0.0022, 1.00],
	["SN20010", "Blue 3", "dither", 2.5, 0.0013, 1.08],
	["SN20011", "Green 3", "dither", 2.7, 0.0019, 1.06],
	#
	["SN20012", "Red 3", "nodither", 2.6, 0.0018, 1.13],
	["SN20013", "Blue 4", "dither", 2.6, 0.0018, 1.07],
	#["SN20014", "Green 4", "dither", 2.7, 0.0019, 1.06],
	#
	["SN20015", "Red 4", "nodither", 2.6, 0.0018, 1.11],
	["SN20016", "Blue 5", "dither", 2.8, 0.0006, 1.11],
	["SN20017", "Green 5", "dither", 2.7, 0.0019, 1.06],
	#
	["SN20018", "Red 5", "nodither", 2.8, 0.0027, 1.16],
	#
	["SN20019", "Blue 6", "dither", 2.3, 0.0008, 1.06],
	["SN20020", "Green 6", "dither", 2.8, 0.0014, 1.11],
	["SN20021", "Red 6", "nodither", 2.7, 0.0020, 1.17],
	#
	["SN20022", "Blue 7", "dither", 2.7, 0.001, 1.11],
	["SN20023", "Green 7", "dither", 2.8, 0.0007, 1.06],
	#["SN20024", "Red 7", "nodither", 2.7, 0.0020, 1.17],
	#
	["SN20025", "Green 8", "dither", 2.5, 0.0016, 1.089]
	["SN20026", "Red 8", "dither", 2.8, 0.0019, 1.16]
	#["SN20002", "Blue 8", "dither", [], [], []]
	]
	
names = []
numbers = []
num = 0
rns = []
dcs = []
gns = []
rns_err = []
dcs_err = []
gns_err = []
colors = []
neutral_color = 'gray'
for ccd in ccd_data:
	if squish == False:
		names.append(ccd[0] + '\n(' + ccd[1] + ')')
	else:
		names.append(ccd[0][-2:])	
	numbers.append(num)
	num += 1
	rns.append(ccd[3][0])
	dcs.append(ccd[4][0])
	gns.append(ccd[5][0])
	
	if uncertainty == True:
		rns_err.append(ccd[3][1])
		dcs_err.append(ccd[4][1])
		gns_err.append(ccd[5][1])

	if 'Blue' in ccd[1]:
		color = 'b'
	elif 'Green' in ccd[1]:
		color = 'g'
	elif 'Red' in ccd[1]:
		color = 'r'
	else:
		color = 'k'
	colors.append(color)
	

perf_rns = []
perf_dcs = []
perf_gns = []
for ccd in perf_data:
	perf_rns.append(ccd[3])
	perf_dcs.append(ccd[4])
	perf_gns.append(ccd[5])

if oneplot:
    #Plot Read Noise
    fig, (ax1, ax2, ax3) = plt.subplots(3)
    ax1.set_title('CCD Read Noise @ -90 C, high gain, 75kHz')
    ax1.set_ylabel('Read Noise (e- RMS)')
    ax1.set_xlabel('')
    ax1.set_xticks([n for n in range(len(numbers))])
    ax1.set_xlim(-0.5, len(numbers)+1)
    #ax1.set_ylim(1,1e6)
    ax1.scatter(numbers, rns, c=colors)
    if uncertainty:
    	ax1.errorbar(numbers, rns, yerr=rns_err, fmt='o', ecolor='k', markersize=0)
    if compare:
    	ax1.scatter(numbers, perf_rns, c=neutral_color)
    ax1.grid()
    ax1.axhline(req_rn, c='purple')
    	
    fig.canvas.draw()
    labels = [item.get_text() for item in ax1.get_xticklabels()]
    for idx,label in enumerate(labels):
    	labels[idx] = ""#names[idx]
    ax1.set_xticklabels(labels)
    	
    #plt.show()
    
    #Plot Dark Current
    #fig,ax = plt.subplots()
    ax2.set_title('CCD Dark Current @ -90 C, high gain, 75kHz')
    ax2.set_ylabel('Dark Current (e-/s/px)')
    ax2.set_xlabel('')
    ax2.set_xticks([n for n in range(len(numbers))])
    ax2.set_xlim(-0.5, len(numbers)+1)
    #ax2.set_ylim(1,1e6)
    ax2.scatter(numbers, dcs, c=colors)
    if uncertainty:
    	ax2.errorbar(numbers, dcs, yerr=dcs_err, fmt='o', ecolor='k', alpha=1.0, markersize=0)
    if compare:
    	ax2.scatter(numbers, perf_dcs, c=neutral_color)
    ax2.grid()
    ax2.axhline(req_dark, c='purple')
    
    fig.canvas.draw()
    labels = [item.get_text() for item in ax2.get_xticklabels()]
    for idx,label in enumerate(labels):
    	labels[idx] = ""#names[idx]
    ax2.set_xticklabels(labels)
    	
    #plt.show()
    
    #Plot Gain
    #fig,ax = plt.subplots()
    ax3.set_title('CCD Gain @ -90 C, high gain, 75kHz')
    ax3.set_ylabel('Gain (e-/ADU)')
    ax3.set_xlabel('')
    ax3.set_xticks([n for n in range(len(numbers))])
    ax3.set_xlim(-0.5, len(numbers)+1)
    #ax3.set_ylim(1,1e6)
    ax3.scatter(numbers, gns, c=colors)
    if uncertainty:
    	ax3.errorbar(numbers, gns, yerr=gns_err, fmt='o', ecolor='k', alpha=1.0, markersize=0)
    if compare:
    	ax3.scatter(numbers, perf_gns, c=neutral_color)
    ax3.grid()
    ax3.axhline(req_gainmin, c='purple', linestyle = 'dashed')
    ax3.axhline(req_gain, c='purple')
    ax3.axhline(req_gainmax, c='purple', linestyle = 'dashed')
    
    fig.canvas.draw()
    labels = [item.get_text() for item in ax3.get_xticklabels()]
    for idx,label in enumerate(labels):
    	labels[idx] = names[idx]
    ax3.set_xticklabels(labels)
    	
    plt.show()
    
else:
    #Plot Read Noise
    fig, ax1 = plt.subplots()
    ax1.set_title('CCD Read Noise @ -90 C, high gain, 75kHz')
    ax1.set_ylabel('Read Noise (e- RMS)')
    ax1.set_xlabel('')
    ax1.set_xticks([n for n in range(len(numbers))])
    ax1.set_xlim(-0.5, len(numbers)+1)
    #ax1.set_ylim(1,1e6)
    ax1.scatter(numbers, rns, c=colors)
    if uncertainty:
    	ax1.errorbar(numbers, rns, yerr=rns_err, fmt='o', ecolor='k', markersize=0)
    if compare:
    	ax1.scatter(numbers, perf_rns, c=neutral_color)
    ax1.grid()
    ax1.axhline(req_rn, c='purple')
    	
    fig.canvas.draw()
    labels = [item.get_text() for item in ax1.get_xticklabels()]
    for idx,label in enumerate(labels):
    	labels[idx] = names[idx]
    ax1.set_xticklabels(labels)
    	
    plt.show()
    
    #Plot Dark Current
    fig, ax2 = plt.subplots()
    ax2.set_title('CCD Dark Current @ -90 C, high gain, 75kHz')
    ax2.set_ylabel('Dark Current (e-/s/px)')
    ax2.set_xlabel('')
    ax2.set_xticks([n for n in range(len(numbers))])
    ax2.set_xlim(-0.5, len(numbers)+1)
    #ax2.set_ylim(1,1e6)
    ax2.scatter(numbers, dcs, c=colors)
    if uncertainty:
    	ax2.errorbar(numbers, dcs, yerr=dcs_err, fmt='o', ecolor='k', alpha=1.0, markersize=0)
    if compare:
    	ax2.scatter(numbers, perf_dcs, c=neutral_color)
    ax2.grid()
    ax2.axhline(req_dark, c='purple')
    
    fig.canvas.draw()
    labels = [item.get_text() for item in ax2.get_xticklabels()]
    for idx,label in enumerate(labels):
    	labels[idx] = names[idx]
    ax2.set_xticklabels(labels)
    	
    plt.show()
    
    #Plot Gain
    fig, ax3 = plt.subplots()
    ax3.set_title('CCD Gain @ -90 C, high gain, 75kHz')
    ax3.set_ylabel('Gain (e-/ADU)')
    ax3.set_xlabel('')
    ax3.set_xticks([n for n in range(len(numbers))])
    ax3.set_xlim(-0.5, len(numbers)+1)
    #ax3.set_ylim(1,1e6)
    ax3.scatter(numbers, gns, c=colors)
    if uncertainty:
    	ax3.errorbar(numbers, gns, yerr=gns_err, fmt='o', ecolor='k', alpha=1.0, markersize=0)
    if compare:
    	ax3.scatter(numbers, perf_gns, c=neutral_color)
    ax3.grid()
    ax3.axhline(req_gainmin, c='purple', linestyle = 'dashed')
    ax3.axhline(req_gain, c='purple')
    ax3.axhline(req_gainmax, c='purple', linestyle = 'dashed')
    
    fig.canvas.draw()
    labels = [item.get_text() for item in ax3.get_xticklabels()]
    for idx,label in enumerate(labels):
    	labels[idx] = names[idx]
    ax3.set_xticklabels(labels)
    	
    plt.show()