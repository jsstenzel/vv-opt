#import argparse
import os
import sys
import shutil
import csv
import fileinput
sys.path.insert(0, "..")

import matplotlib.pyplot as plt
from astropy.io import fits

import numpy as np
import multiprocessing as mp
from SALib.analyze import sobol, fast
import scipy.optimize as optimization

num_workers = mp.cpu_count() 


def sobol_saltelli_file(datafile, var_names, doSijCalc=False, doPlot=False, doPrint=False, doReport=''):
    y = []
    x = []
    with open(datafile) as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',')
        for row in csvreader:
            x.append([float(xi) for xi in row[0:-1]])
            y.append(float(row[-1]))
    #print(x)
    #print(y)

    sobol_saltelli(x, var_names, y, doSijCalc, doPlot, doPrint, doReport)


def sobol_saltelli(var_list, var_names, y_list, doSijCalc=False, doPlot=False, doPrint=False, doReport=''):
    if len(var_list) != len(y_list):
        print("ERROR: sobol_saltelli input problem: mismatching number of model evals")
        print("Length var_list:",len(var_list))
        print("Length y_list:",len(y_list))
        exit()
        
    if len(var_list[0]) != len(var_names):
        print("ERROR: sobol_saltelli input problem: mismatching number of variables")
        print("Length var_list:",len(var_list))
        print("Length var_names:",len(var_names))
        
    X = np.array(var_list)
    Y = np.array(y_list)

    dim = len(var_names)
    problem = {
        'names': var_names,
        'num_vars': dim,
        'bounds': [[0, 1]]*dim,
        'dists': ['norm']*dim
    }

    ###calculate indices
    conf = 0.95
    Si = sobol.analyze(problem, Y, calc_second_order=doSijCalc, conf_level=conf, print_to_console=doPrint, parallel=False, n_processors=None) 
    if doPrint:
        print("Confidence levels on each parameter calculated at ",conf)
    
    if doPlot:
        Si.plot() #wtf this does nothinG? Why do i have it
    
    if doReport != '':
        with open(doReport+'.csv', 'w+', newline='') as csvfile:
            writer = csv.writer(csvfile)
            for i,_ in enumerate(Si['S1']):
                row = [Si['S1'][i], Si['S1_conf'][i], Si['ST'][i], Si['ST_conf'][i]]
                if doSijCalc:
                    row.append([Si['S2'][i], Si['S2_conf'][i]])
                writer.writerow(row)

"""
def plot_gsa_WIP(filename):
    S1, S1_conf, ST, ST_conf = [],[],[],[]
    with open(filename) as csvfile:
        csvreader = csv.reader(csvfile, delimiter=',')
        for row in csvreader:
            S1.append(float(row[0]))
            S1_conf.append(float(row[1]))
            ST.append(float(row[2]))
            ST_conf.append(float(row[3]))
            
    varnames = ['red\nrn', 'red\ndc', 'green\nrn', 'green\ndc', 'blue\nrn', 'blue\ndc', 'frd', 'bg\nbias', 'sl\nbias', 'vph1\nred', 'vph2\nred', 'vph3\nred', 'vph1\ngreen', 'vph2\ngreen', 'vph3\ngreen', 'vph1\nblue', 'vph2\nblue', 'vph3\nblue']
    
    # set width of bar
    barWidth = 0.1
    fig = plt.subplots(figsize =(12, 8))
     
    # set height of bar    
    ii = 0
    bar_heights = [[],[],[],[],[],[],[],[]]
    for spect in bar_heights:
        for var in varnames:
            spect.append(ST[ii])
            ii += 1
     
    # Set position of bar on X axis
    bar_pos = [[],[],[],[],[],[],[],[]]
    bar_pos[0] = np.arange(len(varnames))
    bar_pos[1] = [x + barWidth for x in bar_pos[0]]
    bar_pos[2] = [x + barWidth for x in bar_pos[1]]
    bar_pos[3] = [x + barWidth for x in bar_pos[2]]
    bar_pos[4] = [x + barWidth for x in bar_pos[3]]
    bar_pos[5] = [x + barWidth for x in bar_pos[4]]
    bar_pos[6] = [x + barWidth for x in bar_pos[5]]
    bar_pos[7] = [x + barWidth for x in bar_pos[6]]
    
    colors = ["#005f73","#0a9396","#94d2bd","#e9d8a6","#ee9b00","#ca6702","#bb3e03","#9b2226"]
     
    # Make the plot
    #plt.bar(br1, IT, color ='r', width = barWidth, edgecolor ='grey', label ='IT')
    #plt.bar(br2, ECE, color ='g', width = barWidth, edgecolor ='grey', label ='ECE')
    #plt.bar(br3, CSE, color ='b', width = barWidth, edgecolor ='grey', label ='CSE')
    for ii,_ in enumerate(bar_heights):
        plt.bar(bar_pos[ii], bar_heights[ii], color=colors[ii], width = barWidth, edgecolor='grey', label = str(ii))
     
    # Adding Xticks
    plt.xlabel('Parameters', fontweight ='bold', fontsize = 15)
    plt.ylabel('S_T', fontweight ='bold', fontsize = 15)
    plt.xticks([r + barWidth for r in range(len(varnames))], varnames)
    #plt.tight_layout()
    plt.yscale('log')    
    plt.show()
"""
    
    
if __name__ == '__main__':  
    var_names = ['x1','x2','x3']
    sobol_saltelli_file('testdata.csv', var_names, doSijCalc=True, doPlot=True, doPrint=True, doReport='testreport')