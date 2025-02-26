import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from pandas.plotting import scatter_matrix

def plotmatrix(matrix, names=[], xlim=[], ylim=[], title='', c='#21aad3', rescale=False):
	if len(names)==0:
		names = ["index"+str(i) for i in range(len(matrix))]
	
	data_df = pd.DataFrame(matrix, columns=names)
	
	if rescale:
		data_df = (data_df - data_df.mean())/data_df.std()
	
	hist_kwds = {'color':c}
	axarr = scatter_matrix(data_df, hist_kwds=hist_kwds, alpha=0.2, figsize=(7, 6), diagonal='hist', c=c)
	
	#set all of the subfigures to match each other
	for i in range(len(names)):
		for j in range(len(names)):
			axarr[i,j].set_xlim(xlim[0],xlim[1]) if len(xlim)>0 else 0
			if i != j:
				axarr[i,j].set_ylim(ylim[0],ylim[1]) if len(ylim)>0 else 0
		
	plt.suptitle(title)
	plt.show()
	plt.close()

def covmatrix_heatmap(matrix, names, rescale=True, doPlot=True):
	if len(names)==0:
		names = ["index"+str(i) for i in range(len(matrix))]
	
	data_df = pd.DataFrame(matrix, columns=names)
	
	if rescale:
		data_df = (data_df - data_df.mean())/data_df.std()

	# Compute the covariance matrix
	cov_matrix = data_df.cov()

	# Create a heatmap using Seaborn
	if doPlot:
		sns.heatmap(cov_matrix, annot=False, xticklabels=True, yticklabels=True, cmap='coolwarm')
		plt.show()
	
	return cov_matrix
	
def cov_heatmap(cov_matrix, names, rescale=True, doPlot=True):
	# Create a heatmap using Seaborn
	if doPlot:
		sns.heatmap(cov_matrix, annot=False, xticklabels=names, yticklabels=names, cmap='coolwarm')
		plt.show()
	
	return cov_matrix