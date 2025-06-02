#k-clustering problem for large datasets, and evaluating their BIC and AIC
#This is an attempt to find a good initial ncomp for BN modeling

import os
import sys
import math
import csv
import random
import matplotlib.pyplot as plt
import numpy as np
from sklearn.cluster import MiniBatchKMeans

import matplotlib.cm as cm
from sklearn.cluster import KMeans
from sklearn.datasets import make_blobs
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.metrics import calinski_harabasz_score

sys.path.append('..')
#from inference.goal_based_inference import *
				

#Read and interpret the savefile
def cluster_load_samples(savefile, expected_length, doPrint=False, do_subset=0, doDiagnostic=True):
	filename = savefile if savefile.endswith('.csv') else savefile+'.csv'

	###Make sure the files exist
	if not os.path.isfile(filename):
		print("File",filename,"is missing")
		sys.exit()
		
	###Safely read out all of the samples into matrices
	qyd = []
	
	if doPrint:
		print("Reading the data file",savefile,"...",flush=True)
	
	with open(filename) as csvfile:
		csvreader = csv.reader((line.replace('\0','') for line in csvfile ), delimiter=',')
		for l,row in enumerate(csvreader):
			if len(row) != expected_length:
				if doDiagnostic:
					print("Warning: dropped line",l+1,"(length "+str(len(row))+' expected', str(expected_length)+')',"from",filename)
			elif do_subset == 0 or len(qyd) < do_subset:
				qyd_line = [float(e) for e in row]
				qyd.append(qyd_line)
			else:
				break
	
	if do_subset:
		qyd = qyd[:do_subset]
		
	###Standardize
	data_mean = np.mean(qyd, axis=0)
	data_std = np.std(qyd, axis=0)
	standard_data = [list((qyd_i - data_mean)/data_std) for qyd_i in qyd]

	###Shuffle it once, for good measure
	#random.shuffle(standard_data)

	return standard_data
		
def cluster_criterion_plot(ks, BICs, AICs, inertias, scores):
	###Print, plot
	fig, (ax1, ax2) = plt.subplots(1,2)

	# Plot the first data set on the left y-axis
	"""
	ax1.plot(ks, BICs, color='blue', label='BIC')
	ax1.plot(ks, AICs, color='orange', label='AIC')
	ax1.set_xlabel("Number of clusters")
	ax1.set_xticks(ks)
	ax1.set_ylabel('Information criterion')
	ax1.legend()
	"""
	ax1.plot(ks, scores, color='purple')#, label='BIC')
	ax1.set_xlabel("Number of clusters")
	ax1.set_xticks(ks)
	ax1.set_ylabel('Calinski-Harabasz score')
	#ax1.legend()
	
	ax2.set_xlabel("Number of clusters")
	ax2.set_xticks(ks)
	ax2.plot(ks, inertias, color='green')#, label='inertia')
	ax2.tick_params(axis='y')
	ax2.set_ylabel('Inertia')
	#ax2.legend()

	plt.show()

#NOTE that this AIC and BIC are based on k = # components
#This is different from the use case of training a GMM, which has:
#k = ncomp*((D*D - D)/2 + 2D + 1)-1 for a full covariance matrix
#therefore, both AIC and BIC will have a weaker penality,
#and therefore will predict k_cluster > n_comp
#thus, this method gives kind of an upper bound of what GMM training will show
def cluster_calc_criteria(score, k, d, N):
	#https://stats.stackexchange.com/questions/55147/compute-bic-clustering-criterion-to-validate-clusters-after-k-means
	BIC = -2 * score + k * d * np.log(N)
	AIC = -2 * score + 2*d*k
	return score, BIC, AIC

###Perform k-clustering for k, calculate score, BIC, AIC
def cluster_eval_k(data, k, doPrint=True):
	###Shuffle data, for good measure
	#if doPrint:
	#	print("Shuffling data...",flush=True)
	#random.shuffle(standard_data)
	N = len(data)
	d = len(data[0])
	
	###run kmeans fit
	if doPrint:
		print("Solving k="+str(k),"k-means problem, N="+str(N)+" d="+str(d),"...",flush=True)
	kmeans = MiniBatchKMeans(n_clusters=k, batch_size=1024, n_init='auto')
	kmeans.fit(data)
	
	###get the score
	if doPrint:
		print("Scoring the kmeans...",flush=True)
	scores = kmeans.score(data)
	score = np.mean(scores)
	
	CH_score = calinski_harabasz_score(data, kmeans.labels_)
	
	###return BIC, AIC
	_, BIC, AIC = cluster_calc_criteria(score, k, d, N)
	return BIC, AIC, kmeans.inertia_, CH_score

###run the whole problem for multiple k, then plot
def cluster_eval(savefile, ks, expected_length, do_subset=0, doPrint=True):
	if doPrint:
		print("Getting data from",savefile,'...',flush=True)
	data = cluster_load_samples(savefile, expected_length=expected_length, doPrint=doPrint, do_subset=do_subset, doDiagnostic=False)
	BICs = [None]*len(ks)
	AICs = [None]*len(ks)
	inertias = [None]*len(ks)
	scores = [None]*len(ks)
	
	for i,k in enumerate(ks):
		if doPrint:
			print("Evaluating k="+str(k),'...',flush=True)
		BIC, AIC, inertia, CH_score = cluster_eval_k(data, k=k, doPrint=False)
		print(BIC, AIC, inertia, CH_score)
		BICs[i] = BIC
		AICs[i] = AIC
		inertias[i] = inertia
		scores[i] = CH_score
		
	cluster_criterion_plot(ks, BICs, AICs, inertias, scores)
	
def silhouette_scoring(savefile, n_clusters=2, expected_length=0, do_subset=0, doPrint=True):
	if doPrint:
		print("Getting data from",savefile,'...',flush=True)
	data = cluster_load_samples(savefile, expected_length=expected_length, doPrint=doPrint, do_subset=do_subset, doDiagnostic=False)
	X = np.array(data)

	fig, ax1 = plt.subplots(1, 1)
	fig.set_size_inches(18, 7)
	# The 1st subplot is the silhouette plot
	# The silhouette coefficient can range from -1, 1 but in this example all
	# lie within [-0.1, 1]
	ax1.set_xlim([-0.1, 1])
	# The (n_clusters+1)*10 is for inserting blank space between silhouette
	# plots of individual clusters, to demarcate them clearly.
	ax1.set_ylim([0, len(X) + (n_clusters + 1) * 10])

	# Initialize the clusterer with n_clusters value and a random generator
	# seed of 10 for reproducibility.
	clusterer = KMeans(n_clusters=n_clusters, random_state=10)
	cluster_labels = clusterer.fit_predict(X)

	# The silhouette_score gives the average value for all the samples.
	# This gives a perspective into the density and separation of the formed
	# clusters
	silhouette_avg = silhouette_score(X, cluster_labels)
	print(
		"For n_clusters =",
		n_clusters,
		"The average silhouette_score is :",
		silhouette_avg,
	)

	# Compute the silhouette scores for each sample
	sample_silhouette_values = silhouette_samples(X, cluster_labels)

	y_lower = 10
	
	for i in range(n_clusters):
		# Aggregate the silhouette scores for samples belonging to
		# cluster i, and sort them
		ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]

		ith_cluster_silhouette_values.sort()

		size_cluster_i = ith_cluster_silhouette_values.shape[0]
		y_upper = y_lower + size_cluster_i

		color = cm.nipy_spectral(float(i) / n_clusters)
		ax1.fill_betweenx(
			np.arange(y_lower, y_upper),
			0,
			ith_cluster_silhouette_values,
			facecolor=color,
			edgecolor=color,
			alpha=0.7,
		)

		# Label the silhouette plots with their cluster numbers at the middle
		ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

		# Compute the new y_lower for next plot
		y_lower = y_upper + 10  # 10 for the 0 samples
		
	plt.show()

if __name__ == '__main__':  
	ks = [k for k in range(2,100)]
	cluster_eval("BN_batch_samples", ks, expected_length=10, do_subset=4000000, doPrint=True)
	#for n in [2,3,4,5,6,7,8,9,10]:
	#	silhouette_scoring("BN_batch_samples", n_clusters=n, do_subset=4000000, expected_length=10)
	
"""
Evaluating k=10 ...
200384042.1432969 200376221.34100923
Evaluating k=20 ...
161989230.1230605 161973588.51848516
Evaluating k=30 ...
150768602.70796314 150745140.30110016
Evaluating k=40 ...
152850332.7854778 152819049.57632715
Evaluating k=50 ...
135162092.98006383 135122988.96862552
Evaluating k=60 ...
134048402.48227558 134001477.66854961
Evaluating k=70 ...
143185541.11812785 143130795.50211424
Evaluating k=80 ...
142398170.55252546 142335604.13422418
Evaluating k=90 ...
130200972.76842284 130130585.54783389
Evaluating k=100 ...
125958216.02703139 125880008.00415477
Evaluating k=110 ...
125516120.45433122 125430091.62916695
Evaluating k=120 ...
123927957.07727835 123834107.4498264
Evaluating k=130 ...
137934029.36595672 137832358.93621713
Evaluating k=140 ...
122000877.04480289 121891385.81277563
Evaluating k=150 ...
119884316.94874781 119767004.9144329
Evaluating k=160 ...
121133754.54135247 121008621.70474988
Evaluating k=170 ...
118504242.14559141 118371288.50670117
Evaluating k=180 ...
118019291.20492493 117878516.76374702
Evaluating k=190 ...
121372741.30634673 121224146.06288116
Evaluating k=200 ...
119307504.48440887 119151088.43865564
Evaluating k=210 ...
116648197.16092502 116483960.31288412
Evaluating k=220 ...
119280541.09095117 119108483.44062263
Evaluating k=230 ...
115179568.2649243 114999689.81230809
Evaluating k=240 ...
114940314.28634125 114752615.03143738
Evaluating k=250 ...
115117577.46647623 114922057.4092847
Evaluating k=260 ...
114049756.26420797 113846415.40472879
Evaluating k=270 ...
115495209.86147773 115284048.19971088
Evaluating k=280 ...
112909903.74370465 112690921.27965012
Evaluating k=290 ...
113282132.8470477 113055329.58070552
Evaluating k=300 ...
113308553.27894706 113073929.21031721
Evaluating k=310 ...
111371719.29277901 111129274.42186151
Evaluating k=320 ...
113469454.78278501 113219189.10957985
Evaluating k=330 ...
118830773.3654894 118572686.88999657
Evaluating k=340 ...
112206316.1759645 111940408.89818402
Evaluating k=350 ...
111570697.38556847 111296969.30550033
Evaluating k=360 ...
111768412.6615825 111486863.77922669
Evaluating k=370 ...
114853562.88038915 114564193.19574569
Evaluating k=380 ...
109861728.92770177 109564538.44077064
Evaluating k=390 ...
117416469.26059839 117111457.9713796
Evaluating k=400 ...
110751825.68353543 110438993.59202898
Evaluating k=410 ...
110050417.6105072 109729764.71671309
Evaluating k=420 ...
110115357.60375285 109786883.90767108
Evaluating k=430 ...
109235924.68739602 108899630.18902658
Evaluating k=440 ...
109736038.69868448 109391923.39802739
Evaluating k=450 ...
108692604.38397625 108340668.28103149
Evaluating k=460 ...
109771509.10575788 109411752.20052546
Evaluating k=470 ...
110364339.60742836 109996761.89990827
Evaluating k=480 ...
111621152.40737076 111245753.89756301
Evaluating k=490 ...
108909842.39686699 108526623.08477159
Evaluating k=500 ...
109569365.23784009 109178325.12345701
Evaluating k=510 ...
122063030.77781169 121664169.86114097
Evaluating k=520 ...
109046493.81426682 108639812.09530842
Evaluating k=530 ...
107324905.57248928 106910403.05124323
Evaluating k=540 ...
109655285.1823528 109232961.85881908
Evaluating k=550 ...
107308635.67760715 106878491.55178578
Evaluating k=560 ...
108714817.80971321 108276852.88160418
Evaluating k=570 ...
108511715.71374348 108065929.98334678
Evaluating k=580 ...
106805196.76251014 106351590.22982578
Evaluating k=590 ...
108074074.6572342 107612647.32226218
Evaluating k=600 ...
106909393.52704236 106440145.38978268
Evaluating k=610 ...
106263008.469842 105785939.53029466
Evaluating k=620 ...
118071003.39072834 117586113.64889334
Evaluating k=630 ...
107767983.52578565 107275272.98166299
Evaluating k=640 ...
105095157.09122103 104594625.74481072
Evaluating k=650 ...
106533092.14475022 106024739.99605224
Evaluating k=660 ...
105276827.84358107 104760654.89259541
Evaluating k=670 ...
105092685.37611207 104568691.62283877
Evaluating k=680 ...
105140917.43058422 104609102.87502325
Evaluating k=690 ...
104967413.66447958 104427778.30663095
Evaluating k=700 ...
107748601.46677352 107201145.30663723
Evaluating k=710 ...
105459191.0672189 104903914.10479495
Evaluating k=720 ...
105644189.94202888 105081092.17731726
Evaluating k=730 ...
106936599.38343282 106365680.81643355
Evaluating k=740 ...
104881292.97622484 104302553.6069379
Evaluating k=750 ...
104455317.52222401 103868757.3506494
Evaluating k=760 ...
104461243.81166908 103866862.83780682
Evaluating k=770 ...
105142297.50859708 104540095.73244715
Evaluating k=780 ...
106279407.6512943 105669385.07285672
Evaluating k=790 ...
107276640.02048624 106658796.63976099
Evaluating k=800 ...
105272684.32110542 104647020.13809252
Evaluating k=810 ...
103988453.74193689 103354968.75663632
Evaluating k=820 ...
104717595.59073596 104076289.80314773
Evaluating k=830 ...
104164916.77660055 103515790.18672466
Evaluating k=840 ...
106874168.45268466 106217221.0605211
Evaluating k=850 ...
105322608.23274502 104657840.03829381
Evaluating k=860 ...
103734038.1829721 103061449.18623322
Evaluating k=870 ...
103084290.54731074 102403880.7482842
Evaluating k=880 ...
103908767.64093092 103220537.03961672
Evaluating k=890 ...
104679270.52367435 103983219.1200725
Evaluating k=900 ...
104595542.34122561 103891670.13533609
Evaluating k=910 ...
103149258.13830945 102437565.13013227
Evaluating k=920 ...
103589195.34642212 102869681.53595728
Evaluating k=930 ...
102873080.24179105 102145745.62903854
Evaluating k=940 ...
103360004.21311066 102624848.79807049
Evaluating k=950 ...
103170174.66569968 102427198.44837184
Evaluating k=960 ...
103923026.72601883 103172229.70640334
Evaluating k=970 ...
102336949.47354756 101578331.65164441
Evaluating k=980 ...
103838612.97167516 103072174.34748435
Evaluating k=990 ...
102757473.63606623 101983214.20958775
"""
"""
Evaluating k=1000 ...
111944000.84412444 111156953.6503441
Evaluating k=2000 ...
106127102.15989754 104553007.77233687
Evaluating k=3000 ...
102377360.92028879 100016219.3389478
Evaluating k=4000 ...
102083722.82026845 98935534.04514714
Evaluating k=5000 ...
100895893.96065024 96960657.99174859
Evaluating k=6000 ...
100073777.71354988 95351494.5508679
Evaluating k=7000 ...
100302849.66980253 94793519.31334022
Evaluating k=8000 ...
98657685.95286953 92361308.4026269
Evaluating k=9000 ...
98804465.38818447 91721040.64416151
"""