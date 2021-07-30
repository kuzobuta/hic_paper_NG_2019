'''
Compartment enrichment analysis
Using juicer dump oe data 
Python ver. 2.7
Used in Miura et al., Nat.Genet. 2019 Fig3e
'''

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed
from time import time
from scipy.sparse import coo_matrix

argvs=sys.argv

compartment_file=argvs[1]
oe_file=argvs[2]
out_name=argvs[3]

#CPU cores
#n_jobs=15
n_jobs=2

######################load compartment file##############
file=compartment_file
print("loading "+file)

data = np.loadtxt(file,delimiter="\t",skiprows=0,dtype='string')
data_woX=data[data[:,0]!="chrX",:]

percentiles=[]
for i in range(101):
   percentiles.append([i,np.percentile(np.array(data_woX[:,3],dtype='float'),i)])

percentiles=np.array(percentiles,dtype="float")
pc_percentile= np.array(data_woX[:,3],dtype="float")
tmp=np.array(data_woX[:,3],dtype="float")
for i in range(100):
   if i==99:
      pc_percentile[(tmp >= percentiles[i,1]) & (tmp <= percentiles[i+1,1])]=i+1
   else:
      pc_percentile[(tmp >= percentiles[i,1]) & (tmp < percentiles[i+1,1])]=i+1


data_woX_percentiles=np.transpose(np.vstack([np.transpose(data_woX[:,0:2]),pc_percentile]))

####################load oe file##########################
file2=oe_file
print("loading "+file2)
data_oe = np.loadtxt(file2,delimiter="\t",dtype='string')
data_oe_woX=data_oe[data_oe[:,0]!="chrX",:]

#convert pos to percentile
#Dict
chr=np.char.array(data_woX_percentiles[:,0])
pos=np.char.array(data_woX_percentiles[:,1])
sep=np.char.array(len(data_woX_percentiles)*["_"])
pos_names=chr + sep + pos
pos2percentile=dict(zip(pos_names,data_woX_percentiles[:,2]))

#oe data
chr1=np.char.array(data_oe_woX[:,0])
pos1=np.char.array(data_oe_woX[:,1])
chr2=np.char.array(data_oe_woX[:,2])
pos2=np.char.array(data_oe_woX[:,3])
sep=np.char.array(len(data_oe_woX)*["_"])

chr_pos1=chr1 + sep + pos1
chr_pos2=chr2 + sep + pos2
chr_pos1_id = np.array([pos2percentile[letter] for letter in chr_pos1],dtype="float")
chr_pos2_id = np.array([pos2percentile[letter] for letter in chr_pos2],dtype="float")
idx1= chr_pos1_id.copy()
idx1[chr_pos1_id < chr_pos2_id]=chr_pos2_id[chr_pos1_id < chr_pos2_id]
idx2= chr_pos2_id.copy()
idx2[chr_pos1_id < chr_pos2_id]=chr_pos1_id[chr_pos1_id < chr_pos2_id]

idx1=np.array(100 - idx1,dtype="int64")
idx2=np.array(100 - idx2,dtype="int64")

idx1_c=np.char.array(idx1)
idx2_c=np.char.array(idx2)
sep=np.char.array(len(idx1)*["_"])
idx12_c= idx1_c + sep + idx2_c

oe_data_woX=np.array(data_oe_woX[:,4],dtype="float")

unique_idxes=np.unique(idx12_c)

#Not use this , as it is taking too much long time..
#tmp=[];i=0
#for name in unique_idxes:
#   i=i+1
#   tmp.append([name,np.mean(oe_data_woX[idx12_c==name])])
#   print(i)

#tmp2 = [np.mean(np.array(oe_data_woX[idx12_c==name],dtype="float") for name in unique_idxes]

#def mean(name):
#  return(np.mean(oe_data_woX[idx12_c==name]))
#
#start = time()
#tmp3 = Parallel(n_jobs=-1, verbose=0)([delayed(mean)(oe_data_woX[idx12_c==name]) for name in unique_idxes]) 
#time_tmp3 = time() - start
#no parallel.. 255.15645503997803 secs, still taking long time..

print("Calculate mean, median, std oe")

def mean(name):
  return(np.nanmean(oe_data_woX[idx12_c==name]))

def median(name):
  return(np.nanmedian(oe_data_woX[idx12_c==name]))

def std(name):
  return(np.nanstd(oe_data_woX[idx12_c==name]))

start = time()
tmp_mean = Parallel(n_jobs=n_jobs, verbose=0)([delayed(mean)(name) for name in unique_idxes]) 
tmp_median = Parallel(n_jobs=n_jobs, verbose=0)([delayed(median)(name) for name in unique_idxes])
tmp_std = Parallel(n_jobs=n_jobs, verbose=0)([delayed(std)(name) for name in unique_idxes])
time_analysis = time() - start

print(str(time_analysis)+" secs")
#29.757885932922363 <= fast!!

oe_mean=np.array(tmp_mean,dtype="float").copy()
oe_median=np.array(tmp_median,dtype="float").copy()
oe_std=np.array(tmp_std,dtype="float").copy()

#Sparse to dense matrix
bin1 = np.array(unique_idxes.split("_").tolist(),dtype="int64")[:,0]
bin2 = np.array(unique_idxes.split("_").tolist(),dtype="int64")[:,1]

#out_oe_mean=np.vstack([bin1,bin2,oe_mean])
#np.save(out_name+"_mean_oe_percentile_AtoB.npy",out_oe_mean)
out_oe_mean_med_std= np.vstack([bin1,bin2,oe_mean,oe_median,oe_std])
np.save(out_name+"_mean_med_std_oe_percentile_AtoB.npy",out_oe_mean_med_std)

#size=100
#coo_tri=coo_matrix((oe_mean,(bin1,bin2)), shape=(size,size))
#diag= bin1 != bin2
#bin1_diag = np.array(bin1[diag],dtype="float")
#bin2_diag = np.array(bin2[diag],dtype="float")
#oe_mean_diag = np.array(oe_mean[diag],dtype="float")
#coo_tri_T_wo_diag = coo_matrix((oe_mean_diag,(bin1_diag,bin2_diag)), shape=(size, size)).transpose()
#coo_complete = coo_tri + coo_tri_T_wo_diag
#out = coo_complete.todense()
#A = np.array(out)

###########################Plot compartment enrichment##################
#fig_name="test_contact_enrichment.png"
#fig_name=out_name+'_contract_enrichment_in_compartments.png'
#plt.figure(figsize=(5.5,5))
#plt.gcf.canvas.set_window_title('test Pdzrn3 MEF WT#9 A/B enrichment')
#plt.title(out_name+'\n A/B enrichment')
#plt.imshow(np.log2(A),interpolation='none',cmap='jet')
#cb = plt.colorbar()
#cb.set_label('log2 Mean oe')
#plt.ylabel('Percentile AtoB')
#plt.xlabel('Percentile AtoB')
#plt.savefig(fig_name,dpi=100)





