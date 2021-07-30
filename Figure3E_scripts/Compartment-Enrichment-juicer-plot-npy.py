'''
Compartment enrichment analysis
Using juicer dump oe data 
Python ver. 2.7
Used in Miura et al., Nat.Genet. 2019 Fig3e
Plot
'''

import os
import sys
import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed
from time import time
from scipy.sparse import coo_matrix

argvs=sys.argv

#out_name=argvs[1]
#out_name+"_mean_med_std_oe_percentile_AtoB.npy"
#"../data/compartment_enrichment/CBMS1_d0_ALL_200k_newfilter_IC_cis.hic_200kb_mean_med_std_oe_percentile_AtoB.npy"
npy_file=argvs[1]
out_prefix=argvs[2]

out_name=os.path.basename(npy_file)
out_name=out_name.replace('_200kb_mean_med_std_oe_percentile_AtoB.npy', '')

out_oe_mean=np.load(npy_file)
bin1 = out_oe_mean[0,:]
bin2 = out_oe_mean[1,:]
oe_mean = out_oe_mean[2,:]
oe_med = out_oe_mean[3,:]
oe_std = out_oe_mean[4,:] 

size=100
diag= bin1 != bin2
bin1_diag = np.array(bin1[diag],dtype="float")
bin2_diag = np.array(bin2[diag],dtype="float")

coo_tri=coo_matrix((oe_mean,(bin1,bin2)), shape=(size,size))
oe_mean_diag = np.array(oe_mean[diag],dtype="float")
coo_tri_T_wo_diag = coo_matrix((oe_mean_diag,(bin1_diag,bin2_diag)), shape=(size, size)).transpose()
coo_complete = coo_tri + coo_tri_T_wo_diag
out = coo_complete.todense()
Amean = np.array(out)

coo_tri=coo_matrix((oe_med,(bin1,bin2)), shape=(size,size))
oe_mean_diag = np.array(oe_med[diag],dtype="float")
coo_tri_T_wo_diag = coo_matrix((oe_mean_diag,(bin1_diag,bin2_diag)), shape=(size, size)).transpose()
coo_complete = coo_tri + coo_tri_T_wo_diag
out = coo_complete.todense()
Amed = np.array(out)

###########################Plot compartment enrichment##################
#fig_name=out_name+'_contract_enrichment_in_compartments_scale.png'
fig_name=out_prefix+'.png'
plt.figure(figsize=(5.5,5))
#plt.title(os.path.basename(out_name)+'\nA/B enrichment')
plt.title(out_name+'\nA/B enrichment')
plt.imshow(np.log2(Amean),interpolation='none',vmin=-1.5,vmax=1.5,cmap='jet')
cb = plt.colorbar()
cb.set_label('log2 Mean oe')
plt.ylabel('Percentile AtoB')
plt.xlabel('Percentile AtoB')
plt.savefig(fig_name,dpi=100)
plt.close()

#fig_name=out_name+'_contract_enrichment_in_compartments_med_scale.png'
#plt.figure(figsize=(5.5,5))
#plt.gcf.canvas.set_window_title('test Pdzrn3 MEF WT#9 A/B enrichment')
#plt.title(os.path.basename(out_name)+'\nA/B enrichment')
#plt.imshow(np.log2(Amed),interpolation='none',vmin=-1.0,vmax=1.0,cmap='jet')
#cb = plt.colorbar()
#cb.set_label('log2 Median oe')
#plt.ylabel('Percentile AtoB')
#plt.xlabel('Percentile AtoB')
#plt.savefig(fig_name,dpi=100)
#plt.close()


