import numpy as np
from sklearn.decomposition import PCA,TruncatedSVD
from sklearn.neighbors import NearestNeighbors
import networkx as nx
from fa2 import ForceAtlas2
import matplotlib.pyplot as plt

"""
SPRING like run using scRepliseq data
"""

fname="data/scRepliseq_884cells_Sphase_G1_data_set_filtered.txt"
E=np.loadtxt(fname, dtype="float", delimiter="\t",skiprows=1) 
Z = E.T
pca = PCA()
pca.fit(Z)
Epca=pca.transform(Z)

#Use all PCA components
X=Epca.copy()

#default
k=5
dist_metric="euclidean"
nbrs = NearestNeighbors(n_neighbors=k, metric=dist_metric).fit(X)

knn = nbrs.kneighbors(return_distance=False)
links = set([])
for i in range(knn.shape[0]):
	for j in knn[i,:]:
		links.add(tuple(sorted((i,j))))

knn_graph=knn.copy()
links=list(links)
G = nx.Graph()
G.add_nodes_from(range(E.shape[1]))
G.add_edges_from(links)

#manhattan
k=5
dist_metric="manhattan"
nbrs_man = NearestNeighbors(n_neighbors=k, metric=dist_metric).fit(X)

knn = nbrs_man.kneighbors(return_distance=False)
links = set([])
for i in range(knn.shape[0]):
	for j in knn[i,:]:
		links.add(tuple(sorted((i,j))))

knn_graph=knn.copy()
links=list(links)
G_m = nx.Graph()
G_m.add_nodes_from(range(E.shape[1]))
G_m.add_edges_from(links)

#default non zero
k=5
dist_metric="euclidean"
nbrs = NearestNeighbors(n_neighbors=k, metric=dist_metric).fit(X2)

knn = nbrs.kneighbors(return_distance=False)
links = set([])
for i in range(knn.shape[0]):
	for j in knn[i,:]:
		links.add(tuple(sorted((i,j))))

knn_graph=knn.copy()
links=list(links)
G_non_zero = nx.Graph()
G_non_zero.add_nodes_from(range(E.shape[1]))
G_non_zero.add_edges_from(links)

###forceatlas2 setting###
forceatlas2 = ForceAtlas2(
                          # Behavior alternatives
                          outboundAttractionDistribution=False,  # Dissuade hubs
                          linLogMode=False,  # NOT IMPLEMENTED
                          adjustSizes=False,  # Prevent overlap (NOT IMPLEMENTED)
                          edgeWeightInfluence=1.0,

                          # Performance
                          jitterTolerance=1.0,  # Tolerance
                          barnesHutOptimize=True,
                          barnesHutTheta=2,
                          multiThreaded=False,  # NOT IMPLEMENTED

                          # Tuning
                          scalingRatio=1.0,
                          strongGravityMode=False,
                          gravity=0.05,

                          # Log
                          verbose=True)
        
###set random seed to get the same results###
#default
num_fa2_iter=5000
np.random.seed(0)
pos=np.random.random([E.shape[1],2])        
poslist = np.asarray([pos[i] for i in G.nodes()])
M = nx.to_scipy_sparse_matrix(G, dtype='f', format='lil')
l = forceatlas2.forceatlas2(M, pos=poslist, iterations=num_fa2_iter)
positions = dict(zip(G.nodes(), l))
positions = np.array([positions[i] for i in sorted(positions.keys())])
out_positions=positions.copy()
#output will be almost same as scRepliseq_SPRING_pos.txt
np.savetxt('data/scRepliseq_884cells_Sphase_G1_data_set_filtered_PCA_knn5_pos.txt', out_positions)



