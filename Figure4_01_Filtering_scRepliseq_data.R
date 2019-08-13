####Filtering scRepliseq data####

countZeros=function(x){
  z=sum(x==0)
  return(z)
}

load("data/scRepliseq_884cells_Sphase_G1_data_set.Rdata") #As scRepliseq_884cells_Sphase_G1_data_set
CountZero_bins_wG1=apply(scRepliseq_884cells_Sphase_G1_data_set[,4:dim(scRepliseq_884cells_Sphase_G1_data_set)[2]],1,countZeros)
#Filtering genomic positions No data (0) in > 50 samples
scRepliseq_884cells_Sphase_G1_data_set_filter=scRepliseq_884cells_Sphase_G1_data_set[CountZero_bins_wG1<=50,]

#1 copy to 1, 2 copies to 2, no data to 0 
#data without genomic position
g=scRepliseq_884cells_Sphase_G1_data_set_filter[,4:dim(scRepliseq_884cells_Sphase_G1_data_set_filter)[2]]
g2=g
g[g2==-1]=1
g[g2==1]=2

write.table(g,
            "data/scRepliseq_884cells_Sphase_G1_data_set_filtered.txt",
            row.names=F,col.names=T,sep="\t",quote=F)

