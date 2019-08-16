#Plot of SPRING result

library(msir)
library(fields)
library(matlab)


load("data/scRepliseq_884cells_Sphase_G1_data_set.Rdata") #As scRepliseq_884cells_Sphase_G1_data_set
sample_sizes=read.table("data/scRepliseq_884cells_Sphase_G1_data_set_numbers.txt",sep="\t",header=T)

repscores=function(x){
  score=sum(x==1)/sum(x!=0)
  score
}
all_srt_r_scores=apply(scRepliseq_884cells_Sphase_G1_data_set[,4:dim(scRepliseq_884cells_Sphase_G1_data_set)[2]],2,repscores)

spring_out=read.table("data/scRepliseq_SPRING_pos.txt")

#########colors#######
cols1_2=c(rep(rainbow(8)[1],sample_sizes[1,2]),
          rep(rainbow(8)[2],sample_sizes[2,2]),
          rep(rainbow(8)[3],sample_sizes[3,2]),
          rep(rainbow(8)[4],sample_sizes[4,2]),
          rep(rainbow(8)[5],sample_sizes[5,2]),
          rep(rainbow(8)[6],sample_sizes[6,2]),
          rep(rainbow(8)[7],sample_sizes[7,2]),
          rep("black",sample_sizes[8,2]),
          rep("grey",sample_sizes[9,2]))

s=as.numeric(sample_sizes[,2])
pos=rbind(c(1,s[1]),
          c(1+sum(s[1:1]),sum(s[1:2])),
          c(1+sum(s[1:2]),sum(s[1:3])),
          c(1+sum(s[1:3]),sum(s[1:4])),
          c(1+sum(s[1:4]),sum(s[1:5])),
          c(1+sum(s[1:5]),sum(s[1:6])),
          c(1+sum(s[1:6]),sum(s[1:7])),
          c(1+sum(s[1:7]),sum(s[1:8])),
          c(1+sum(s[1:8]),sum(s[1:9])))

############save the data of loess#############
loess_line_list=list()
for (i in 1:8){
  id_sets=c(pos[i,1]:pos[i,2])
  df=data.frame(a=spring_out$V1[id_sets],
                b=spring_out$V2[id_sets])
  l <- loess(b ~ a,df)
  loess_line_list[[i]] = data.frame(x=df$a[order(df$a)],
                                    y=predict(l,df$a[order(df$a)]))
}

####################Mark the scRT data within IQR & repliscores#########
sample_names=colnames(scRepliseq_884cells_Sphase_G1_data_set)[4:dim(scRepliseq_884cells_Sphase_G1_data_set)[2]]
set_names=c("d0","d2","d3","d4","d5","d6","d7","EpiSC","G1")
all_filter_names=NULL
x=scRepliseq_884cells_Sphase_G1_data_set
for (i in 1:8){
  id_sets=c(pos[i,1]:pos[i,2])
  tmp_data=x[,id_sets+3]
  tmp_names=sample_names[id_sets]
  rep_scores_select=all_srt_r_scores[id_sets]
  df=data.frame(a=spring_out$V1[id_sets],
                b=spring_out$V2[id_sets])
  l <- loess(b ~ a,df)
  med=quantile(df$b - predict(l,df$a),c(0.25,0.75))
  iqr=IQR(df$b - predict(l,df$a))
  #print(c(med[1] - 1.5*iqr,med[2] + 1.5*iqr))
  tUp=med[2] + 1.5*iqr;tDwn=med[1] - 1.5*iqr
  tmp_names_2=cbind(set_names[i],tmp_names,NA,NA,NA,NA) #set, sample name, G1, IQR, Repliscore, Man
  #apply IQR filtering (1st)
  tmp_names_2[which((df$b - predict(l,df$a))<=tDwn | 
                      (df$b - predict(l,df$a))>=tUp),4]=1
  #apply 5, 95% rep score (2nd)
  tmp_names_2[which(rep_scores_select<=0.05 |
                    rep_scores_select>=0.95),5]=1
  #Man filtering after 1st, 2nd filtering
  tmp_data_select=tmp_data[,which((df$b - predict(l,df$a)) > tDwn & 
                                   (df$b - predict(l,df$a)) < tUp  &
                                   rep_scores_select>0.05 &
                                   rep_scores_select<0.95)]
  tmp_data_select_r=rep_scores_select[which((df$b - predict(l,df$a)) > tDwn & 
                                              (df$b - predict(l,df$a)) < tUp  &
                                              rep_scores_select>0.05 &
                                              rep_scores_select<0.95)]
    
  tmp_data_select_srt=tmp_data_select[,order(as.vector(tmp_data_select_r))]
  tmp_data_select_srt[tmp_data_select_srt==0]=NA
  tmp_data_select_srt[tmp_data_select_srt==-1]=0
  tmp_data_select_d=as.matrix(dist(t(tmp_data_select_srt),method="manhattan"))
  n=dim(tmp_data_select_d)[1]
  if(i==1 | i==7){
  ###############################
  med_d=NULL
  for (j in 1:dim(tmp_data_select_d)[1]){
      if (j <= 1+5){
        idx=1:11;idx=idx[idx!=j]
        med_d=c(med_d,median(tmp_data_select_d[j,idx]))
      }else if (j >= (dim(tmp_data_select_d)[1]-5)){
        idx=(dim(tmp_data_select_d)[1]-10):dim(tmp_data_select_d)[1];idx=idx[idx!=j]
        med_d=c(med_d,median(tmp_data_select_d[j,idx]))
      }else{
        idx=(j-5):(j+5);idx=idx[idx!=j]
        med_d=c(med_d,median(tmp_data_select_d[j,idx]))
      }
  }
  ###############################
  filter_IQR=med_d > (quantile(med_d,0.75)+1.5*IQR(med_d))
  out_cells=colnames(tmp_data_select_d)[med_d > (quantile(med_d,0.75)+1.5*IQR(med_d))]
  tmp_data_select_d_f=tmp_data_select_d[!filter_IQR,!filter_IQR]
  
  par(mfrow=c(3,2))
  par(mar=c(3,3,3,3))
  image.plot(t(tmp_data_select_d),main=paste0(set_names[i],"\nManhattan distance","\nn=",n),cex.main=0.8)
  image.plot(t(tmp_data_select_d_f),main=paste0(set_names[i],"\nAfter filtering by Manhattan distance\nn=",dim(tmp_data_select_d_f)[1]),cex.main=0.8)
  barplot(tmp_data_select_r[order(tmp_data_select_r)],
          ylim=c(0,1),border="grey",space=0,xaxt='n',xlab="repliscores",
          main=paste0(set_names[i]," Repliscores\n","n=",n))
  rp=tmp_data_select_r[order(tmp_data_select_r)]
  barplot(rp[!filter_IQR],
          ylim=c(0,1),border="grey",space=0,xaxt='n',xlab="repliscores",
          main=paste0(set_names[i]," Repliscores\n","n=",length(rp[!filter_IQR])))
  x1=tmp_data_select_r[order(tmp_data_select_r)]
  y1=med_d
  plot(x1,y1,pch=20,
       main=paste0(set_names[i]," 3rd quantile + 1.5xIQR filter"),cex.main=0.8,
       xlab="Repliscores",ylab="Manhattan distance")
  points(x1[filter_IQR],y1[filter_IQR],pch=20,col="red")

  man_filters=out_cells
  for (man_filter in man_filters){
    tmp_names_2[grep(man_filter,tmp_names_2[,2]),6]=1
  }
  }
  colnames(tmp_names_2)=c("set","library_id","G1_cells","IQR_filter","out_5_95_rep","Man_filter")
  if(i==1){
    all_filter_names=tmp_names_2
  }else{
    all_filter_names=rbind(all_filter_names,tmp_names_2)
  }
}

i=9
id_sets=c(pos[i,1]:pos[i,2])
tmp_data=x[,id_sets+3]
tmp_names=sample_names[id_sets]
rep_scores_select=all_srt_r_scores[id_sets]
tmp_names_2=cbind(set_names[i],tmp_names,1,NA,NA,NA)
all_filter_names=rbind(all_filter_names,tmp_names_2)

all_filter_scores_names=cbind(names(all_srt_r_scores),as.vector(all_srt_r_scores),all_filter_names)
colnames(all_filter_scores_names)[1:2]=c("lib_id","repscores_woX")
#identical(all_filter_scores_names[,1],all_filter_scores_names[,4])
#[1] TRUE

out_all_filter_scores_names=all_filter_scores_names[,c(3,1,2,5:8)]
write.table(out_all_filter_scores_names,"data/All_Sample_Info_filtering.txt",sep="\t",col.names=T,row.names=F,quote=F)


