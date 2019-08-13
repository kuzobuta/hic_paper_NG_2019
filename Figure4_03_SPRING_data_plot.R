#Plot of SPRING result

library(msir)
library(fields)
library(matlab)

dir="/Users/hisashimiura/Documents/RIKEN_2019/190722_HiC_paper_custom_scripts/"
setwd(dir)

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

#########Figure 4c & 4d#######
pdf("Figures/Figure4c.pdf")
par(mar=c(5,4,3,6))
plot(spring_out$V1,spring_out$V2,pch=21,
     bg=cols1_2,cex=1.0,
     main="SPRING k=5",
     xlab="X",ylab="Y")
par(xpd=TRUE)
legend(x=par()$usr[2],y=par()$usr[4],legend=sample_sizes$Sample_set,
       pch=20,col=c(rainbow(8)[1:7],"black","grey"),cex=0.7,xjust=-.1)
dev.off()

pdf("Figures/Figure4d.pdf")
col4_srt=jet.colors(101)[ceiling(all_srt_r_scores*100)+1]
par(mar=c(5,4,3,6))
plot(spring_out$V1,spring_out$V2,pch=21,
     col="black",bg=col4_srt,cex=1.0,
     main="SPRING k=5 vs Repliscores",
     xlab="X",ylab="Y")
image.plot(legend.only=TRUE, zlim= range(0,100),col = jet.colors(101)) 
dev.off()

#########Figure 4e#######
#Plot data EpiSC vs each time points#
#save the data of loess#
loess_line_list=list()
for (i in 1:8){
  id_sets=c(pos[i,1]:pos[i,2])
  df=data.frame(a=spring_out$V1[id_sets],
                b=spring_out$V2[id_sets])
  l <- loess(b ~ a,df)
  loess_line_list[[i]] = data.frame(x=df$a[order(df$a)],
                                    y=predict(l,df$a[order(df$a)]))
}

pdf("Figures/Figure4e.pdf")
#white : outlier cells
#grey  : repliscore > 0.05 & < 0.95

color_maps=c(rainbow(8)[1:7],"black")

par(mfrow=c(2,2))
for (i in 1:8){
  plot(spring_out$V1,spring_out$V2,pch=21,
       col="black",bg=cols1_2,cex=1.5,type="n",
       xlab="X",ylab="Y",
       main=sample_sizes[i,1])
  id_sets=c(pos[i,1]:pos[i,2])
  id_sets_EpiSCs=c(pos[8,1]:pos[8,2])
  id_sets_wG1=c(id_sets,pos[9,1]:pos[9,2])
  rep_scores_select=all_srt_r_scores[id_sets]
  points(spring_out$V1[id_sets],
         spring_out$V2[id_sets],pch=21,
         col="black",
         bg=cols1_2[id_sets],
         cex=1.5)
  
  df=data.frame(a=spring_out$V1[id_sets],
                b=spring_out$V2[id_sets])
  l <- loess(b ~ a,df)
  lines(df$a[order(df$a)],predict(l,df$a[order(df$a)]),col=color_maps[i],lwd=2,lty=1)
  med=quantile(df$b - predict(l,df$a),c(0.25,0.75))
  iqr=IQR(df$b - predict(l,df$a))
  print(c(med[1] - 1.5*iqr,med[2] + 1.5*iqr))
  tUp=med[2] + 1.5*iqr;tDwn=med[1] - 1.5*iqr
  df2 = df[which((df$b - predict(l,df$a))>tDwn & 
                   (df$b - predict(l,df$a))<tUp),]
  l2 <- loess(b ~ a,df2)
  #lines(df2$a[order(df2$a)],predict(l2,df2$a[order(df2$a)]),lwd=2,col="black")
  df3 = df[which((df$b - predict(l,df$a))<=tDwn | 
                   (df$b - predict(l,df$a))>=tUp),]
  points(df3$a,
         df3$b,pch=21,
         col="black",
         bg="white",
         cex=1.5)
  points(df$a[rep_scores_select<0.05 | rep_scores_select>0.95],
         df$b[rep_scores_select<0.05 | rep_scores_select>0.95],
         pch=21,
         col="black",
         bg="grey",
         cex=1.5)
  for (j in grep(i,seq(1,8),invert = T)){
    if(j==8){
      lines(loess_line_list[[j]]$x,loess_line_list[[j]]$y,
            lwd=2,col=color_maps[j])
    }else{
      lines(loess_line_list[[j]]$x,loess_line_list[[j]]$y,
            lwd=1,col=color_maps[j],lty=2)  
    }
  }
  total_cells_wo_Grey=length(df$a[rep_scores_select > 0.05 & rep_scores_select < 0.95])
  outliner_cells=sum(((df$b - predict(l,df$a))<=tDwn | (df$b - predict(l,df$a))>=tUp ) &
                       (rep_scores_select > 0.05 & rep_scores_select < 0.95))
  legend("bottomright",
         paste(outliner_cells,"/",total_cells_wo_Grey,
               "(",round(outliner_cells/total_cells_wo_Grey*100,1),")",sep=""),
         pch=21,col="black")
}
dev.off()


