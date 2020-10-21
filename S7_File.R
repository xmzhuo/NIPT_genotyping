rm(list=ls())
#Xinming Zhuo, xmzhuo@gmail.com
# R_script_to_compare_haplotype_of_NIPT_cells_with_parents
#dependency "tcltk", "ape"
# this script take samcluster.csv files from S3_file.sh for comparison in a trio or duo samples
# it take one ore more files to plot
# the plot and new csv file will saved to the the same location and name as last file
#only for pair wise comparison


#set how many top read to compare default is 10
top_num=10 #number to reads to display
threshold=0.02 #filter for noise level 0-1
cut_height=2 #threshold for cutting the distance between two sample, the higher the more strict

# do not change scrpts below this line #
##############################################################################################

library(tcltk)

filter_set=matrix(c("CSV",".csv"),1,2,byrow=TRUE)


args_quest=tk_choose.files(default = "", caption = "Select Target cluster files to compare (REQUIRED, at lease one)",
                     multi = TRUE, filters = filter_set, index = 1)

setwd(dirname(args_quest[1]))

args=c(args_quest)

##############################################################################################

arg_num=length(args)

cut_num=2+arg_num

address<- args[1]
da=read.table(address[1],skip = 1,sep=",",fill=T)
da[,5]=1
min_num=min(top_num, nrow(da))
db=da[1:min_num,]

#combine reads in each cluster in each sample
da_sum=matrix(NA,max(da[,3]),5)
for (sum_num in 1:max(da[,3])){
  
  da_sum1=da[da[,3]==sum_num,]
  da_sum[sum_num,1]=as.character(da_sum1[da_sum1[,2]==max(da_sum1[,2]),1])[1]
  da_sum[sum_num,2]=sum(da_sum1[,2])
  da_sum[sum_num,4]=sum(da_sum1[,4])
  da_sum[sum_num,3]=sum_num
  da_sum[sum_num,5]=da[1,5]
  
}
db_sum=da_sum

for (i in 2:arg_num){
  address[i]<-args[i]
  da=read.table(address[i],skip = 1,sep=",",fill=T)
  da[,5]=i
  min_num=min(top_num, nrow(da))
  db=rbind(db,da[1:min_num,])
  
  da_sum=matrix(NA,max(da[,3]),5)
  for (sum_num in 1:max(da[,3])){
    
    da_sum1=da[da[,3]==sum_num,]
    da_sum[sum_num,1]=as.character(da_sum1[da_sum1[,2]==max(da_sum1[,2]),1])[1]
    da_sum[sum_num,2]=sum(da_sum1[,2])
    da_sum[sum_num,4]=sum(da_sum1[,4])
    da_sum[sum_num,3]=sum_num
    da_sum[sum_num,5]=da[1,5]
    
  }
  
  db_sum=rbind(db_sum,da_sum)
  
  
}


#reset cut_height according to SNP number
SNP_num=nchar(as.character(db[1,1]))
if (SNP_num<20) {cut_height=cut_height^(SNP_num/20)}#1+cut_height^(SNP_num/20)/2}

cluster_name=gsub("\\_cluster.*","",gsub(".*.sam_","",basename(args[1])))

##################make cluster 
db=db[db[,4]>=threshold,]
db_a_1=gsub("[*]","",db[,1])
d<-adist(db_a_1,costs=list(ins=0.5,del=0.5,sub=1))
#d<-adist(db[,1])
rownames(d)=db[,1]
hc<-hclust(as.dist(d))

db[,6]=cutree(hc,h=cut_height+SNP_num/30)
df1=rbind(matrix(NA,arg_num,ncol(db)),db)
df1[1:arg_num,5]=1:arg_num
df1[1:arg_num,1]=address
dff<-data.frame(df1)
#creat new folder
address_folder=paste(dirname(args[1]),"/pair",sep="")
dir.create(address_folder,showWarnings =F)
write.csv(dff,file=paste(address_folder,"_compareall.csv",sep=""),row.names=FALSE)

db_sum=db[db_sum[,4]>=threshold,]
db_sum=db_sum[complete.cases(db_sum),]
db_sum_a_1=gsub("[*]","",db_sum[,1])
d_sum<-adist(db_sum_a_1,costs=list(ins=0.5,del=0.5,sub=1))
#d_sum<-adist(db_sum[,1])
rownames(d_sum)=db_sum[,1]
hc_sum<-hclust(as.dist(d_sum))

db_sum[,6]=cutree(hc_sum,h=cut_height+SNP_num/30)
df1_sum=rbind(matrix(NA,arg_num,ncol(db_sum)),db_sum)
df1_sum[1:arg_num,5]=1:arg_num
df1_sum[1:arg_num,1]=address
dff_sum<-data.frame(df1_sum)

write.csv(dff_sum,file=paste(address_folder,"_compareall_sum.csv",sep=""),row.names=FALSE)


library(ape)
#making cluster plot with bar plot

###########make png ###############
png(filename = paste(address_folder,"_compareall.png",sep=""), width = 6, height = 6+32*arg_num/72, units = "in", res = 300)

par(oma = c(arg_num, 1, 0.5, 1))

layout(matrix(c(1,1,2,2,2,2,2,2,3,3,3,3,3,3), nrow = 7, ncol = 2, byrow = TRUE))
#making bar plot according to cluster
par(mar = c(0, 1, 2, 1))

barplot(db[hc$order,2],las=3,cex.names=1,col=db[hc$order,5],border=NA,log="y",ylim=c(10,1000))
abline(h=c(50,100,500,1000),col='wheat')
par(mar = c(0.5, 1.5, 0.5, 1.5))

ph <- as.phylo(hc)  
plot(ph,tip.color=db[,5],direction = "downwards",cex=0.6, font=2,adj=0,
     underscore=TRUE,label.offset=0.1,family="mono",las=3) #hang = -1)
axisPhylo(2, las = 1)

rect.hclust(hc,h=cut_height+SNP_num/30,border=7)

#making bar plot according to cluster
par(mar = c(0, 1, 0.2, 1))
barplot(db[hc$order,4],las=3,cex.names=1,col=db[hc$order,5],border=NA)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", basename(address), xpd = TRUE, horiz = FALSE, 
       inset = c(0,0), bty = "n", pch=15,col = min(db[,5]):max(db[,5]), 
       cex = 0.75)

dev.off()

png(filename = paste(address_folder,"_compareall_sum.png",sep=""), width = 6, height = 6+32*arg_num/72, units = "in", res = 300)

par(oma = c(arg_num, 1, 0.5, 1))

layout(matrix(c(1,1,2,2,2,2,2,2,3,3,3,3,3,3), nrow = 7, ncol = 2, byrow = TRUE))
#making bar plot according to cluster
par(mar = c(0, 1, 2, 1))

barplot(as.numeric(db_sum[hc_sum$order,2]),las=3,cex.names=1,col=db_sum[hc_sum$order,5],border=NA,log="y",ylim=c(10,1000))
abline(h=c(50,100,500,1000),col='wheat')
par(mar = c(0.5, 1.5, 0.5, 1.5))

ph_sum <- as.phylo(hc_sum)  
plot(ph_sum,tip.color=db_sum[,5],direction = "downwards",cex=0.6, font=2,adj=0,
     underscore=TRUE,label.offset=0.1,family="mono",las=3) #hang = -1)
axisPhylo(2, las = 1)

rect.hclust(hc_sum,h=cut_height+SNP_num/30,border=7)#h=2*cut_height,border=7)#k=cut_num,border=7,which=c(1:cut_num,1))

#making bar plot according to cluster
par(mar = c(0, 1, 0.2, 1))
barplot(as.numeric(db_sum[hc_sum$order,4]),las=3,cex.names=1,col=db_sum[hc_sum$order,5],border=NA)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", basename(address), xpd = TRUE, horiz = FALSE, 
       inset = c(0,0), bty = "n", pch=15,col = min(db_sum[,5]):max(db_sum[,5]), 
       cex = 0.75)

dev.off()



#################################making cluster file one by one

arg_num_1=length(args)

ROC_table=matrix(NA,arg_num*(arg_num-1)/2,5)
table_row=0

#make a 2 way table
ROC_table2=matrix(NA,arg_num_1+2,arg_num_1+2)
ROC_table2[1,3:(arg_num_1+2)]=basename(args)
ROC_table2[2,3:(arg_num_1+2)]=1:arg_num_1
ROC_table2[3:(arg_num_1+2),1]=basename(args)
ROC_table2[3:(arg_num_1+2),2]=1:arg_num_1



for (quest_num in 1:(arg_num_1-1)) {
  
  for (quest_num2 in (quest_num+1):arg_num_1) {
    
    #cluster 1 by 1
    db_1=db[is.element(db[,5],c(quest_num,quest_num2)),]
    
    db_1_1=gsub("[*]","",db_1[,1])
    d_1<-adist(db_1_1,costs=list(ins=0.5,del=0.5,sub=1))
    
    rownames(d_1)=db_1[,1]
    hc_1<-hclust(as.dist(d_1))
    
    db_1[,6]=cutree(hc_1,h=cut_height+SNP_num/30)#k=5)
    df1_1=rbind(matrix(NA,3,ncol(db_1)),db_1)
    
    df1_1[1,5]="2"
    df1_1[1,4]="2"
    df1_1[1,1]=args[quest_num2]
    
    df1_1[2,5]="1"
    df1_1[2,4]="1"
    df1_1[2,1]=args[quest_num]
     
    dff_1<-data.frame(df1_1)
    
    write.csv(dff_1,file=paste(address_folder,"/",quest_num,"_",quest_num2,".pair.csv",sep=""),row.names=FALSE)
    
    table_row=table_row+1
    
    ROC_table[table_row,1:2]=gsub("\\.fastq.*","",basename(args[c(quest_num,quest_num2)]))
    #color label
    db_temp=db_1
    db_temp[!is.element(db_1[,5],c(quest_num,quest_num2)),6]="8"
    db_temp[db_1[,5]==quest_num,6]="1"
    db_temp[db_1[,5]==quest_num2,6]="2"
    
    #gnerate cluster types for each
    db_temp1=sort(unique(db_1[db_temp[,6]==1 & db_temp[,4]>=5*threshold,6]))
    db_temp2=sort(unique(db_1[db_temp[,6]==2 & db_temp[,4]>=5*threshold,6]))
    db_tempc=if(identical(db_temp1,db_temp2)) {0} else {1}
    ROC_table[table_row,3]=toString(db_temp1)
    ROC_table[table_row,4]=toString(db_temp2)
    ROC_table[table_row,5]=db_tempc
    
    h=cut_height+SNP_num/30
    colnames(ROC_table)=c("Sample_A","Sample_B",
                          paste(cluster_name,h,"A_Hap",sep="-"),
                          paste(cluster_name,h,"B_Hap",sep="-"),
                          paste(cluster_name,h,sep="-"))
    ROC_table2[quest_num2+2,quest_num+2]=db_tempc  
    
  }
 
}

write.csv(ROC_table,file=paste(address_folder,"/",cluster_name,"_ROC_pair.csv",sep=""),row.names=FALSE)

write.csv(ROC_table2,file=paste(address_folder,"/",cluster_name,"_ROC_table.csv",sep=""),row.names=FALSE)

