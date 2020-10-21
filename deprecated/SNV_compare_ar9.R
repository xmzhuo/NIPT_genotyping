rm(list=ls())
#Xinming Zhuo, xzhuo@bcm.edu

#this script take vcf files
#The VCF have cutoff of min 50 read depth, min 10 read each strand, 0.05 to call homozygous
#run install gplots at the first time
#install.packages("gplots")
#allow intra case comparion to mother
#need all case initiated with "NIPT" to trim text
# report True as different to mother as in first number, indel as indel in second number
#add amplicon chr info to the table
#add HID-Y as special at the end of table
#make allele between 0.05 to 0.2 as lowercase
#double check the min coverage of minor allele
#incorporate coverage data and alternative allele annotations
#fix problem of list unwanted het mom in short list
#fix problem of extreme low het_call cases
#adding graphic QC
#color to display QC, 1-Black (0-4), 2-Blue(5-9), 3-Green(10-14),4-Green(15-19),5-Orange(20-24),6-Orange(25-29),7-Red(30+)
#shorten the text in scatter plot
#fix the bugs when sample has only one amplicon


min_coverage=200 #the minmum cut off for counting SNV
min_amp_QC=50 #coverage threshold to do the QC with counting amplicon number
het_call=c(0.05,0.20) #cutoff of calling het >0.20 and hom <0.05, anything in between will be in lowercase 
het_min_cov=20 #min coverage for calling allele
truediff_num=2 #the minimum cut off to call a fetus condition 1
truediff_per=0.05 #the minimum percentage to call fetus condition 2

# do not change scrpts below this line #
##############################################################################################
library(gplots) 
library(tcltk)

ref_alt = tk_choose.files(default = "", caption = "Select One referance alternative allele info",
                           multi = FALSE, index = 1)

#vals <- varEntryDialog(vars=c('Variable1', 'Variable2'))
#tk_choose.files(default = "", caption = "Select files",multi = TRUE, filters = NULL, index = 1)
filter_set=matrix(c(".call.csv"),1,2,byrow=TRUE)

arg_mother=tk_choose.files(default = "", caption = "Select One Reference 1 (Mother) cluster file (REQUIRED)",
                           multi = FALSE, filters = filter_set, index = 1)

setwd(dirname(arg_mother))


arg_father=tk_choose.files(default = "can be empty", caption = "Select One Reference 2 (Father) cluster file (OPTIONAL)",
                           multi = FALSE, filters = filter_set, index = 1)

#arg_father=NULL

args_quest=tk_choose.files(default = "", caption = "Select Target cluster files to compare (REQUIRED, at lease one)",
                     multi = TRUE, filters = filter_set, index = 1)



args=c(arg_father,arg_mother,args_quest)



#write your target file name here
#windows
#args=c("C:/Users/xzhuo/Documents/bcm/Beaudet lab/NIPT/NGS/11092016_qun/ampliseq/RE__plugin_results_finished/TSVC_variants_IonXpress_001_743M_G1000.vcf")
#       "C:/Users/xzhuo/Documents/bcm/Beaudet lab/NIPT/NGS/08302016 batch/IonXpress_031_R_2016_08_29_NIPT590M_gDNA.fastq.sam_chrYq11cluster.csv",
#       "C:/Users/xzhuo/Documents/bcm/Beaudet lab/NIPT/NGS/08302016 batch/IonXpress_095_R_2016_08_29_NIPT591M_gDNA.fastq.sam_chrYq11cluster.csv",
#       "C:/Users/xzhuo/Documents/bcm/Beaudet lab/NIPT/NGS/08302016 batch/IonXpress_096_R_2016_08_29_NIPT601M_gDNA.fastq.sam_chrYq11cluster.csv")
#linux
#args=c("/home/xzhuo/Documents/08092016/IonXpress_096_R_2016_08_09_Jlup_0.5_HLA.fastq.sam_chr15q24cluster.csv"
#       ,"/home/xzhuo/Documents/08092016/IonXpress_095_R_2016_08_09_Jlup_1_HLA.fastq.sam_chr15q24cluster.csv")
##############################################################################################
SNP_refalt=read.table(ref_alt,sep=",",head=TRUE,fill=T)

#creat a table for QC report
daf_QC_sum=matrix(NA,length(args),6)
colnames(daf_QC_sum)=c("Sample_ID","Amp_num","QC_score","Dif_SNP","Sus_SNP","Y_SNP")

for (i in 1:length(args)){
  
  address<- args[i]
  da=read.table(address[1],sep=",",head=TRUE,fill=T)
  da=da[da[,15]!="",]
  
  da_QC=da[da[,3]>=min_amp_QC,]
  #count total SNP pass threshold
  if(is.null(nrow(da_QC))){
    
    if(length(da_QC)>0){
      daf_QC_SNP_num = 1
    }else{
      daf_QC_SNP_num = 0
    }
  }else{
    daf_QC_SNP_num = nrow(da_QC)
  }
  
  #count total amplicons
  daf_QC_amp_num=sum(table(da_QC$chr)>0)+
    sum(da_QC$chr %in% "HID")-sum(sum(da_QC$chr %in% "HID")>0)+
    sum(da_QC$chr %in% "HID-Y")-sum(sum(da_QC$chr %in% "HID-Y")>0)
  
  daf_QC_sum[i,2]=daf_QC_amp_num #amplicon number
  daf_QC_sum[i,1]=gsub("\\..*","",gsub(".*NIPT","",basename(args[i]))) #sample name
  daf_QC_sum[i,3]=floor(daf_QC_amp_num/5) #give QC score
  daf_QC_sum[i,6]=sum(da_QC$chr %in% "HID-Y") #Y SNP number
  
  ##start counting for min_coverage in comparison
  da=da[da[,3]>=min_coverage,]
  
  #extract the fraction info
  da_frac=t(apply(da[,5:10], 1, sort))
  
  sort_col<-function(x){
    #gsub("\\_.*","",colnames(sort(x)))
    gsub("\\_|[a-z]","",colnames(sort(x)))
  }
  
  da_call=da[,5:10]
  #colnames(da_call)=c("A","C","G","T","I","D")
  #reassign the het call according to new het call value
  if(nrow(da_call)>0){
    
    for (k in 1:nrow(da_call)){
      da_call_temp=sort_col(da_call[k,])
      
      if(k==1){
        da_calla=da_call_temp
      }else{
        da_calla=rbind(da_calla,da_call_temp)
      }
    }
    
    if(nrow(da_call)>1){
      da_call1=da_calla[,ncol(da_calla)]
      da_call2=da_calla[,ncol(da_calla)-1]
      da_call2[da_frac[,ncol(da_frac)-1]<het_call[1]]=da_call1[da_frac[,ncol(da_frac)-1]<het_call[1]]
      
    }else{
      da_call1=da_calla[length(da_calla)]
      da_call2=da_calla[length(da_calla)-1]
      da_call2[da_frac[ncol(da_frac)-1]<het_call[1]]=da_call1[da_frac[,ncol(da_frac)-1]<het_call[1]]
      
    }
    
    da_callT=paste(da_call1,da_call2,sep="")
    
  }
  
  da_GT=da[,c(13,14,15)]
  da_GT[,1]=da[,1]
  da_GT_cov=da[,3]
  
  #call the genotype
  da_GT[,2]=nchar(as.character(da[,13]))
  
  if (nrow(da_GT)>0){
    #define 0/0 GT
    da_GT[da[,14]=="=",2]=0
    #define 1/1 GT
    da_GT[da_GT[,2]==1,2]=-1
    #define 0/1 GT
    da_GT[da_GT[,2]>=2,2]=1
    #define 1/1 GT
    da_GT[da_GT[,2]==-1,2]=2
    
    #take the first two allele
    #da_GT_3=da_GT[,3]
    da_GT_3=da_callT #reassign the new call first two allele according to minimum het call value
    da_GT_3_1=substr(da_GT_3,1,1)
    da_GT_3_2=substr(da_GT_3,2,2)
    #change the second allele to lowercase if <0.2
    da_GT_3_2[da_frac[,5]>=het_call[1] & da_frac[,5]<het_call[2]]=tolower(da_GT_3_2[da_frac[,5]>=het_call[1] & da_frac[,5]<het_call[2]])
    #ignore minor allele if <10 read, consider as homozygous
    da_GT_3_2[da_frac[,5]*da_GT_cov<het_min_cov]=da_GT_3_1[da_frac[,5]*da_GT_cov<het_min_cov]
    
    da_GT[,3]=paste(da_GT_3_1, da_GT_3_2, sep="")
    
    #da_GT[,3]=substr(da_GT[,3],1,2)
    
    #define function strSort to sort string
    strSort <- function(x)
      sapply(lapply(strsplit(x, NULL), sort), paste, collapse="")
    
    #sort the called genotype to unify the order
    da_GT[,4]= strSort(da_GT[,3])
    
    #da_GT_chr=matrix(rapply(strsplit(da_GT[,3],""),sort),nrow(da_GT),byrow=TRUE)
    
    #define the identity of sample
    da_GT[,5]=i
    #return the GT_3
    #da_GT[,3]=da_GT_3
    da_GT[,3]=da[,11]
    
    #assign coverage data for each sample
    da_GT[,2]=da[,3]
    
    #pileup data
    if(i==1){
      da_GT_all=da_GT
    } else{
      da_GT_all=rbind(da_GT_all,da_GT)
    }
  
  }
  
}

length_args=length(args)
address_all=args

if (length(arg_father)==0){
  #reserve color 1 for father
  da_GT_all[,5]=as.numeric(da_GT_all[,5])+1
  length_args=length_args+1
  #reserve name 1 for father
  address_all=c("",args)
  
  daf_QC_sum=rbind(0,daf_QC_sum)
  daf_QC_sum[1,1]="NA"
}


#consolidate all  
#get all the sample names in short , detect NIPT in front of every sample
sample_names=gsub("\\..*","",gsub(".*NIPT","",basename(address_all)))

#Outer join: merge(x = df1, y = df2, by = "CustomerId", all = TRUE)
#Left outer: merge(x = df1, y = df2, by = "CustomerId", all.x = TRUE)
#Right outer: merge(x = df1, y = df2, by = "CustomerId", all.y = TRUE)
#Cross join: merge(x = df1, y = df2, by = NULL)
dap=da_GT_all[da_GT_all[,5]==1,]
dam=da_GT_all[da_GT_all[,5]==2,]
#colname_temp=c("GT","chr","call","ID")
colnames(dap)=c("rs",paste(sample_names[1],"Coverage",sep=" - "),"chr",paste(sample_names[1],"call",sep=" - "),paste(sample_names[1],"ID",sep=" - "))
colnames(dam)=c("rs",paste(sample_names[2],"Coverage",sep=" - "),"chr",paste(sample_names[2],"call",sep=" - "),paste(sample_names[2],"ID",sep=" - "))

for(j in 3:length_args){
  address<- args[j]
  daf_temp=da_GT_all[da_GT_all[,5]==j,]
  colnames(daf_temp)=c("rs",paste(sample_names[j],"Coverage",sep=" - "),"chr",paste(sample_names[j],"call",sep=" - "),paste(sample_names[j],"ID",sep=" - "))
  
  if (j==3) {
    daf=daf_temp
  } else{
    #merged.data.frame one by one
    daf=merge(daf,daf_temp,by=c("chr","rs"), all=TRUE)
  }
  
}

#merge with mother keep only fetal rs, ID 2
#da_allmerge=merge(dam,daf,by="rs", all.y=TRUE)

#merge with mother keep only fetal and mother shared rs, ID 2
da_allmerge=merge(dam,daf,by=c("chr","rs"), all=FALSE)

#merge with father keep only fetal and mother rs, ID 1
if (length(arg_father)==0){
  dap=dam
  dap[,c(2,4)]=NA
  dap[,5]=1
}
da_allmerge=merge(dap,da_allmerge,by=c("chr","rs"), all.y=TRUE)




#########################################################################################
#need to tell difference 
da_all_m=da_allmerge[,7]
da_all_m=toupper(da_all_m)
#da_all_m=matrix(unlist(strsplit(da_allmerge[,8],"")),nrow(da_allmerge),byrow=TRUE)
#da_all_p=matrix(unlist(strsplit(da_allmerge[,4],"")),nrow(da_allmerge),byrow=TRUE)

#creat vector for filling stat colname
da_stat=matrix(0,length_args,byrow=TRUE)

for (j in c(1,2:length_args)){
  
  da_all_f_temp=da_allmerge[,j*3+1]
  da_all_f_temp=toupper(da_all_f_temp)
  da_all_f_temp[is.na(da_all_f_temp)]="NN"
  da_all_f=matrix(unlist(strsplit(da_all_f_temp,"")),nrow(da_allmerge),byrow=TRUE)
  
  for (k in 1:nrow(da_all_f)) {
    #da_all_f[k,1]=grepl(da_all_f[k,1],da_all_m[k])
    #da_all_f[k,2]=grepl(da_all_f[k,2],da_all_m[k])
    da_all_f_temp[k]=!(grepl(da_all_f[k,1],da_all_m[k]) & grepl(da_all_f[k,2],da_all_m[k]))
  }
  #da_all_f[da_all_f==TRUE]=0
  #da_all_f[da_all_f==FALSE]=1
  da_all_f_temp[is.na(da_allmerge[,j*3+1])]=NA
  #count all negative include indel
  da_all_f_false=table(da_all_f_temp)["FALSE"]
  #mask indel present
  da_all_f_temp[(da_all_f[,1]=="I"|da_all_f[,1]=="D")|(da_all_f[,2]=="I"|da_all_f[,2]=="D")]="Indel"
  da_allmerge[,j*3+2]=da_all_f_temp
  
  da_all_f_true=table(da_all_f_temp)["TRUE"] #sum(da_all_f_temp, na.rm=TRUE)
  #da_all_f_false=table(da_all_f_temp)["FALSE"]
  da_all_f_all=sum(!is.na(da_all_f_temp))
  
  da_all_f_judge = if(da_all_f_true>=truediff_num & da_all_f_true/da_all_f_all>truediff_per & !is.na(da_all_f_true)) {"Fetus"} else {"Uninfo"}
  
  da_stat[j]=paste(da_all_f_judge," : ",
                   da_all_f_true,"/",da_all_f_all," - ",
                   da_all_f_all-da_all_f_false,"/",da_all_f_all,sep="")
  colnames(da_allmerge)[j*3+2]=da_stat[j]
  
  daf_QC_sum[j,4]=if(is.na(da_all_f_true)){0}else{da_all_f_true}
  daf_QC_sum[j,5]=if(is.na(da_all_f_all-da_all_f_false)){0}else{da_all_f_all-da_all_f_false}
  
}

#add ChrY back
ChrY="HID-Y"
#ChrY="HLA-B"
daf_Y=daf[daf[,1]==ChrY, ]

if(nrow(daf_Y)>0) {
  dam_Y=dam[dam[,3]==ChrY, ]
  da_Y_merge_temp=merge(dam_Y,daf_Y,by=c("chr","rs"), all.y=TRUE)
  
  dap_Y=dap[dap[,3]==ChrY, ]
  da_Y_merge=merge(dap_Y,da_Y_merge_temp,by=c("chr","rs"), all.y=TRUE)
  
  da_Y_merge[,3*(1:length_args)+2]=NA
  
  colnames(da_Y_merge)=colnames(da_allmerge)
  da_allmerge=rbind(da_allmerge,da_Y_merge)
}

#incoporate ref snp info
da_allsum = merge(SNP_refalt[,2:4],da_allmerge,by.x="SNP",by.y="rs", all.y=TRUE)

da_allsum=cbind(da_allsum[,4],da_allsum[,-4])

colnames(da_allsum)[1]="Amplicon"
da_allsum=da_allsum[order(da_allsum$Amplicon),]
#########################
#da_allsum = da_allmerge[,-(4*(1:length_args)-1)]
#da_allsum = da_allmerge

da_allsum[is.na(da_allsum)]=""
#da_allsum[da_allsum==0]="0|0"
#da_allsum[da_allsum==1]="0|1"
#da_allsum[da_allsum==2]="1|1"

write.csv(da_allsum,file=paste(gsub("\\..*","",arg_mother),"_case.csv",sep=""),row.names=TRUE)

#only keep homozygous mother or empty
#da_allsum_s=da_allsum[da_allsum[,6]!="0|1"|da_allsum[,6]=="",]

da_allsum_s=da_allsum[substr(da_allsum[,9],1,1)==substr(da_allsum[,9],2,2)|da_allsum[,9]=="",]

#only keep father differ from mother
da_allsum_s=da_allsum_s[da_allsum_s[,7]==TRUE | da_allsum_s[,7]=="",]
#da_allsum_s=da_allsum_s[,-(3*(1:length_args))]
write.csv(da_allsum_s,file=paste(gsub("\\..*","",arg_mother),"_S_case.csv",sep=""),row.names=TRUE)

#add Y SNP number to QC_sum report
daf_QC_sum[j,4]=as.numeric(daf_QC_sum[j,4])+as.numeric(daf_QC_sum[j,6])
daf_QC_sum[j,5]=as.numeric(daf_QC_sum[j,5])+as.numeric(daf_QC_sum[j,6])
write.csv(daf_QC_sum,file=paste(gsub("\\..*","",arg_mother),"_cases_QC.csv",sep=""),row.names=TRUE)

#color to display QC, 1-Black (0-4), 2-Blue(5-9), 3-Green(10-14),4-Green(15-19),5-Orange(20-24),6-Orange(25-29),7-Red(30+)
colorcode=as.numeric(daf_QC_sum[,3])+1
colorcode[colorcode>7]=7

jpeg(paste(gsub("\\..*","",arg_mother),"_QC.jpg",sep=""), width = 6, height = 6, units = "in", res = 600)
par(oma = c(5, 1, 1, 1))
par(mfrow=c(2,1))
par(mar=c(0,1,1,1))
plot(daf_QC_sum[,4], col=c("Black","Violet","Blue","Green","Yellow","Orange","Red")[colorcode],
     ylim=c(10,40),ylab="Dif_SNP",xlab="",xaxt='n',yaxs='i',
     main=paste(gsub("\\_S.*","",daf_QC_sum[2,1]),"_QC",sep=""))

par(mar=c(1.5,1,0,1))
plot(daf_QC_sum[,4], col=c("Black","Blue","Green","Green","Orange","Orange","Red")[colorcode],
     ylim=c(-0.5,10),ylab="Dif_SNP",xlab="",xaxt='n',yaxs='i')
axis(1,at=seq_along(daf_QC_sum[,4]),labels=sub("_S.*_R","_R",(gsub("\\_00.*","",daf_QC_sum[,1]))),las=2,cex.axis=0.5)
#threshold for different
abline(h=2,col="Green")
abline(h=3,col="Orange")
abline(h=4,col="Red")
dev.off()

jpeg(paste(gsub("\\..*","",arg_mother),"_QC2.jpg",sep=""), width = 6, height = 6, units = "in", res = 600)
par(oma = c(5, 1, 1, 1))
plot(daf_QC_sum[,2], daf_QC_sum[,4], ylab="Dif_SNP",xlab="Amp_Num",main=paste(gsub("\\_S.*","",daf_QC_sum[2,1]),"_QC",sep=""),
     xlim=c(0,40),ylim=c(0,30))
text(as.numeric(daf_QC_sum[,2]), as.numeric(daf_QC_sum[,4]),
     labels=sub(".*M-G","G",sub("_S.*_R","_R",(gsub("\\_00.*","",daf_QC_sum[,1])))),
     las=2,cex=0.5,pos=3)
abline(a=0,b=0.1,col="grey")
abline(a=1,b=0.1,col="wheat")
dev.off()
