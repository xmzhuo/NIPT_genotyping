options(warn =-1)
# this program read sam files 
# it also subset the read sequence according to bed file
# Makeing a cluster for the most populus reads with similarity
# vectorized to speed up
#length limit float with amplicon lenght
#make folder for each amplicon
#adjust super short reads min_length filter  line 68-73

#adjust for illumina 150 se or pe
#1c1pe modified to fit pe on both end, introduce adjustable base quality and map quality filtering 
#missing as *, deletion as _, padding as N, low score base as U

rm(list=ls())

#options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

address <- args[1]
SNPsite <- args[2]

read_num <- args[3]
read_len <- args[4]

read_num = as.numeric(as.character(read_num))
read_len = as.numeric(as.character(read_len)) #setting of read length of illumina

print(paste(read_num, "*",read_len,"bp", sep=""))

#to be hided
#address="~/Documents/summary/test_pe/30423-Gdna_S17_L001_R_001.fastq.gz.c.sam"
#SNPsite="~/Documents/multiplex_v1h.fa.SNP.bed"


min_length_m=read_len+10 #min illumina read + 10 padding N
sample_size=5000 #
primter_length=25 #adjust display number of SNP 
min_read=10
min_type=10
MAPQ_min=30
BaseQ_min=13 #95% time is correct

# load SNP sites
SNPSiteAll=read.table(SNPsite,skip=1,sep="")
#count how many amplicons in SNP file
SNP_chr=table(SNPSiteAll[,1])

#read sam files
t=read.table(address,skip = nrow(SNP_chr),sep="",fill=T,row.names=NULL)

#remove empty rowst
t1=t[complete.cases(t),]
#remove unmapped reads
t1=t1[t1[,3]!="*",]
t1[,14:16]=""
#calculate read length
t1[,17]=nchar(as.character(t1[,10]))

#filter read length
t1=t1[t1[,17]>=min_length_m,]
#filter low MAPQ
t1=t1[t1[,5]>=MAPQ_min,]

t_all=t1
rm(t)


for (chr_num in 1:nrow(SNP_chr)){
  print(c(chr_num, SNP_chr[chr_num]))
  
  #one amplicon at a time, extract reads and SNPsites
  t1=t_all[t_all[,3]==rownames(SNP_chr)[chr_num],]
  
  SNPSite=SNPSiteAll[SNPSiteAll[,1]==rownames(SNP_chr)[chr_num],]
  
  SNPSite_start=SNPSite[,2]
  SNPSite_stop=SNPSite[,3]
  #get the min length according to each amplicon
  min_length=min(SNPSite_stop[length(SNPSite_stop)]-SNPSite_start[1])
  min_length_SNP=min_length #set min_length to display SNP
  #set the min_length_SNP at least 150-25bp
  if (min_length_SNP<(read_num*read_len)-primter_length){
    min_length_SNP=(read_num*read_len)-primter_length
  }
  
  row_num_t1=nrow(t1)
  print(row_num_t1)
  
  #adjusting the min_length with median length
  if (row_num_t1>0){
    if (median(t1[,17],na.rm=TRUE)<135) {
      min_length=min_length_m
      #min_length_SNP=mean(c(min_length_SNP,min_length))
      min_length_SNP=sqrt(min_length_SNP * min_length)
    }
  }
  
  
  #filter read length according to each chromosome
  #t1=t1[t1[,17]>min_length,] turn off filter for illumina
  
  row_num_t1=nrow(t1)
  print(row_num_t1)
  
  # exclude the case of HID HID-Y and Disease
  if (!(rownames(SNP_chr)[chr_num] %in% c("HID","HID-Y","Disease"))){
    
    #downsize sample size
    if (row_num_t1>sample_size) t1=t1[sample(row_num_t1,size=sample_size,replace=FALSE),]
    
    row_num_t1d=nrow(t1)
    print(row_num_t1d)
    
    #vectorized the dataframe
    t1_chr=t1[,3]
    t1_pos=t1[,4]
    t1_cig=t1[,6]
    t1_seq=t1[,10]
    t1_seq_col=t1[,14]
    t1_seq_ex=t1[,15]
    t1_seq_exnum=t1[,16]
    t1_seq_score=t1[,11]
    
    #minimu 10 reads to wrok with
    if (row_num_t1d>min_read) {
      #start extracting sequence according to current amplicon
      for (i in 1:row_num_t1d){
        #get all the numbers and operations
        cigar_num=gsub("([aA-zZ]+)", ",", as.character(t1_cig[i]))
        cigar_chr=gsub("([0-9]+)", "", as.character(t1_cig[i]))
        cigar_num1=strsplit(cigar_num,",")
        cigar_chr1=strsplit(cigar_chr,"")
        cigar_num1=unlist(cigar_num1)
        cigar_chr1=unlist(cigar_chr1)
        
        #cigar indext start and stop
        cigar_num2=cigar_num1
        cigar_num2[cigar_chr1=="D"]=0
        ci_stop=cumsum(cigar_num2)
        ci_start=ci_stop-as.numeric(cigar_num2)+1
        
        start_pos=ci_start[cigar_chr1=="M"]
        stop_pos=ci_stop[cigar_chr1=="M"]
        
        #identify the padding SNP position in ref, include 5bp of buffer zone of each end
        amp_match_len=sum(as.numeric(cigar_num1[cigar_chr1 %in% c("M","D")]))
        Padding_start=t1_pos[i]+read_len-5
        Padding_stop=t1_pos[i]+amp_match_len-1-read_len+5
        #identify the padding SNP in SNP list
        PadSNP_start=which(SNPSite_start >=Padding_start & SNPSite_start <=Padding_stop)
        PadSNP_stop=which(SNPSite_stop >=Padding_start & SNPSite_stop <=Padding_stop)
        
        #get base quality score
        temp_seq_score=as.numeric(charToRaw(as.character(t1_seq_score[i])))-33
        temp_seq=unlist(strsplit(as.character(t1_seq[i]),""))
        temp_seq[temp_seq_score<BaseQ_min] = "U" #replace low quality (<20) base with U
        temp_seq=paste(temp_seq,sep="",collapse="")
        #extract sequence (match sequence only)
        #t2=substr(rep(t1_seq[i], length(start_pos)), start_pos, stop_pos)
        t2=substr(rep(temp_seq, length(start_pos)), start_pos, stop_pos)
        
        t3=matrix("",1,length(cigar_chr1))
        #insert deletion marks
        t3[cigar_chr1=="M"]=t2
        for (j in 1:length(cigar_chr1)){
          if (cigar_chr1[j]=="D") t3[j]=paste(replicate(cigar_num1[j],"_"),collapse="")
        }
        
        
        #t3[cigar_chr1=="D"]="_"
        t4=paste(t3,sep="",collapse="")
        t1_seq_col[i]=t4
        
        #SNP start and end
        start_num=1
        end_num=length(SNPSite_stop)
        #end_num=length(SNPSite_stop[SNPSite_start<(SNPSite_stop[1]+min_length_SNP-primter_length)])
        
        #only show dsignated number of SNP: end_num-start_num+1
        t5=substr(rep(t4, end_num-start_num+1), 
                  (SNPSite_start[start_num:end_num]-as.numeric(as.character(t1_pos[i]))+1), 
                  (SNPSite_stop[start_num:end_num]-as.numeric(as.character(t1_pos[i]))+1))
        #insert * for missing read in certain SNP position
        t6=replace(t5, t5=="", "*")
        
        #replace nt call in padding region with "N"
        t6_temp=unlist(strsplit(as.character(t6),""))
        t6_temp[PadSNP_start]="N"
        t6=paste(t6_temp,sep="",collapse="")
      
        t1_seq_ex[i]=paste(t6,sep="",collapse="")
        # count actual length w/o *  in the extracted sequence
        t1_seq_exnum=nchar(paste(t5,sep="",collapse=""))
        
      }
      
      
      a=sort(table(as.character(t1_seq_ex)),decreasing = TRUE)
      if(!is.null(nrow(a))){
        if (nrow(a)>2) {
          #Levenshtein Distance
          top_num=min(nrow(a),min_type)
          d<-adist(rownames(a[1:top_num]))
          rownames(d)=rownames(a[1:top_num])
          hc<-hclust(as.dist(d))
          
          #dataframe of top  20 reads and grouping
          df<-data.frame(a[1:top_num],cutree(hc,k=max(2,ceiling(top_num/3))))
          
          #percentage of each group
          df[,3]=df[,1]/sum(df[,1])
          
          
          #write.csv(df,file=paste(address,"_",rownames(SNP_chr)[chr_num],"cluster.csv",sep=""),row.names=TRUE)
          #save file into new folder
          address_folder=paste(dirname(address),"/",rownames(SNP_chr)[chr_num],sep="")
          dir.create(address_folder,showWarnings =F)
          write.csv(df,file=paste(address_folder,"/",basename(address),"_",rownames(SNP_chr)[chr_num],"_cluster.csv",sep=""),row.names=TRUE)
          
          #making cluster plot with bar plot
          png(paste(address_folder,"/",basename(address),"_",rownames(SNP_chr)[chr_num],".png",sep=""),
              width = 6, height = 6, units = "in", res = 300)
          
          par(oma = c(1, 1, 1, 1))
          par(mfrow=c(2,1))
          par(mar = c(0, 1.5, 0.5, 1.5))
          plot(hc,cex=0.6,main=paste(basename(address),"_",rownames(SNP_chr)[chr_num],sep=""), 
               cex.main=0.75, axes=T, hang=-1,family = "mono", font = 2)
          rect.hclust(hc,k=max(2,ceiling(top_num/3)))
          #dev.new(width=5, height=4)
          #making bar plot according to cluster
          par(mar = c(0, 1, 0.2, 1))
          barplot(a[hc$order],las=3,cex.names=1)
          par(new = TRUE)
          #add second graph of percentage overlay
          barplot(df[hc$order,3],axes=F,xlab=NA, ylab=NA,las=3,cex.names=1)
          axis(4, las=3)
          
          dev.off()
        }
      }
      
      
    }
    
  } 
  
  
}


