options(warn =-1)

# R_script_for_genotype_calling
# this script is used by S3_file.sh 
# bwa-mem samtools bamreadcounts report format as Nextgene
#make different folder for amplicon
#use rs for annotation

rm(list=ls())

library(gdata) 

args <- commandArgs(trailingOnly = TRUE)
print(args)
name <- args[1]
SNPsite <- args[2]

#to increase sensitivity, you can modify min_coveage and min_call from 100 and 0.05 to 20 and 0.01
min_coverage=100    #threshold  for coverage
min_call=0.05        #threshold to call SNP  

#input filename

M_gDNA<-read.csv(name,sep="",FALSE,colClasses=rep("character",200))
ref_position=29911000 #

head1=c("chr","position","reference_base","depth","base","count","avg_mapping_quality","avg_basequality","avg_se_mapping_quality","num_plus_strand","num_minus_strand","avg_pos_as_fraction","avg_num_mismatches_as_fraction","avg_sum_mismatch_qualities","num_q2_containing_reads","avg_distance_to_q2_start_in_q2_reads","avg_clipped_length","avg_distance_to_effective_3p_end ","base","count","avg_mapping_quality","avg_basequality","avg_se_mapping_quality","num_plus_strand","num_minus_strand","avg_pos_as_fraction","avg_num_mismatches_as_fraction","avg_sum_mismatch_qualities","num_q2_containing_reads","avg_distance_to_q2_start_in_q2_reads","avg_clipped_length","avg_distance_to_effective_3p_end ","base","count","avg_mapping_quality","avg_basequality","avg_se_mapping_quality","num_plus_strand","num_minus_strand","avg_pos_as_fraction","avg_num_mismatches_as_fraction","avg_sum_mismatch_qualities","num_q2_containing_reads","avg_distance_to_q2_start_in_q2_reads","avg_clipped_length","avg_distance_to_effective_3p_end ","base","count","avg_mapping_quality","avg_basequality","avg_se_mapping_quality","num_plus_strand","num_minus_strand","avg_pos_as_fraction","avg_num_mismatches_as_fraction","avg_sum_mismatch_qualities","num_q2_containing_reads","avg_distance_to_q2_start_in_q2_reads","avg_clipped_length","avg_distance_to_effective_3p_end ","base","count","avg_mapping_quality","avg_basequality","avg_se_mapping_quality","num_plus_strand","num_minus_strand","avg_pos_as_fraction","avg_num_mismatches_as_fraction","avg_sum_mismatch_qualities","num_q2_containing_reads","avg_distance_to_q2_start_in_q2_reads","avg_clipped_length","avg_distance_to_effective_3p_end ","base","count","avg_mapping_quality","avg_basequality","avg_se_mapping_quality","num_plus_strand","num_minus_strand","avg_pos_as_fraction","avg_num_mismatches_as_fraction","avg_sum_mismatch_qualities","num_q2_containing_reads","avg_distance_to_q2_start_in_q2_reads","avg_clipped_length","avg_distance_to_effective_3p_end")

M_gDNA1<-M_gDNA
M_gDNA=M_gDNA[(M_gDNA1[,5]=="="),]

colnames(M_gDNA)<-head1
#extract hit chr amplicon
M_gDNA_chr=table(M_gDNA[,1])
M_gDNA_chr=M_gDNA_chr[M_gDNA_chr>0]

# load SNP sites
SNPSiteAll=read.table(SNPsite,skip=1,sep="")
#count how many amplicons in SNP file
SNP_chr=table(SNPSiteAll[,1])

sum_report<-matrix(0,nrow(M_gDNA),25)
colnames(sum_report)=c("position","rs","Gene","CDS","chr","ref","coverage","Score","A","","C","","G","","T","","Ins",
                       "","Del","","","","","","")
sum_report[,1]=M_gDNA[,2]
sum_report[,2]=M_gDNA[,2]
sum_report[,5]=as.character(M_gDNA[,1])
sum_report[,6]=as.character(M_gDNA[,3])
sum_report[,7]=M_gDNA[,4]
sum_report[,9]=M_gDNA[,20]
sum_report[,11]=M_gDNA[,34]
sum_report[,13]=M_gDNA[,48]
sum_report[,15]=M_gDNA[,62]

#calculate ins
for (j in 1:nrow(M_gDNA)){
  for (i in 1:((ncol(M_gDNA)-88)/14)){
    if (!is.null(M_gDNA[j,((i-1)*14+89)])){
      if (M_gDNA[j,((i-1)*14+89)]!=""){
        if (regexpr("-",as.character(M_gDNA[j,(i-1)*14+89]))>0) sum_report[j,19]=as.numeric(sum_report[j,19])+as.numeric(if (is.na(M_gDNA[j,(i-1)*14+90])) 0 else M_gDNA[j,(i-1)*14+90])
        else sum_report[j,17]=as.numeric(sum_report[j,17])+as.numeric(if (is.na(M_gDNA[j,(i-1)*14+90])) 0 else M_gDNA[j,(i-1)*14+90])
      } else 0
    }else 0
    
  }
}

#add rs info
for (chr_num in 1:length(M_gDNA_chr)){
  
  sum_report_chrnum=sum_report[sum_report[,5]==names(M_gDNA_chr)[chr_num],]
  SNPSiteAll_chrnum=SNPSiteAll[SNPSiteAll[,1]==names(M_gDNA_chr)[chr_num],]
  
  sum_report_chrnum=na.omit(sum_report_chrnum)
  
  
  if (length(SNPSiteAll_chrnum)>1 ) {
    
    if( length(sum_report_chrnum)>1) {
      
      if( length(nrow(sum_report_chrnum))>0) {
        for (pos_n in 1:nrow(sum_report_chrnum)){
          
          if(sum_report_chrnum[pos_n,1] %in% SNPSiteAll_chrnum[,2]){
            #when SNVs
            #sum_report_chrnum[pos_n,2]=as.character(SNPSiteAll_chrnum[SNPSiteAll_chrnum[,3]==sum_report_chrnum[pos_n,1],4])
            sum_report_chrnum[pos_n,2] = as.character(SNPSiteAll_chrnum[SNPSiteAll_chrnum[,2]==sum_report_chrnum[pos_n,1],4])
            
          }else{
            #when multi-base 
            closest=which.min(abs(as.numeric(sum_report_chrnum[pos_n,1])-as.numeric(SNPSiteAll_chrnum[,2])))
            sum_report_chrnum[pos_n,2]=paste(SNPSiteAll_chrnum[closest,4],"+",sum_report_chrnum[pos_n,1],sep="")
            
          }
          
          }
                
        sum_report[sum_report[,5] %in% names(M_gDNA_chr)[chr_num],2]=sum_report_chrnum[,2]
        
      } else{
        sum_report_chrnum[2]=as.character(SNPSiteAll_chrnum[SNPSiteAll_chrnum[,3] %in% sum_report_chrnum[1],4])
        sum_report[sum_report[,5]==names(M_gDNA_chr)[chr_num],2]=sum_report_chrnum[2]
      }
      
    }
    
  }
  
}

sum_report=rbind(matrix(NA,3,25),colnames(sum_report),sum_report)

sum_report[1,1]<-c("bwa-samtools-bamreadcount")
sum_report[1,2]<-name

#save file as formatted
write.csv(sum_report,file=paste(name,".sumreport.csv",sep=""),row.names=FALSE)

### call the SNP change
library(tools)

#setting to compare difference of two sample, threshold to determine if a site is het or hom
sum_diff=0.8          
sum_same=0.1
sqsum_diff=0.5
sqsum_same=0.05

DNA1=sum_report

DNA1a=DNA1[5:nrow(DNA1),c(2,6,7,8,9,11,13,15,17,19,5,1)]

#calculate fraction of A C G T Ins and Del
DNA1b=matrix(0,nrow(DNA1a),ncol(DNA1a))
DNA1b[,1]=as.character(DNA1a[,1])
DNA1b[,2]=as.character(DNA1a[,2])
DNA1b[,3]=as.numeric(as.character(DNA1a[,3]))
DNA1b[,4]=as.numeric(as.character(DNA1a[,4]))
DNA1b[,5]=as.numeric(as.character(DNA1a[,5]))/as.numeric(as.character(DNA1a[,3]))
DNA1b[,6]=as.numeric(as.character(DNA1a[,6]))/as.numeric(as.character(DNA1a[,3]))
DNA1b[,7]=as.numeric(as.character(DNA1a[,7]))/as.numeric(as.character(DNA1a[,3]))
DNA1b[,8]=as.numeric(as.character(DNA1a[,8]))/as.numeric(as.character(DNA1a[,3]))
DNA1b[,9]=as.numeric(as.character(DNA1a[,9]))/as.numeric(as.character(DNA1a[,3]))
DNA1b[,10]=as.numeric(as.character(DNA1a[,10]))/as.numeric(as.character(DNA1a[,3]))
DNA1b[,11]=as.character(DNA1a[,11])
DNA1b[,12]=as.numeric(DNA1a[,12])
colnames(DNA1b)=c("rs","ref","coverage","score","A_frac","C_frac","G_frac","T_frac","Ins_frac","Del_frac","chr","pos")

DNA1c=DNA1b

#call DNA1
DNA_call1 <- matrix(0,nrow(DNA1c),3)
DNA_name<-c("A","C","G","T","I","D")

#call the variants
for (i in 1:nrow(DNA1c))  {
  DNA_sort<-sort.int(as.numeric(as.character(DNA1c[i,5:10])),decreasing = TRUE,index.return = TRUE)
  DNA_sort$ix[DNA_sort$x > min_call]
  DNA_call1[i,1]<- paste(DNA_name[DNA_sort$ix[DNA_sort$x > min_call]],collapse="")
  DNA_ref <- DNA_call1[i,1]==DNA1c[i,2]
  DNA_call1[i,2]<- if (DNA_call1[i,1]==as.character(DNA1c[i,2])) "=" else paste(DNA1c[i,2],">",DNA_call1[i,1],collapse="")
  #DNA_call1[i,3]<- DNA_call1[i,1]
  DNA_call1[i,3]=if (nchar(DNA_call1[i,1])==1) paste(DNA_call1[i,1],DNA_call1[i,1],sep="") else DNA_call1[i,1]
}

DNA1c=cbind(DNA1c,DNA_call1)
#filter min coverage

#save fraction result
write.csv(DNA1c,file=paste(name,".call.csv",sep=""),row.names=FALSE)

############################
#plot every amplicon

for (chr_num1 in 1:length(M_gDNA_chr)){
  
  sum_report_chr=sum_report[5:nrow(sum_report),]
  sum_report_chr=sum_report_chr[sum_report_chr[,5]==names(M_gDNA_chr)[chr_num1],]
  DNA1c_chr=DNA1c[DNA1c[,11]==names(M_gDNA_chr)[chr_num1],]
  
  if( length(sum_report_chr)>0) {
    if(length(nrow(sum_report_chr))>0){
      
      #jpeg(paste(name,"_",names(M_gDNA_chr)[chr_num1],".jpg",sep=""), width = 6, height = 6, units = "in", res = 300)
      dir.create(paste(dirname(name),"/",names(M_gDNA_chr)[chr_num1],sep=""),showWarnings =F)
      jpeg(paste(dirname(name),"/",names(M_gDNA_chr)[chr_num1],"/",basename(name),"_",names(M_gDNA_chr)[chr_num1],".jpg",sep=""),
           width = 6, height = 6, units = "in", res = 300)
      
      par(oma = c(5, 1, 1, 1))
      par(mfrow=c(2,1))
      par(mar=c(0.5,1,1,1))
      #reads plot
      barplot(t(sum_report_chr[,c(9,11,13,15,17,19)]),
              col=c('green','blue','orange','red','grey','wheat'),
              main=paste(basename(name),"_",names(M_gDNA_chr)[chr_num1],sep=""), cex.main=0.75,
              las=3,cex.names=0.5)
      abline(h=min_coverage,col='violet',lwd=0.5)
      
      par(mar=c(1.5,1,0.5,1))
      #fraction plot
      barplot(t(DNA1c_chr[,c(5,6,7,8,9,10)]),
              col=c('green','blue','orange','red','grey','wheat'),
              cex.main=0.75,
              names.arg=DNA1c_chr[,1],las=3,cex.names=0.75)
      abline(h=c((1:10)/10),col='white')
      
      par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
      legend("bottom", c('A','C','G','T','ins','del'), xpd = TRUE, horiz = TRUE, 
             inset = c(0,0), bty = "n", pch=15,col = c('green','blue','orange','red','grey','wheat'), 
             cex = 1)
      dev.off()
      
      #save.image()
      
    }
    
  }
  
}


#plot all SNP together
sum_report_chr=sum_report[5:nrow(sum_report),]

DNA1c_chr=DNA1c

if( length(sum_report_chr)>0) {
  
  jpeg(paste(name,".jpg",sep=""), width = 6, height = 6, units = "in", res = 600)
  
  par(oma = c(5, 1, 1, 1))
  par(mfrow=c(2,1))
  par(mar=c(0.5,1,1,1))
  #reads plot
  barplot(t(sum_report_chr[,c(9,11,13,15,17,19)]),
          col=c('green','blue','orange','red','grey','wheat'),
          main=paste(basename(name),sep=""), cex.main=0.75,border=NA,
          las=3,cex.names=0.5)
  abline(h=min_coverage,col='violet',lwd=0.5)
  
  par(mar=c(1.5,1,0.5,1))
  #fraction plot
  barplot(t(DNA1c_chr[,c(5,6,7,8,9,10)]),
          col=c('green','blue','orange','red','grey','wheat'),border=NA,
          cex.main=0.75,
          names.arg=DNA1c_chr[,1],las=3,cex.names=0.25)
  abline(h=c((1:10)/10),col='white',lwd=0.5)
  
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", c('A','C','G','T','ins','del'), xpd = TRUE, horiz = TRUE, 
         inset = c(0,0), bty = "n", pch=15,col = c('green','blue','orange','red','grey','wheat'), 
         cex = 1)
  dev.off()
  
}
