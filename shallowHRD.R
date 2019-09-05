##### SCRIPT : shallowHRD #####

args = commandArgs(trailingOnly=TRUE)

NAMEEE = args[1]
path = args[2]
hg19_cyto = args[3]

# NAMEEE = "D153R01"
# path = "/home/sternlab/Desktop/creation_nice_pipeline_shallowHRD/D153R01"
# hg19_cyto = "/home/sternlab/Desktop/creation_nice_pipeline_shallowHRD/cytoBand_adapted_fr_hg19.csv"

# setwd(paste(path, "/", NAMEEE,"/", sep = ""))
# 
# pPath<-paste(path, "/", NAMEEE,"/", sep = "")
# outputPath<-paste(pPath)
# inputPath<-paste(path, "/", NAMEEE,"/", sep = "")

setwd(paste(path, sep = ""))

pPath<-paste(path, sep = "")
outputPath<-paste(pPath)
inputPath<-paste(path, sep = "")


##### Avoid error problem in command line #####

continue_on_error <- function() { 
  print("on going...") 
}
options(error=continue_on_error) 


##### Cytoband_hg19 #####

# cytobAnnot<-"cytoBand_adapted_fr_hg19.csv"
cyt_Annot<-read.csv(hg19_cyto,header=T,fill=T)

cyt_Annot<-cyt_Annot[,1:6]
for(k in 2:dim(cyt_Annot)[[1]]){
  str1<-as.numeric(cyt_Annot[k-1,2:3])
  str2<-as.numeric(cyt_Annot[k,2:3])
  if(sum(abs(str1-str2))==0){cyt_Annot[k,1]<-0}
}

cyt_Annot[1:(dim(cyt_Annot)[[1]]-1),1]<-cyt_Annot[1:(dim(cyt_Annot)[[1]]-1),1]+cyt_Annot[2:dim(cyt_Annot)[[1]],1]
cyt_Annot[dim(cyt_Annot)[[1]],1]<-1
cyt_Annot<-cyt_Annot[-which(cyt_Annot[,1]==0),]
cyt_Annot[,1]<-cyt_Annot[,6]-cyt_Annot[,5]
cyt_Annot[,2]<-cyt_Annot[,2]+(cyt_Annot[,3]-1)/2

for(k in 2:dim(cyt_Annot)[[1]]){
  if(cyt_Annot[k,2]==cyt_Annot[k-1,2]){cyt_Annot[k,1]<-cyt_Annot[k,6]-cyt_Annot[k-1,5];cyt_Annot[k-1,1]<-0}
}

cyt_Annot<-cyt_Annot[-which(cyt_Annot[,1]==0),]
arm_len<-cyt_Annot[,1:2]
chr_len<-arm_len
chr_len[,2]<-trunc(chr_len[,2])
for(k in 2:dim(chr_len)[[1]]){
  if(chr_len[k,2]==chr_len[k-1,2]){chr_len[k,1]==chr_len[k-1,1]+chr_len[k,1]
    chr_len[k-1,1]<-0
  }
}
chr_len<-chr_len[-which(chr_len[,1]==0),]
chr_len<-chr_len[-dim(chr_len)[[1]],]
GG<-sum(chr_len[,1])
chr_p<-round(chr_len[,1]/GG,3)


##### LIBRARIES #####

# library("ggplot2")
# library("RColorBrewer")
# require("scales")
library("ggpubr")
library("gridExtra")
library("DescTools")


##### FUNCTIONS #####
## Function Local minima ##

localMinima <- function(x) { 
  y <- diff(c(.Machine$integer.max, x)) < 0L # diff un à un des y +  < 0
  rle(y)$lengths                             # number of chained TRUE or FALSE
  y <- cumsum(rle(y)$lengths)                # somme cumulé aux différents points : 1 7 10 13... 
  y <- y[seq.int(1L, length(y), 2L)]         # 
  if (x[[1]] == x[[2]]) {                    # retire des éléments de y
    y <- y[-1]
  }
  y                                          # retourne les entiers qui marchent 
}


## Function readSegmFile ##

readSegmFile<-function(segFileName){
  aa<-scan(segFileName,nlines=1,what="character",quiet = T)
  k<-1
  
  while((length(grep("#",aa))!=0)&&(grep("#",aa)==1)){
    aa<-scan(segFileName,skip=k,nlines=1,what="character",quiet = T)
    k<-k+1
  }
  
  tmp<-read.delim(segFileName,skip=k-1, header=T,sep="\t")
  tmp<-data.matrix(tmp)
  if(k==2){
    results<-getInfoSegm(scan(segFileName,nlines=1,what="character",quiet = T))
    homoConst<-results$homoConst
    sdH<-results$sdH	
    p_BAF<-results$p_BAF
    q_LRR<-results$q_LRR
    Delta<-results$Delta
  }else{homoConst<-1;sdH<-NA;p_BAF<-NA;q_LRR<-NA;Delta<-c(2,NA)}
  
  out<-list(tmp=tmp,homoConst=homoConst,sdH=sdH,p_BAF=p_BAF,q_LRR=q_LRR,Delta=Delta)
}


## getSegmentID : gather in the same name segment of the length  ##

getSegmentID<-function(THR,tmp,c_chr,c_cn,c_conf){
  
  tmp[,c_conf]<-0;tmp[1,c_conf]<-1
  
  for(k in 2:dim(tmp)[[1]]){
    
    if(tmp[k,c_chr]==tmp[k-1,c_chr] & abs(tmp[k,c_cn]-tmp[k-1,c_cn])<THR){
      
      tmp[k,c_conf]<-tmp[k-1,c_conf]
      
    }else{tmp[k,c_conf]<-tmp[k-1,c_conf]+1}
  } 
  
  out<-tmp

}


## shrinkReprTMP

shrinkReprTMP<-function(tmp,c_posS,c_posE,c_cn,c_conf){
  
  for(k in 1:(dim(tmp)[[1]]-1)){
    
    if(tmp[k,c_conf]==tmp[k+1,c_conf]){
      
      tmp[k+1,c_posS]<-tmp[k,c_posS]
      
      w<-c(tmp[k+1,c_posE]-tmp[k+1,c_posS],tmp[k,c_posE]-tmp[k,c_posS])
      
      tmp[k+1,c_cn]<-weighted.mean(c(tmp[k+1,c_cn],tmp[k,c_cn]),w)
      
      tmp[k,1]<-0
    }
  }
  
  tt<-which(tmp[,1]==0);if(length(tt)!=0){tmp<-tmp[-tt,]}
  out<-tmp

}


## Function getInfoSegm

getInfoSegm<-function(infoString){
  
  tt<-unlist(strsplit(infoString,";"))
  hh<-grep("homoConst=",tt)
  if(length(hh)==1){homoConst<-as.numeric(sub("homoConst=","",tt[hh]))}else{homoConst<-1}
  hh<-grep("sdH=",tt)
  if(length(hh)==1){sdH<-as.numeric(sub("sdH=","",tt[hh]))}else{sdH<-NA}
  
  hh<-grep("p=",tt)
  if(length(hh)==1){p_BAF<-as.numeric(sub("p=","",tt[hh]))}else{p_BAF<-NA}
  
  hh<-grep("q=",tt)
  if(length(hh)==1){q_LRR<-as.numeric(sub("q=","",tt[hh]))}else{q_LRR<-NA}
  
  hh<-grep("2LR=",tt)
  if(length(hh)==1){Delta<-as.numeric(sub("2LR=","",tt[hh]))}else{Delta<-NA}
  Delta<-c(2,Delta)
  
  
  if((is.na(homoConst))||(!is.numeric(homoConst))){homoConst<-1
  }else{if((homoConst<0.8)||(homoConst>1)){homoConst<-1}
  }
  if((is.na(sdH))||(!is.numeric(sdH))){sdH<-NA
  }else{if((sdH<0.05)||(sdH>1)){sdH<-NA}
  }
  
  out<-list(homoConst=homoConst,sdH=sdH,p_BAF=p_BAF,q_LRR=q_LRR,Delta=Delta)
  
}


## SmoothBreaks_SNP

smoothBreaksLength<-function(THR,Length,tmp,c_ind,c_chr,c_posS,c_posE,c_cn,c_CN,c_conf){

  tmp<-getSegmentID(THR=THR,tmp=tmp,c_chr=c_chr+1,c_cn=c_cn,c_conf=c_conf)
  tmp<-shrinkReprTMP(tmp=tmp,c_posS=c_posS,c_posE=c_posE,c_cn=c_cn,c_conf=c_conf)
  
  ss<-1
  for(ss in Length){
    
    tt<-which(round((tmp[,c_posE]-tmp[,c_posS])/10^6,1)<=ss)
    
    if(length(tt)>0){tmp<-tmp[-tt,]}
    
    
    for(k in 1:(dim(tmp)[[1]]-1)){
      
      
      if(tmp[k,c_chr+1]==tmp[k+1,c_chr+1] & abs(tmp[k,c_cn]-tmp[k+1,c_cn])<THR ){
        
        tmp[k+1,c_posS]<-tmp[k,c_posS]
        
        w<-c(tmp[k+1,c_posE]-tmp[k+1,c_posS],tmp[k,c_posE]-tmp[k,c_posS])
        
        tmp[k+1,c_cn]<-weighted.mean(c(tmp[k+1,c_cn],tmp[k,c_cn]),w)
        
        tmp[k,c_ind]<-0
        
      }
    }
    
    tt<-which(tmp[,c_ind]==0)
    if(length(tt)>0){tmp<-tmp[-tt,]}
  }
  
  
  tmp<-getSegmentID(THR=THR,tmp=tmp,c_chr=c_chr+1,c_cn=c_cn,c_conf=c_conf)
  tmp<-shrinkReprTMP(tmp=tmp,c_posS=c_posS,c_posE=c_posE,c_cn=c_cn,c_conf=c_conf)
  
  out<-tmp
  
}


## breakSmoothToLST

breakSmoothToLST<-function(THR,tmp,c_ind,c_chr,c_posS,c_posE,c_cn,c_CN,c_conf){
  
  
  tmp<-getSegmentID(THR=THR,tmp=tmp,c_chr=c_chr+1,c_cn=c_cn,c_conf=c_conf)
  tmp<-shrinkReprTMP(tmp=tmp,c_posS=c_posS,c_posE=c_posE,c_cn=c_cn,c_conf=c_conf)
  
  tmp[,c_ind]<-seq(1,dim(tmp)[[1]])
  
  BreaksSmall<-0
  
  FL<-T
  while(FL){
    
    tmp[,c_ind]<-seq(1,dim(tmp)[[1]])
    
    kk<-which(round((tmp[,c_posE]-tmp[,c_posS])/10^6,1)<3) # segments inférieurs à 3Mb
    
    #		kk<-union(kk,which(tmp[,c_bin]<50))
    
    if(length(kk)>0){
      
      kk<-kk[order((tmp[kk,c_posE]-tmp[kk,c_posS])/10^6)] # ordonne les segments inférieurs à 3Mb
      
      
      for(k in 1:length(kk)){
        
        if(kk[k]==1){                                        # un seul segment < 3Mb
          if(tmp[kk[k]+1,c_chr+1]==tmp[kk[k],c_chr+1]){      # si segment d'après même chr_arm
            if(tmp[kk[k]+1,c_ind]!=0){                       # si l'index d'après différent de 0
              tmp[kk[k],c_ind]<-0                            # alors on met un index 0   
            }
          }                          
        }
        
        else{
          if(kk[k]==dim(tmp)[[1]]){                                                # tous les segments concernés
            if(tmp[kk[k]-1,c_chr+1]==tmp[kk[k],c_chr+1] & tmp[kk[k]-1,c_ind]!=0){  # si chrom_arm d'avant même et index différent de 0  
              tmp[kk[k],c_ind]<-0}                                                 # on met index =0
          }
          else{
            
            # si chr_arm d'avant et d'après même chr_arm et index différent de 0, on met index = 0
            if(tmp[kk[k]-1,c_chr+1]==tmp[kk[k]+1,c_chr+1] & tmp[kk[k]-1,c_ind]!=0 & tmp[kk[k]+1,c_ind]!=0){tmp[kk[k],1]<-0}
            
            # si chr_arm avant même chr_arm et après différent, et que avant index différent de 0, on met index = 0
            if(tmp[kk[k]-1,c_chr+1]==tmp[kk[k],c_chr+1] & tmp[kk[k]+1,c_chr+1]!=tmp[kk[k],c_chr+1] & tmp[kk[k]-1,c_ind]!=0){tmp[kk[k],c_ind]<-0}
            
            # si chr_arm après même chr_arm et avant différent et que après index différent de 0, on met index = 0 
            if(tmp[kk[k]+1,c_chr+1]==tmp[kk[k],c_chr+1] & tmp[kk[k]-1,c_chr+1]!=tmp[kk[k],c_chr+1] & tmp[kk[k]+1,c_ind]!=0){tmp[kk[k],c_ind]<-0}
            
            # si chr_arm après et avant différent, alors on met 0 
            if(tmp[kk[k]+1,c_chr+1]!=tmp[kk[k],c_chr+1] & tmp[kk[k]-1,c_chr+1]!=tmp[kk[k],c_chr+1]){tmp[kk[k],c_ind]<-0}
            
          }
        }
      } # assignation de 0 à tout segments inférieurs à 0 possibles ceux possibles 
      
      tt<-which(tmp[,c_ind]==0) # tous les segments inférieurs à 3 Mb qui ont été assigné à 0
      
      if(length(tt)>0){tmp <- tmp[-tt,];BreaksSmall<-BreaksSmall+length(tt)} # s'il y a des segments dans ce cas, retiré de tmp
      # Breaksmall prends le nombre de segment dans ce cas
      
      for(k in 1:(dim(tmp)[[1]]-1)){ # parcours nouveau tmp
        
        if(round((tmp[k+1,c_posS]-tmp[k,c_posE])/10^6,1)<3){    # arrondi de la fin du segment d'après moins début de celui d'avant 
          # inférieur à 3Mb
          
          if(tmp[k,c_chr+1]==tmp[k+1,c_chr+1] & abs(tmp[k,c_cn]-tmp[k+1,c_cn])<THR){ # si segment et après même chr_arm, et diff < THR
            
            tmp[k+1,c_posS]<-tmp[k,c_posS]          # segment start k+1 prend start k d'avant
            
            w<-c(tmp[k+1,c_posE]-tmp[k+1,c_posS],tmp[k,c_posE]-tmp[k,c_posS]) # c(distance fin k+1 début k+1, taille k)
            
            tmp[k+1,c_cn]<-weighted.mean(c(tmp[k+1,c_cn],tmp[k,c_cn]),w) # ratio_median = moyenne pondéré en fonction facteur taille
            
            tmp[k,c_ind]<-0                                              # index de k = 0
            
          }
          
          
        }
      }
      
      tt<-which(tmp[,c_ind]==0);if(length(tt)>0){tmp<-tmp[-tt,]} # retire segment ainsi lier
      
      
    }else{FL<-F}  
  }
  
  
  
  tmp<-getSegmentID(THR=THR,tmp=tmp,c_chr=c_chr+1,c_cn=c_cn,c_conf=c_conf)
  tmp<-shrinkReprTMP(tmp=tmp,c_posS=c_posS,c_posE=c_posE,c_cn=c_cn,c_conf=c_conf)
  
  tmp[1,c_ind]<-BreaksSmall
  
  out<-tmp
  
  
}


## LST_control

LST_control<-function(THR,lenBIN,lenMB,tmp,c_ind,c_chr,c_posS,c_posE,c_cn,c_CN,c_conf){

  tmp[,c_ind]<-seq(1,dim(tmp)[[1]])
  
  
  WC<-tmp[,c_CN]*0
  
  j<-5
  
  for(k in 2:dim(tmp)[[1]]){
    
    if(tmp[k,c_chr+1]==tmp[k-1,c_chr+1] & round((tmp[k,c_posS]-tmp[k-1,c_posE])/10^6,1)<3){
      
      if(round((tmp[k,c_posE]-tmp[k,c_posS])/10^6,0)>=lenMB & round((tmp[k-1,c_posE]-tmp[k-1,c_posS])/10^6,0)>=lenMB){
        
        
        if(abs(tmp[k,c_cn]-tmp[k-1,c_cn])>THR*coefficient){WC[k]<-1}
        
      }

    }
  }
  
  
  out<-WC
  
  
}


## ShortIntestBreaksSmooth

ShortIntestBreaksSmooth<-function(THR,lenMB,tmp,c_ind=c_ind,c_chr=c_chr,c_posS=c_posS,c_posE=c_posE,c_cn=c_cn,c_CN=c_CN,c_conf=c_conf){
  
  
  segSize<-round((tmp[,c_posF]-tmp[,c_posS])/10^6,1)
  
  
  for(k in 2:(dim(tmp)[[1]]-1)){
    
    LengthOut<-max(3,segSize[k]*2)
    
    if(tmp[k+1,c_chr]==tmp[k-1,c_chr] & segSize[k]<lenMB & segSize[k-1]>LengthOut & segSize[k+1]>LengthOut){
      
      
      if(abs(tmp[k+1,c_cn]-tmp[k-1,c_cn])<THR){
        
        tmp[k+1,c_posS]<-tmp[k-1,c_posS]
        
        w<-c(tmp[k+1,c_posE]-tmp[k+1,c_posS],tmp[k-1,c_posE]-tmp[k-1,c_posS])
        
        tmp[k+1,c_cn]<-weighted.mean(c(tmp[k+1,c_cn],tmp[k-1,c_cn]),w)
        
        tmp[c(k-1,k),c_ind]<-0
        
        
      }
    }
  }
  
  
  tt<-which(tmp[,c_ind]==0);if(length(tt)>0){tmp<-tmp[-tt,]}
  
  tmp[1,c_ind]<-length(tt)
  
  out<-tmp
  
  
}


##### Fast ControlFREEC gathering #####

dataTable <- read.table(paste(pPath,"/",NAMEEE,".bam_ratio.txt", sep = ""), header = TRUE)
dataTable = dataTable[,-5]

dataTable = dataTable[which(!dataTable$Chromosome == "X"),] # remove X
dataTable = dataTable[which(!dataTable$Chromosome == "Y"),] # remove Y

dataTable[,1] = as.numeric(as.character(dataTable[,1]))

dataTable = dataTable[order(dataTable[,1]),]

L = dim(dataTable)[1]

size_window = dataTable[2,2] - dataTable[1,2]

dataTable = cbind(1:L, dataTable[,1], dataTable[,2], dataTable[,2]+size_window-1, 
                  dataTable[,3], dataTable[,4]) 

dataTable = as.data.frame(dataTable)
colnames(dataTable) = c("feature", "chromosome", "start", "end", "ratio", "ratio_median")

dataTable = dataTable[which(!dataTable$ratio == -1),]
dataTable = dataTable[which(!dataTable$ratio_median == -1),]

dataTable[,5] = log2(dataTable[,5])
dataTable[,6] = log2(dataTable[,6])

ratio <- data.frame(dataTable)

short_size_window = size_window/1000

write.table(ratio, file = paste("Ratio_",NAMEEE,"_",short_size_window,"kb.tsv", sep = "") , row.names = FALSE, col.names = TRUE, sep = "\t")

X=read.table(paste("Ratio_",NAMEEE,"_",short_size_window,"kb.tsv", sep = ""), header=TRUE, sep = "\t")            

X=X[,-1] # feature
X=X[,-4] # ratio

## Remove spurious regions based on telomere and centromere from UCSC

X <- X[which(!(X[,1]==1 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==1 & X[,2]>=121535434 & X[,3]<=124535434)),]
X <- X[which(!(X[,1]==1 & X[,2]>=249240621 & X[,3]<=249250621)),]
X <- X[which(!(X[,1]==2 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==2 & X[,2]>=92326171 & X[,3]<=95326171)),]
X <- X[which(!(X[,1]==2 & X[,2]>=243189373 & X[,3]<=243199373)),]
X <- X[which(!(X[,1]==3 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==3 & X[,2]>=90504854 & X[,3]<=93504854)),]
X <- X[which(!(X[,1]==3 & X[,2]>=198012430 & X[,3]<=198022430)),]
X <- X[which(!(X[,1]==4 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==4 & X[,2]>=49660117 & X[,3]<=52660117)),]
X <- X[which(!(X[,1]==4 & X[,2]>=191144276 & X[,3]<=191154276)),]
X <- X[which(!(X[,1]==5 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==5 & X[,2]>=46405641 & X[,3]<=9405641)),]
X <- X[which(!(X[,1]==5 & X[,2]>=180905260 & X[,3]<=180915260)),]
X <- X[which(!(X[,1]==6 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==6 & X[,2]>=58830166 & X[,3]<=61830166)),]
X <- X[which(!(X[,1]==6 & X[,2]>=171105067 & X[,3]<=171115067)),]
X <- X[which(!(X[,1]==7 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==7 & X[,2]>=58054331 & X[,3]<=61054331)),]
X <- X[which(!(X[,1]==7 & X[,2]>=159128663 & X[,3]<=159138663)),]
X <- X[which(!(X[,1]==8 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==8 & X[,2]>=43838887 & X[,3]<=46838887)),]
X <- X[which(!(X[,1]==8 & X[,2]>=146354022 & X[,3]<=146364022)),]
X <- X[which(!(X[,1]==9 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==9 & X[,2]>=47367679 & X[,3]<=50367679)),]
X <- X[which(!(X[,1]==9 & X[,2]>=141203431 & X[,3]<=141213431)),]
X <- X[which(!(X[,1]==10 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==10 & X[,2]>=39254935 & X[,3]<=42254935)),]
X <- X[which(!(X[,1]==10 & X[,2]>=135524747 & X[,3]<=135534747)),]
X <- X[which(!(X[,1]==11 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==11 & X[,2]>=51644205 & X[,3]<=54644205)),]
X <- X[which(!(X[,1]==11 & X[,2]>=134996516 & X[,3]<=135006516)),]
X <- X[which(!(X[,1]==12 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==12 & X[,2]>=34856694 & X[,3]<=37856694)),]
X <- X[which(!(X[,1]==12 & X[,2]>=133841895 & X[,3]<=133851895)),]
X <- X[which(!(X[,1]==13 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==13 & X[,2]>=16000000 & X[,3]<=19000000)),]
X <- X[which(!(X[,1]==13 & X[,2]>=115159878 & X[,3]<=115169878)),]
X <- X[which(!(X[,1]==14 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==14 & X[,2]>=16000000 & X[,3]<=19000000)),]
X <- X[which(!(X[,1]==14 & X[,2]>=107339540 & X[,3]<=107349540)),]
X <- X[which(!(X[,1]==15 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==15 & X[,2]>=17000000 & X[,3]<=20000000)),]
X <- X[which(!(X[,1]==15 & X[,2]>=102521392 & X[,3]<=102531392)),]
X <- X[which(!(X[,1]==16 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==16 & X[,2]>=35335801 & X[,3]<=38335801)),]
X <- X[which(!(X[,1]==16 & X[,2]>=90344753 & X[,3]<=90354753)),]
X <- X[which(!(X[,1]==17 & X[,2]>=22263006 & X[,3]<=25263006)),]
X <- X[which(!(X[,1]==18 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==18 & X[,2]>=15460898 & X[,3]<=18460898)),]
X <- X[which(!(X[,1]==18 & X[,2]>=78067248 & X[,3]<=78077248)),]
X <- X[which(!(X[,1]==19 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==19 & X[,2]>=24681782 & X[,3]<=27681782)),]
X <- X[which(!(X[,1]==19 & X[,2]>=59118983 & X[,3]<=59128983)),]
X <- X[which(!(X[,1]==20 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==20 & X[,2]>=26369569 & X[,3]<=29369569)),]
X <- X[which(!(X[,1]==20 & X[,2]>=63015520 & X[,3]<=63025520)),]
X <- X[which(!(X[,1]==21 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==21 & X[,2]>=11288129 & X[,3]<=14288129)),]
X <- X[which(!(X[,1]==21 & X[,2]>=48119895 & X[,3]<=48129895)),]
X <- X[which(!(X[,1]==22 & X[,2]>=0 & X[,3]<=10000)),]
X <- X[which(!(X[,1]==22 & X[,2]>=13000000 & X[,3]<=16000000)),]
X <- X[which(!(X[,1]==22 & X[,2]>=51294566 & X[,3]<=51304566)),]


## Add chromosome arm 

options(show.error.messages = FALSE)

X$chr_arm <- rep(0,nrow(X))
colnames(X) <- c("chr", "start", "end", "ratio_median", "chr_arm")
X=X[c("chr", "chr_arm", "start", "end", "ratio_median")]

X[,3]=as.numeric(as.character(X[,3])) # start
X[,4]=as.numeric(as.character(X[,4])) # end

tt <- X[,1] == 1 & X[,3] < 125000000 & X[,4] < 125000000
X[tt,2] = 1
tt <- X[,1] == 1 & X[,3] > 125000000 & X[,4] > 125000000
X[tt,2] = 2
tt <- X[,1] == 1 & X[,3] < 125000000 & X[,4] > 125000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 125000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 125000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])  
rownames(X) <- NULL 

tt <- X[,1] == 2 & X[,3] < 93300000 & X[,4] < 93300000
X[tt,2] = 1
tt <- X[,1] == 2 & X[,3] > 93300000 & X[,4] > 93300000
X[tt,2] = 2
tt <- X[,1] == 2 & X[,3] < 93300000 & X[,4] > 93300000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 93300000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 93300001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 3 & X[,3] < 91000000 & X[,4] < 91000000
X[tt,2] = 1
tt <- X[,1] == 3 & X[,3] > 91000000 & X[,4] > 91000000
X[tt,2] = 2
tt <- X[,1] == 3 & X[,3] < 91000000 & X[,4] > 91000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 91000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 91000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 4 & X[,3] < 50400000 & X[,4] < 50400000
X[tt,2] = 1
tt <- X[,1] == 4 & X[,3] > 50400000 & X[4] > 50400000
X[tt,2] = 2
tt <- X[,1] == 4 & X[,3] < 50400000 & X[,4] > 50400000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 50400000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 50400001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 5 & X[,3] < 48400000 & X[,4] < 48400000
X[tt,2] = 1
tt <- X[,1] == 5 & X[,3] > 48400000 & X[,4] > 48400000
X[tt,2] = 2
tt <- X[,1] == 5 & X[,3] < 48400000 & X[,4] > 48400000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 48400000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 48400001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 6 & X[,3] < 61000000 & X[,4] < 61000000
X[tt,2] = 1
tt <- X[,1] == 6 & X[,3] > 61000000 & X[,4] > 61000000
X[tt,2] = 2
tt <- X[,1] == 6 & X[,3] < 61000000 & X[,4] > 61000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 61000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 61000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 7 & X[,3] < 59900000 & X[,4] < 59900000
X[tt,2] = 1
tt <- X[,1] == 7 & X[,3] > 59900000 & X[,4] > 59900000
X[tt,2] = 2
tt <- X[,1] == 7 & X[,3] < 59900000 & X[,4] > 59900000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 59900000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 59900001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 8 & X[,3] < 45600000 & X[,4] < 45600000
X[tt,2] = 1
tt <- X[,1] == 8 & X[,3] > 45600000 & X[,4] > 45600000
X[tt,2] = 2
tt <- X[,1] == 8 & X[,3] < 45600000 & X[,4] > 45600000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 45600000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 45600001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 9 & X[,3] < 49000000 & X[,4] < 49000000
X[tt,2] = 1
tt <- X[,1] == 9 & X[,4] > 49000000 & X[,4] > 49000000
X[tt,2] = 2
tt <- X[,1] == 9 & X[,3] < 49000000 & X[,4] > 49000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 49000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 49000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 10 & X[,3] < 40200000 & X[,4] < 40200000
X[tt,2] = 1
tt <- X[,1] == 10 & X[,3] > 40200000 & X[,4] > 40200000
X[tt,2] = 2
tt <- X[,1] == 10 & X[,3] < 40200000 & X[,4] > 40200000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 40200000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 40200000, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 11 & X[,3] < 53700000 & X[,4] < 53700000
X[tt,2] = 1
tt <- X[,1] == 11 & X[,3] > 53700000 & X[,4] > 53700000
X[tt,2] = 2
tt <- X[,1] == 11 & X[,3] < 53700000 & X[,4] > 53700000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 53700000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 53700001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 12 & X[,3] < 35800000 & X[,4] < 35800000
X[tt,2] = 1
tt <- X[,1] == 12 & X[,3] > 35800000 & X[,4] > 35800000
X[tt,2] = 2
tt <- X[,1] == 12 & X[,3] < 35800000 & X[,4] > 35800000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 35800000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 35800001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 13 & X[,3] < 17900000 & X[,4] < 17900000
X[tt,2] = 1
tt <- X[,1] == 13 & X[,3] > 17900000 & X[,4] > 17900000
X[tt,2] = 2
tt <- X[,1] == 13 & X[,3] < 17900000 & X[,4] > 17900000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 17900000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 17900001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 14 & X[,3] < 17600000 & X[,4] < 17600000
X[tt,2] = 1
tt <- X[,1] == 14 & X[,3] > 17600000 & X[,4] > 17600000
X[tt,2] = 2
tt <- X[,1] == 14 & X[,3] < 17600000 & X[,4] > 17600000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 17600000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 17600001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 15 & X[,3] < 19000000 & X[,4] < 19000000
X[tt,2] = 1
tt <- X[,1] == 15 & X[,3] > 19000000 & X[,4] > 19000000
X[tt,2] = 2
tt <- X[,1] == 15 & X[,3] < 19000000 & X[,4] > 19000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 19000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 19000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 16 & X[,3] < 36600000 & X[,4] < 36600000
X[tt,2] = 1
tt <- X[,1] == 16 & X[,3] > 36600000 & X[,4] > 36600000
X[tt,2] = 2
tt <- X[,1] == 16 & X[,3] < 36600000 & X[,4] > 36600000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 36600000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 36600001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 17 & X[,3] < 24000000 & X[,4] < 24000000
X[tt,2] = 1
tt <- X[,1] == 17 & X[,3] > 24000000 & X[,4] > 24000000
X[tt,2] = 2
tt <- X[,1] == 17 & X[,3] < 24000000 & X[,4] > 24000000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 24000000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 24000001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 18 & X[,3] < 17200000 & X[,4] < 17200000
X[tt,2] = 1
tt <- X[,1] == 18 & X[,3] > 17200000 & X[,4] > 17200000
X[tt,2] = 2
tt <- X[,1] == 18 & X[,3] < 17200000 & X[,4] > 17200000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 17200000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 17200001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 19 & X[,3] < 26500000 & X[,4] < 26500000
X[tt,2] = 1
tt <- X[,1] == 19 & X[,3] > 26500000 & X[,4] > 26500000
X[tt,2] = 2
tt <- X[,1] == 19 & X[,3] < 26500000 & X[,4] > 26500000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 26500000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 26500001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 20 & X[,3] < 27500000 & X[,4] < 27500000
X[tt,2] = 1
tt <- X[,1] == 20 & X[,3] > 27500000 & X[,4] > 27500000
X[tt,2] = 2
tt <- X[,1] == 20 & X[,3] < 27500000 & X[,4] > 27500000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 27500000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 27500001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 21 & X[,3] < 13200000 & X[,4] < 13200000
X[tt,2] = 1
tt <- X[,1] == 21 & X[,3] > 13200000 & X[,4] > 13200000
X[tt,2] = 2
tt <- X[,1] == 21 & X[,3] < 13200000 & X[,4] > 13200000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 13200000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 13200001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

tt <- X[,1] == 22 & X[,3] < 14700000 & X[,4] < 14700000
X[tt,2] = 1
tt <- X[,1] == 22 & X[,3] > 14700000 & X[,4] > 14700000
X[tt,2] = 2
tt <- X[,1] == 22 & X[,3] < 14700000 & X[,4] > 14700000
X[tt,2] = 1
B = X[tt,4]
X[tt,4] = 14700000
X=rbind(X[(1:which(tt==TRUE)),], c(X[tt,1], 2, 14700001, B, X[tt,5]) , X[-(1:which(tt==TRUE)),])
rownames(X) <- NULL 

X[,2]<-X[,1]+(X[,2]-1)/2      # tranform into 1 & 1.5

tt=which(X$ratio_median==-Inf)

if (length(tt) >= 1){
  X=X[-tt,]  
}

X=as.matrix(X)

L = dim(X)[1]
A = matrix(0, ncol=6, nrow=L)

i=1    # ligne pour regarder dans X
c=1    # ligne pour stocker bonne valeur dans A

while (i < L){
  if(X[i,2]==X[i+1,2]){             # Egalité chrom_arm ligne i et i+1 
    if (X[i,5]==X[i+1,5]){          # Egalité ratio_median ligne i et i+1
      n=1
      while (X[i,2]==X[i+n,2] & X[i,5]==X[i+n,5] & i+n < L){
        n=n+1
      }
      if(i+n == L){
        A[c,]=c(X[i,1], X[i,2], X[i,3], X[i+n,4], X[i,5], X[i+n,4]-X[i,3]+1) 
        i=i+n
        c=c+1
      }
      else{
        A[c,]=c(X[i,1], X[i,2], X[i,3], X[i+n-1,4], X[i,5], X[i+n-1,4]-X[i,3]+1) 
        i=i+n
        c=c+1
      }
    }
    else{
      A[c,]=c(X[i,1], X[i,2], X[i,3], X[i,4], X[i,5], X[i,4]-X[i,3]+1)
      i=i+1
      c=c+1
    }
  }
  else{
    A[c,]=c(X[i,1], X[i,2], X[i,3], X[i,4], X[i,5], X[i,4]-X[i,3]+1)
    i=i+1                                        # avancement de i+1 si chr_arm différents
    c=c+1                         
  }  
}

A=subset(A, A[,1] != 0)
rownames(A) <- NULL 
colnames(A) <- c("chr", "chr_arm", "start", "end", "ratio_median", "size")

write.table(A, file = paste(NAMEEE,"_medianratio_gatheredR_automatized_final.txt", sep = ""), sep = "\t", row.names = FALSE)

options(show.error.messages = TRUE)

##### Find Threshold #####

X = read.table(paste(NAMEEE,"_medianratio_gatheredR_automatized_final.txt", sep = ""), sep = "\t", header = TRUE)
X = X[which(X[,6] > 2999999),]

X = X[which(X[,6] > ((quantile(X[,6])[[4]] - quantile(X[,6])[[2]])/2)),]   # find the best size for the density plot

L=dim(X)[1]
X=data.matrix(X)

test=c()
i=1
for (i in 1:L){
  v=i+1
  for (y in v:L-1){
    test = c(test, abs(X[i,5] - X[y,5]))}
}  

minx = localMinima(density(test)$y)[2]
Threshold = density(test)$x[minx]

test = as.data.frame(test)

A = data.frame(xstart = c(0, 0.185544223315942), xend = c(0.185544223315942, 0.3))

ploty <- ggplot(test, aes(x = test)) + geom_density() +
  geom_vline(aes(xintercept = Threshold), color = "blue") + 
  ggtitle(paste("Density pairwise difference large segment")) + xlab("difference") +
  theme(plot.title = element_text(hjust = 0.5, size = 13),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=13),
        axis.text.y = element_text(size=15),
        axis.title.y = element_blank(),
        panel.background = element_blank())
ggsave(paste(NAMEEE,"_THR", sep = ""), plot = ploty, device = "jpeg", width = 20, height = 10)


##### Reading and initialisation #####

segFiles <- list.files(inputPath, pattern="automatized_final.txt",full.names=T) 

fileNames <- gsub(inputPath,"",segFiles) 
fileNames <- gsub("proc_","",fileNames) 
fileNames <- gsub(".txt","",fileNames) 
fileNames <- gsub("/","",fileNames) 

nSample<-1
nSample<-length(segFiles)

coefficient = 1.0

ColNamesTMP<-c("index", "chr","chr_arm","posStart", "posEnd", "ratio")

c_ind<-1; c_chr<-2;c_posS<-4; c_posE<-5; c_cn<-6; c_conf<-8

LLBB<-c(3,4,5,6,7,8,9,10,11) # what size LST

TAB<-matrix(0,nSample,20)              # matrix to stock results (colonne 20)
rownames(TAB)<-seq(1:nSample)          # Nombre de rangée = nombre de sample
colnames(TAB)<-seq(1,20)               # colnames de 1 à 20

colnames(TAB)[1:6]<-c("p_BAF","q_LRR","2copyLRR","DNA_ind","All_breaks","less3Mb_breaks")

ii<-1

rownames(TAB)[ii]<-fileNames[ii] # place the name of file in TAB rowname de TAB 

results<-readSegmFile(segFileName=segFiles[ii])  # lecture du fichier ( bcq cols inutiles)
tmp<-results$tmp                                 # Stocke partie résultat interressant 

## ajoute les index et segments

tmp <- cbind(seq(1,dim(tmp)[[1]]),tmp) # add index
colnames(tmp)[1]<-c("index")

tmp <- cbind(tmp,rep(0,nrow(tmp)))# add segment
colnames(tmp)[8]<-c("level")

write.table(tmp, file = paste(NAMEEE,"_",short_size_window,"kb_I", sep=""), sep = "\t", row.names = FALSE)  # gather


tt<-which(tmp[,c_chr+1]>23.6)                            # exclusion chr24 et plus 
if(length(tt)>0){tmp<-tmp[-tt,]}

tt<-which(tmp[,c_chr+1]==21)                             # exclusion partie short arm chr21
if(length(tt)>0){tmp<-tmp[-tt,]}

tt<-which(tmp[,c_chr+1]==22)                             # exclusion short arm chr22
if(length(tt)>0){tmp<-tmp[-tt,]}

write.table(tmp, file = paste(NAMEEE,"_",short_size_window,"kb_II", sep=""), sep = "\t", row.names = FALSE) # no chr_amr p 21 22 et plus de 24

tmp[,7] = tmp[,5] - tmp[,4] + 1

tmp_3mb = tmp[which(tmp[,7] > 2999999),] # segment de plus de 3 Mb
cop_tmp = tmp_3mb                       # copie du fichier des segments de plus de 3Mb
THR=Threshold

if (THR > 0.45){
  THR = 0.28}

if (THR < 0.025){
  THR = 0.025
}

level=1


##### Gather big segment (> 3Mb) #####
options(show.error.messages = FALSE)

while (dim(cop_tmp)[1] > 1){  
  
  line_largest_index = which.max(cop_tmp[,7])             # ligne avec segment le plus grand
  largest_index = cop_tmp[line_largest_index,1]           # index associé (change quand on remove)
  ratio_largest = cop_tmp[line_largest_index,6]           # median_ratio segment
  ratio_largest = ratio_largest + 0.00001                 # avoid picking same value
  
  close_index = cop_tmp[which.min(abs(cop_tmp[,6] - ratio_largest)),1]
  closest_index = c(largest_index, close_index)
  
  ## Find all value that matches the two closest index fit with THR 
  
  L=dim(cop_tmp)[1]
  
  for (i in 1:L){
    validation = 0
    n = length(closest_index)
    for (g in 1:n){
      if (abs(cop_tmp[i,6] - cop_tmp[which((cop_tmp[,1]==closest_index[g])),6]) < THR){   # 
        validation = validation + 1}
    }
    if (validation == n){
      closest_index = c(closest_index, cop_tmp[i,1])}}
  
  closest_index = unique(closest_index)
  
  ##  Assign level based on the index of close segment in genome    
  
  n = length(closest_index)
  for (i in 1:n){
    tmp_3mb[which((tmp_3mb[,1] == closest_index[i])),8] = level
    cop_tmp <- cop_tmp[which(!(cop_tmp[,1]==closest_index[i])),]
  }
  dim(tmp_3mb)
  dim(cop_tmp)
  level = level + 1
} 

##  Gather with level

options(show.error.messages = TRUE)

L1 = dim(tmp_3mb)[1]

A = matrix(0, ncol=8, nrow=L1)

i=1
c=1   # 2
while (i<L1+1){
  if (i==L1){
    A[c,]=c(tmp_3mb[i,1], tmp_3mb[i,2], tmp_3mb[i,3], tmp_3mb[i,4], tmp_3mb[i,5], tmp_3mb[i,6], tmp_3mb[i,5]-tmp_3mb[i,4]+1, tmp_3mb[i,8])    # segment pas sur deux lignes, a stocker seul 
    i=i+1
    c=c+1}
  else{
    if(tmp_3mb[i,3] == tmp_3mb[i+1,3]){                                           # Egalité de chr_arm
      if (tmp_3mb[i,8] == tmp_3mb[i+1,8]){                                        # Egalite de level
        n=1
        somme=tmp_3mb[i,6]
        while (tmp_3mb[i,8] == tmp_3mb[i+n,8] && tmp_3mb[i,3] == tmp_3mb[i+n,3]){
          somme = somme + tmp_3mb[i+n,6]
          n = n + 1
          if(i+n == L1+1){
            break}}
        A[c,]=c(tmp_3mb[i,1], tmp_3mb[i,2], tmp_3mb[i,3], tmp_3mb[i,4], tmp_3mb[i+n-1,5], (somme)/n, tmp_3mb[i+n-1,5]- tmp_3mb[i,4]+1, tmp_3mb[i,8])
        i=i+n
        c=c+1}
      else{                                             # pas d'égalité de niveau, on stocke directement la ligne
        A[c,]=c(tmp_3mb[i,1], tmp_3mb[i,2], tmp_3mb[i,3], tmp_3mb[i,4], tmp_3mb[i,5], tmp_3mb[i,6], tmp_3mb[i,5]-tmp_3mb[i,4]+1, tmp_3mb[i,8])    # segment pas sur deux lignes, a stocker seul 
        i=i+1
        c=c+1}}
    else{                                              # pas egalite des chromosomes entre i et i+1, on stocke la ligne i
      A[c,]=c(tmp_3mb[i,1], tmp_3mb[i,2], tmp_3mb[i,3], tmp_3mb[i,4], tmp_3mb[i,5], tmp_3mb[i,6], tmp_3mb[i,5]-tmp_3mb[i,4]+1, tmp_3mb[i,8])
      i=i+1                                                              
      c=c+1}
  }  
}

A=subset(A, A[,1] != 0)

rownames(A) <- NULL 
colnames(A) <- c("index", "chr", "chr_arm", "start", "end", "ratio_median", "size", "level")

tmp_3mb=A
write.table(tmp_3mb, file = paste(NAMEEE,"_",short_size_window,"kb_III", sep=""), sep = "\t", row.names = FALSE) 


##### Reput small segment & Smoothing #####

options(show.error.messages = FALSE)

tmp_0.1_3mb=tmp[which(tmp[,7] > 99999),]
tmp_0.1_3mb=tmp_0.1_3mb[which(tmp_0.1_3mb[,7] < 2999999),]  # pb ?

L_3mb=dim(tmp_3mb)[1]
l_0.1_3mb=dim(tmp_0.1_3mb)[1]

values = unique(tmp_0.1_3mb[,3][!tmp_0.1_3mb[,3] %in% tmp_3mb[,3]]) # check all chrm_arm that are not in 
length=length(values)

for (missing_chr_arm in values){
  c=1
  L_3mb=dim(tmp_3mb)[1]
  while(c < L_3mb+1){
    if (missing_chr_arm > max(tmp_3mb[,3])){
      tmp_3mb=rbind(tmp_3mb, tmp_0.1_3mb[which(tmp_0.1_3mb[,3] == missing_chr_arm),])}
    else if (missing_chr_arm < tmp_3mb[c,3]){
      if (c==1){
        tmp_3mb=rbind(tmp_0.1_3mb[which(tmp_0.1_3mb[,3] == missing_chr_arm),] , tmp_3mb[1:L_3mb,])
        c=L_3mb+1}
      else{
        tmp_3mb=rbind(tmp_3mb[1:(c-1),], tmp_0.1_3mb[which(tmp_0.1_3mb[,3] == missing_chr_arm),] , tmp_3mb[c:L_3mb,])}
      c=L_3mb+1}
    else{
      c=c+1
    }
  }
}

for (missing_chr_arm in values){
  tmp_0.1_3mb <- tmp_0.1_3mb[which(!(tmp_0.1_3mb[,3]==missing_chr_arm)),]}

L_3mb=dim(tmp_3mb)[1]
l_0.1_3mb=dim(tmp_0.1_3mb)[1]

i=1
c=1

while (i < l_0.1_3mb + 1){          # <=> Ligne arrêt non traitée && parcours segments entre 0.1 et 3 mb
  c=1                                                                            # reset parcours plus 3mb
  L_3mb=dim(tmp_3mb)[1]
  while (c < L_3mb+1){                                                             # parcour segments de plus de 3mb rassemblé
    #print(i, tmp_0.1_3mb[i,]) # debug pb ligne
    if (tmp_0.1_3mb[i,3] == tmp_3mb[c,3]){                                            # même chr_arm
      if (tmp_0.1_3mb[i,5] < tmp_3mb[c,4]){                         ############### 1 - Avant tout segs du chr_arm
        if (c==1){ # dès la première ligne, mais pas à relier
          if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) > THR){          # depassement pas de réunion 
            tmp_3mb=rbind(c(tmp_0.1_3mb[i,1],tmp_0.1_3mb[i,2],tmp_0.1_3mb[i,3],tmp_0.1_3mb[i,4],tmp_0.1_3mb[i,5], 
                            tmp_0.1_3mb[i,6],tmp_0.1_3mb[i,5]-tmp_0.1_3mb[i,4]+1, tmp_0.1_3mb[i,8]) , tmp_3mb[c:L_3mb,])
            i=i+1
            c=L_3mb+1}
          else{                     # réunion
            tmp_3mb=rbind(c(tmp_0.1_3mb[i,1], tmp_0.1_3mb[i,2], tmp_0.1_3mb[i,3], tmp_0.1_3mb[i,4], tmp_3mb[c,5], 
                            tmp_3mb[c,6],tmp_3mb[c,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c,8]), tmp_3mb[(c):L_3mb,])
            i=i+1
            c=L_3mb+1}
        }
        else{
          if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) > THR){                     # depassement THR, pas réunion
            tmp_3mb=rbind(tmp_3mb[1:(c-1),], c(tmp_0.1_3mb[i,1],tmp_0.1_3mb[i,2],tmp_0.1_3mb[i,3],tmp_0.1_3mb[i,4],tmp_0.1_3mb[i,5], 
                                               tmp_0.1_3mb[i,6],tmp_0.1_3mb[i,5]-tmp_0.1_3mb[i,4]+1, tmp_0.1_3mb[i,8]) , tmp_3mb[c:L_3mb,])
            i=i+1
            c=L_3mb+1}
          else{                                                                # reunion avec grand segment
            tmp_3mb=rbind(tmp_3mb[1:(c-1),], c(tmp_0.1_3mb[i,1], tmp_0.1_3mb[i,2], tmp_0.1_3mb[i,3], tmp_0.1_3mb[i,4], tmp_3mb[c,5], 
                                               tmp_3mb[c,6],tmp_3mb[c,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c,8]), tmp_3mb[(c+1):L_3mb,])
            i=i+1
            c=L_3mb+1}
        }
      }
      
      else if (tmp_0.1_3mb[i,4] >= tmp_3mb[c,4] && tmp_0.1_3mb[i,5] <= tmp_3mb[c,5]){    ########## 2 - A l'intérieur segs
        if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) > THR){
          if(c+1>L_3mb){          # A la fin                          
            tmp_3mb=rbind(tmp_3mb[1:(c-1),] , c(tmp_3mb[c,1],tmp_3mb[c,2] ,tmp_3mb[c,3], tmp_3mb[c,4], tmp_0.1_3mb[i,4]-1,tmp_3mb[c,6],tmp_0.1_3mb[i,4]-1-tmp_3mb[c,4]+1,tmp_3mb[c,8]) ,
                          tmp_0.1_3mb[i,] , c(tmp_3mb[c,1],tmp_3mb[c,2],tmp_3mb[c,3], tmp_0.1_3mb[i,5] +1 , tmp_3mb[c,5],tmp_3mb[c,6], tmp_3mb[c,5] - (tmp_0.1_3mb[i,5]+1) + 1,tmp_3mb[c,8]))
            i=i+1
            c=L_3mb+1}
          else{
            if (tmp_0.1_3mb[i,4] == tmp_3mb[c,4]){               # égalité début des deux segments
              tmp_3mb=rbind(tmp_3mb[1:(c-1),] , tmp_0.1_3mb[i,] , c(tmp_3mb[c,1],tmp_3mb[c,2],tmp_3mb[c,3], tmp_0.1_3mb[i,5] +1 , tmp_3mb[c,5],tmp_3mb[c,6], tmp_3mb[c,5] - (tmp_0.1_3mb[i,5]+1) + 1,tmp_3mb[c,8]) ,
                            tmp_3mb[(c+1):L_3mb,])
              i=i+1
              c=L_3mb+1}
            else if (tmp_0.1_3mb[i,5] == tmp_3mb[c,5]) {         # égalité fin des deux segments
              tmp_3mb=rbind(tmp_3mb[1:(c-1),] , c(tmp_3mb[c,1],tmp_3mb[c,2],tmp_3mb[c,3], tmp_3mb[i,4], tmp_0.1_3mb[c,4] - 1, tmp_3mb[c,6], tmp_0.1_3mb[c,4] - 1 - (tmp_3mb[i,4]) + 1,  tmp_3mb[c,8]), tmp_0.1_3mb[i,] , 
                            tmp_3mb[(c+1):L_3mb,])
              i=i+1
              c=L_3mb+1}
            else {               # classique au milieu
              tmp_3mb=rbind(tmp_3mb[1:(c-1),] , c(tmp_3mb[c,1], tmp_3mb[c,2], tmp_3mb[c,3], tmp_3mb[c,4], tmp_0.1_3mb[i,4]-1, tmp_3mb[c,6], tmp_0.1_3mb[i,4]-1-tmp_3mb[c,4]+1,tmp_3mb[c,8]) ,
                            tmp_0.1_3mb[i,] , c(tmp_3mb[c,1],tmp_3mb[c,2],tmp_3mb[c,3], tmp_0.1_3mb[i,5] +1 , tmp_3mb[c,5],tmp_3mb[c,6], tmp_3mb[c,5] - (tmp_0.1_3mb[i,5]+1) + 1,tmp_3mb[c,8]) ,
                            tmp_3mb[(c+1):L_3mb,])
              i=i+1
              c=L_3mb+1}
          }}
        else{                 # dans le fit, on ne fait rien
          i=i+1
          c=L_3mb+1}
      }
      
      else if (tmp_0.1_3mb[i,4] >= tmp_3mb[c,5] && tmp_0.1_3mb[i,5] <= tmp_3mb[c+1,4] && tmp_3mb[c,3] == tmp_3mb[c+1,3]){     ########## 3 - entre deux segs de même chr_arm
        if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) <= THR && abs(tmp_0.1_3mb[i,6] - tmp_3mb[c+1,6]) <= THR){       # deux in_range THR
          if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) <= abs(tmp_0.1_3mb[i,6] - tmp_3mb[c+1,6])){     # c plus proche
            tmp_3mb=rbind(tmp_3mb[1:(c-1),] , c(tmp_3mb[c,1], tmp_3mb[c,2], tmp_3mb[c,3], tmp_3mb[c,4], tmp_0.1_3mb[i,5], tmp_3mb[c,6],
                                                tmp_0.1_3mb[i,5]-tmp_3mb[c,4]+1,tmp_3mb[c,8]) , tmp_3mb[(c+1):L_3mb,])  
            i=i+1
            c=L_3mb+1}
          else {                      # c+1 plus proche
            if (c+2 > L_3mb){
              tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_3mb[c+1,1], tmp_3mb[c+1,2], tmp_3mb[c+1,3], tmp_0.1_3mb[i,4], tmp_3mb[c+1,5], tmp_3mb[c+1,6],
                                              tmp_3mb[c+1,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c+1,8])) 
              i=i+1
              c=L_3mb+1}
            else{
              tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_3mb[c+1,1], tmp_3mb[c+1,2], tmp_3mb[c+1,3], tmp_0.1_3mb[i,4], tmp_3mb[c+1,5], tmp_3mb[c+1,6],
                                              tmp_3mb[c+1,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c+1,8]) , tmp_3mb[(c+2):L_3mb,])  
              i=i+1
              c=L_3mb+1}
          }
        }
        else if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) <= THR){                   # seulement c
          tmp_3mb=rbind(tmp_3mb[1:(c-1),] , c(tmp_3mb[c,1], tmp_3mb[c,2], tmp_3mb[c,3], tmp_3mb[c,4], tmp_0.1_3mb[i,5], tmp_3mb[c,6],
                                              tmp_0.1_3mb[i,5]-tmp_3mb[c,4]+1,tmp_3mb[c,8]) , tmp_3mb[(c+1):L_3mb,])  
          i=i+1
          c=L_3mb+1}
        else if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c+1,6]) <= THR){                  # seulement c+1
          if (c+2 > L_3mb){
            tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_3mb[c+1,1], tmp_3mb[c+1,2], tmp_3mb[c+1,3], tmp_0.1_3mb[i,4], tmp_3mb[c+1,5], tmp_3mb[c+1,6],
                                            tmp_3mb[c+1,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c+1,8]))  
            i=i+1
            c=L_3mb+1}
          else{
            if (tmp_0.1_3mb[i+1,5] < tmp_3mb[c+1,4] && tmp_0.1_3mb[i+1,3] == tmp_3mb[c+1,3]){  # il reste des petits segments à caller après
              tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_0.1_3mb[i,1], tmp_0.1_3mb[i,2], tmp_0.1_3mb[i,3], tmp_0.1_3mb[i,4], tmp_0.1_3mb[i,5], tmp_0.1_3mb[i,6],
                                              tmp_0.1_3mb[i,5]-tmp_0.1_3mb[i,4]+1,tmp_0.1_3mb[i,8]) , tmp_3mb[(c+1):L_3mb,])  
              i=i+1
              c=L_3mb+1}
            else{
              tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_3mb[c+1,1], tmp_3mb[c+1,2], tmp_3mb[c+1,3], tmp_0.1_3mb[i,4], tmp_3mb[c+1,5], tmp_3mb[c+1,6],
                                              tmp_3mb[c+1,5]-tmp_0.1_3mb[i,4]+1,tmp_3mb[c+1,8]) , tmp_3mb[(c+2):L_3mb,])  
              i=i+1
              c=L_3mb+1}
          }
        }
        else{                                      
          tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_0.1_3mb[i,1],tmp_0.1_3mb[i,2],tmp_0.1_3mb[i,3],tmp_0.1_3mb[i,4],
                                          tmp_0.1_3mb[i,5],tmp_0.1_3mb[i,6],tmp_0.1_3mb[i,7],tmp_0.1_3mb[i,8]) , tmp_3mb[(c+1):L_3mb,])
          i=i+1
          c=L_3mb+1}
      }
      else if (tmp_0.1_3mb[i,4] >= tmp_3mb[c,5] && (tmp_3mb[c,3] != tmp_3mb[c+1,3] | is.null(tmp_3mb[c+1,3]))){     ############## 4 - après dernier segs chr_arm
        if (abs(tmp_0.1_3mb[i,6] - tmp_3mb[c,6]) <= THR){
          if (c+1 > L_3mb){
            tmp_3mb=rbind(tmp_3mb[1:(c-1),] , c(tmp_3mb[c,1],tmp_3mb[c,2],tmp_3mb[c,3],tmp_3mb[c,4],
                                                tmp_0.1_3mb[i,5],tmp_3mb[c,6],tmp_0.1_3mb[i,5]-tmp_3mb[c,4]+1,tmp_3mb[c,8])) 
            i=i+1
            c=L_3mb+1}
          else{
            tmp_3mb=rbind(tmp_3mb[1:(c-1),] , c(tmp_3mb[c,1],tmp_3mb[c,2],tmp_3mb[c,3],tmp_3mb[c,4],
                                                tmp_0.1_3mb[i,5],tmp_3mb[c,6],tmp_0.1_3mb[i,5]-tmp_3mb[c,4]+1,tmp_3mb[c,8]) , tmp_3mb[(c+1):L_3mb,])
            i=i+1
            c=L_3mb+1}}
        else{
          if (c+1 > L_3mb){
            tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_0.1_3mb[i,1],tmp_0.1_3mb[i,2],tmp_0.1_3mb[i,3],tmp_0.1_3mb[i,4],
                                            tmp_0.1_3mb[i,5],tmp_0.1_3mb[i,6],tmp_0.1_3mb[i,7],tmp_0.1_3mb[i,8]))
            i=i+1
            c=L_3mb+1}
          else{
            tmp_3mb=rbind(tmp_3mb[1:c,] , c(tmp_0.1_3mb[i,1],tmp_0.1_3mb[i,2],tmp_0.1_3mb[i,3],tmp_0.1_3mb[i,4],
                                            tmp_0.1_3mb[i,5],tmp_0.1_3mb[i,6],tmp_0.1_3mb[i,7],tmp_0.1_3mb[i,8]) , tmp_3mb[(c+1):L_3mb,])
            i=i+1
            c=L_3mb+1}}}
      
      else {
        c=c+1
      }
    }    
    else{                                            
      c=c+1
    }
  }
}

options(show.error.messages = TRUE)

write.table(tmp_3mb, file = paste(NAMEEE,"_",short_size_window,"kb_IV", sep=""), sep = "\t", row.names = FALSE)

tmp_3mb<-breakSmoothToLST(THR,tmp_3mb,c_ind=c_ind,c_chr=c_chr,c_posS=c_posS,c_posE=c_posE,c_cn=c_cn,c_conf=c_conf)

write.table(tmp_3mb, file = paste(NAMEEE,"_",short_size_window,"kb_V", sep=""), sep = "\t", row.names = FALSE)


##### Graphe final segmentation diag #####

test_data_frame = as.data.frame(tmp_3mb)
colnames(test_data_frame) <- make.unique(names(test_data_frame))

test_ordered = test_data_frame[order(test_data_frame$ratio_median),]
test_ordered$num_line <- seq.int(nrow(test_ordered))


test_ploty_CN_level <- ggplot(test_ordered, aes(x = ratio_median, y = num_line)) + 
  geom_point(shape = 1, size = 1, color = "#0072B2", na.rm = TRUE) + geom_line() + ggtitle("Final segmentation diagnostic") + 
  xlab("ratio median") + xlim(-2.5,2.5) + 
  theme(plot.title = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.y = element_blank(),
        panel.background = element_blank())
suppressWarnings(ggsave(paste(NAMEEE,"_test_ploty_CN_level", sep = ""), plot = test_ploty_CN_level, device = "jpeg", width = 20, height = 10))


##### Call LSTs #####

LSTs_data_frame <- as.data.frame(matrix(0, ncol = 2, nrow = 9))
LSTs_data_frame[,1] = c(3:11)
colnames(LSTs_data_frame) <- c("Size_LST", "Number_LST")

for (i in (3:11)){
  WC<-LST_control(THR,lenBIN=500,lenMB=i,tmp_3mb,c_ind=c_ind,c_chr=c_chr,c_posS=c_posS,c_posE=c_posE,c_cn=c_cn,c_conf=c_conf)
  LSTs_data_frame[i-2,2] = sum(WC[,1])}
write.table(LSTs_data_frame, file = paste(NAMEEE,"_LSTs", sep=""), sep = "\t", row.names = FALSE)

# pour graphe 10Mb 

WC<-LST_control(THR,lenBIN=500,lenMB=10,tmp_3mb,c_ind=c_ind,c_chr=c_chr,c_posS=c_posS,c_posE=c_posE,c_cn=c_cn,c_conf=c_conf)


test=cbind(tmp_3mb,WC[,1])
colnames(test) <- c("index", "chr","chr_arm", "start", "end", "ratio_median", "size", "level", "WC")

L=dim(test)[1]

i=1
while(i + 1 < L){
  L=dim(test)[1]
  if (test[i,9] == 0 & test[i+1,9] != 1){
    test=test[-i,]}
  else{
    i=i+1}}

l=dim(test)[1]
if (is.null(l) == FALSE){
  test=test[-l,]
}

write.table(test, file = paste(NAMEEE,"_",short_size_window,"kb_VI", sep=""), sep = "\t", row.names = FALSE)



##### Graphe different steps LSTs calling procedure #####

dataTable <- read.table(paste("Ratio_",NAMEEE,"_",short_size_window,"kb.tsv", sep = ""), header=TRUE)
B <- data.frame(dataTable)
B=B[,-1]
B=B[,-5]
colnames(B) <- c("chr", "start", "end", "readcount")
attach(B)

df=merge(aggregate(start ~ chr, B, min), aggregate(end ~ chr, B, max))[seq(from=1, to=22, by = 2),]
df=as.data.frame(df)
detach(B)


## DEBUT : Normalised and corrected readcount filtered 

dataTable <- read.table(paste("Ratio_",NAMEEE,"_",short_size_window,"kb.tsv", sep = ""), header=TRUE)
B <- data.frame(dataTable)
B=B[,-1] # remove feature
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

Z <- ggplot() +
  geom_rect(data=df, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), fill = "grey80") +
  geom_point(data=B, aes(x = start, y = ratio), size=0.1, na.rm=TRUE) + ylim(-2, 2) +
  scale_x_continuous(expand = c(0, 0)) + xlab("chromosomes") +
  facet_grid(~chr, scales = "free_x", space = "free_x", switch = "x")

Z <- Z + theme(axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.text.y = element_blank(),
               axis.title.y = element_blank(),
               # axis.text.y = element_blank(),
               # axis.ticks.y = element_blank(),
               panel.spacing = unit(0, "lines"),
               strip.text.x = element_blank(),
               line = element_blank(),
               panel.background = element_blank())

suppressWarnings(ggsave(paste("Normalised_Read_Count",NAMEEE,sep=""), plot = Z, device = "jpeg", width = 23, height = 13))

## FIN : Normalised and corrected readcount filtered


##   DEBUT : LSTs pathway

dataTable <- read.table(paste("Ratio_",NAMEEE,"_",short_size_window,"kb.tsv", sep = ""), header=TRUE)
B <- data.frame(dataTable)
B=B[,-1] # remove feature
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

dataTable <- read.table(paste(NAMEEE,"_",short_size_window,"kb_I", sep=""), header=TRUE)
C <- data.frame(dataTable)

I <- ggplot() +
  geom_rect(data=df, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), fill = "grey85") +
  geom_point(data=B, aes(x = start, y = ratio), size=0.1, color= "grey60") +
  geom_segment(data=C, aes(x=start, xend=end, y=ratio_median, yend=ratio_median), size = 3, color = "#990000") +
  ylim(-2, 2) +   scale_x_continuous(expand = c(0, 0)) + xlab("chromosomes") +
  facet_grid(~chr, scales = "free_x", space = "free_x", switch = "x")

I <- I + theme(axis.title.x = element_text(size=20),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_text(size=20),
               # axis.text.y = element_blank(),
               # axis.ticks.y = element_blank(),
               panel.spacing = unit(0, "lines"),
               strip.text.x = element_blank(),
               line = element_blank(),
               panel.background = element_blank())

suppressWarnings(ggsave(paste("1_",NAMEEE,sep=""), plot = I, device = "jpeg", width = 23, height = 13))


dataTable <- read.table(paste("Ratio_",NAMEEE,"_",short_size_window,"kb.tsv", sep = ""), header=TRUE)
B <- data.frame(dataTable)
B=B[,-1] # remove feature
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

dataTable <- read.table(paste(NAMEEE,"_",short_size_window,"kb_II", sep=""), header=TRUE)
C <- data.frame(dataTable)

II <- ggplot() +
  geom_rect(data=df, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), fill = "grey85") +
  geom_point(data=B, aes(x = start, y = ratio), size=0.1, color= "grey60") +
  geom_segment(data=C, aes(x=start, xend=end, y=ratio_median, yend=ratio_median), size = 3, color = "#990000") +
  ylim(-2, 2) +   scale_x_continuous(expand = c(0, 0)) + xlab("chromosomes") +
  facet_grid(~chr, scales = "free_x", space = "free_x", switch = "x")

II <- II + theme(axis.title.x = element_text(size=20),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.y = element_text(size=20),
                 # axis.text.y = element_blank(),
                 # axis.ticks.y = element_blank(),
                 panel.spacing = unit(0, "lines"),
                 strip.text.x = element_blank(),
                 line = element_blank(),
                 panel.background = element_blank())

suppressWarnings(ggsave(paste("2_",NAMEEE,sep=""), plot = II, device = "jpeg", width = 23, height = 13))


dataTable <- read.table(paste("Ratio_",NAMEEE,"_",short_size_window,"kb.tsv", sep = ""), header=TRUE)
B <- data.frame(dataTable)
B=B[,-1] # remove feature
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

dataTable <- read.table(paste(NAMEEE,"_",short_size_window,"kb_III", sep=""), header=TRUE)
C <- data.frame(dataTable)

III <- ggplot() +
  geom_rect(data=df, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), fill = "grey85") +
  geom_point(data=B, aes(x = start, y = ratio), size=0.1, color= "grey60") +
  geom_segment(data=C, aes(x=start, xend=end, y=ratio_median, yend=ratio_median), size = 3, color = "#990000") +
  ylim(-2, 2) +   scale_x_continuous(expand = c(0, 0)) + xlab("chromosomes") +
  facet_grid(~chr, scales = "free_x", space = "free_x", switch = "x")

III <- III + theme(axis.title.x = element_text(size=20),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.title.y = element_blank(),
                   axis.text.y = element_text(size=20),
                   # axis.text.y = element_blank(),
                   # axis.ticks.y = element_blank(),
                   panel.spacing = unit(0, "lines"),
                   strip.text.x = element_blank(),
                   line = element_blank(),
                   panel.background = element_blank())

suppressWarnings(ggsave(paste("3_",NAMEEE,sep=""), plot = III, device = "jpeg", width = 23, height = 13))



dataTable <- read.table(paste("Ratio_",NAMEEE,"_",short_size_window,"kb.tsv", sep = ""), header=TRUE)
B <- data.frame(dataTable)
B=B[,-1] # remove feature
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

dataTable <- read.table(paste(NAMEEE,"_",short_size_window,"kb_IV", sep=""), header=TRUE)
C <- data.frame(dataTable)

IV <- ggplot() +
  geom_rect(data=df, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), fill = "grey85") +
  geom_point(data=B, aes(x = start, y = ratio), size=0.1, color= "grey60") +
  geom_segment(data=C, aes(x=start, xend=end, y=ratio_median, yend=ratio_median), size = 3, color = "#990000") +
  ylim(-2, 2) +   scale_x_continuous(expand = c(0, 0)) + xlab("chromosomes") +
  facet_grid(~chr, scales = "free_x", space = "free_x", switch = "x")

IV <- IV + theme(axis.title.x = element_text(size=20),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.y = element_text(size=20),
                 # axis.text.y = element_blank(),
                 # axis.ticks.y = element_blank(),
                 panel.spacing = unit(0, "lines"),
                 strip.text.x = element_blank(),
                 line = element_blank(),
                 panel.background = element_blank())

suppressWarnings(ggsave(paste("4_",NAMEEE,sep=""), plot = IV, device = "jpeg", width = 23, height = 13))


dataTable <- read.table(paste("Ratio_",NAMEEE,"_",short_size_window,"kb.tsv", sep = ""), header=TRUE)
B <- data.frame(dataTable)
B=B[,-1] # remove feature
colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

dataTable <- read.table(paste(NAMEEE,"_",short_size_window,"kb_V", sep=""), header=TRUE)
C <- data.frame(dataTable)


closest_higlight = Closest(B$start, 30302805)[1]

higlight_CCNE1 = B[B$chr == 19 & B$start == closest_higlight,]

data.segm = data.frame(x=higlight_CCNE1[1,2], y=higlight_CCNE1[1,5] + 1.4, xend=higlight_CCNE1[1,2], yend=higlight_CCNE1[1,5] + 0.05, chr = 19)

data.text = data.frame(x=higlight_CCNE1[1,2], y=higlight_CCNE1[1,5] + 1.5, chr = 19, label = "CCNE1")


V <- ggplot() +
  geom_rect(data=df, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), fill = "grey85") +
  geom_point(data=B, aes(x = start, y = ratio), size=0.1, color= "grey60", na.rm=TRUE) +
  geom_segment(data=C, aes(x=start, xend=end, y=ratio_median, yend=ratio_median), size = 3, color = "#990000", na.rm = TRUE) +
  geom_segment(data=data.segm, mapping=aes(x=x, y=y, xend=xend, yend=yend), 
               arrow=arrow(unit(0.30,"cm"), angle = 20), size=0.6, color="black", inherit.aes = FALSE, na.rm=TRUE) +
  ylim(-2, 2) +   scale_x_continuous(expand = c(0, 0)) + ggtitle(NAMEEE) +
  facet_grid(~chr, scales = "free_x", space = "free_x", switch = "x") + 
  geom_text(data = data.text, mapping = aes(x = x, y = y, label = label), size = 3.3, inherit.aes = TRUE, na.rm=TRUE) +
  geom_point(data = higlight_CCNE1, aes(x=start, y=ratio_median), color = "orange", size = 2, na.rm=TRUE)


V <- V + theme(plot.title = element_text(hjust = 0.5),
               axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.ticks.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.y = element_text(size=15),
               # axis.text.y = element_blank(),
               # axis.ticks.y = element_blank(),
               # text = element_text(size = 70),
               panel.spacing = unit(0, "lines"),
               strip.text.x = element_blank(),
               line = element_blank(),
               panel.background = element_blank())

V = ggplotGrob(x = V)
V$layout$clip = "off"

suppressWarnings(ggsave(paste("5_",NAMEEE,sep=""), plot = V, device = "jpeg", width = 23, height = 13))

##  FIN : LSTs pathway


##    DEBUT : LSTs called
if (sum(WC[,1]) != 0){

  dataTable <- read.table(paste("Ratio_",NAMEEE,"_",short_size_window,"kb.tsv", sep = ""), header=TRUE)
  B <- data.frame(dataTable)
  B=B[,-1] # remove feature
  colnames(B) <- c("chr", "start", "end", "ratio", "ratio_median")

  dataTable <- read.table(paste(NAMEEE,"_",short_size_window,"kb_VI", sep=""), header=TRUE)
  C <- data.frame(dataTable)

  closest_higlight = Closest(B$start, 30302805)[1]
  
  higlight_CCNE1 = B[B$chr == 19 & B$start == closest_higlight,]

  data.segm = data.frame(x=higlight_CCNE1[1,2], y=higlight_CCNE1[1,5] + 1.4, xend=higlight_CCNE1[1,2], yend=higlight_CCNE1[1,5] + 0.05, chr = 19)

  data.text = data.frame(x=higlight_CCNE1[1,2], y=higlight_CCNE1[1,5] + 1.5, chr = 19, label = "CCNE1")

  VI <- ggplot() +
    geom_rect(data=df, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf), fill = "grey85") +
    geom_point(data=B, aes(x = start, y = ratio), size=0.1, color= "grey60", na.rm=TRUE) +
    geom_segment(data=C, aes(x=start, xend=end, y=ratio_median, yend=ratio_median), size = 3, color = "#006600", na.rm=TRUE) +
    geom_segment(data=data.segm, mapping=aes(x=x, y=y, xend=xend, yend=yend), 
               arrow=arrow(unit(0.30,"cm"), angle = 20), size=0.6, color="black", inherit.aes = FALSE, na.rm=TRUE) +
    ylim(-2, 2) +   scale_x_continuous(expand = c(0, 0)) + ggtitle(NAMEEE) +
    facet_grid(~chr, scales = "free_x", space = "free_x", switch = "x") +
    geom_point(data = higlight_CCNE1, aes(x=start, y=ratio_median), color = "orange", size = 2, na.rm=TRUE) +
    geom_text(data = data.text, mapping = aes(x = x, y = y, label = label), size = 4, inherit.aes = TRUE, na.rm=TRUE)


  VI <- VI + theme(plot.title = element_text(hjust = 0.5),
                 axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.y = element_text(size=15),
                 # axis.text.y = element_blank(),
                 # axis.ticks.y = element_blank(),
                 panel.spacing = unit(0, "lines"),
                 strip.text.x = element_blank(),
                 line = element_blank(),
                 panel.background = element_blank())

  VI = ggplotGrob(x = VI)
  VI$layout$clip = "off"

  suppressWarnings(ggsave(paste("6_",NAMEEE,sep=""), plot = VI, device = "jpeg", width = 23, height = 13))
}

##### Output final figure #####
## graphe

graphe = V
if (sum(WC[,1]) != 0){
  graphe = VI
}

## ploidy controlFREEC 

ploidy <- read.table(paste(NAMEEE,".bam_info.txt", sep = ""), header=TRUE)

ploidy_control = ploidy[10,][2]
colnames(ploidy_control) <- NULL
rownames(ploidy_control) <- NULL
ploidy_test = as.integer(as.character(ploidy_control[[1]]))

## Threhold

Threshold = round(Threshold, 3)

## THR (corrected)

THR = round(THR, 3)

## number LSTs

number_LSTs = sum(WC[,1])

## CCNE1 amplification

CN_plus_baseline = round(higlight_CCNE1$ratio_median/THR,3)

## CCNE1 amplification diagnostic

if (CN_plus_baseline >= 4){
  CCNE1_diag = "Amplification : not BRCAness ? (>= 4)"
} else if (CN_plus_baseline >= 3){
  CCNE1_diag = "Gain [3,4["
} else{
  CCNE1_diag = "Not informative (<3)"
}

## BRCAness

if (number_LSTs >= 18){
  BRCAness = "Yes (>= 18)"
} else if (number_LSTs >= 16){
  BRCAness = "Borderline (16/17)"
} else{
  BRCAness = "No (<= 15)"
}

## Table

a = c("Ploidy ControlFREEC", "Threshold value","Threshold corrected","Number LSTs 10Mb", "CCNE1 CN to baseline", "CCNE1 evidence", "BRCAness")
b = c(ploidy_test, Threshold, THR, number_LSTs, CN_plus_baseline, CCNE1_diag, BRCAness)
A = data.frame(a,b)

colnames(A) <- NULL
row.names(A) <- NULL

testy <- qplot(1:10, 1:10, geom = "blank") + 
  theme(line = element_blank(), text = element_blank(), panel.background = element_blank()) +
  annotation_custom(grob = tableGrob(A, rows = NULL))

## Summary figure

summary_plot <- suppressWarnings(ggarrange(graphe,
                          ggarrange(ploty, test_ploty_CN_level, testy, labels = c("B", "C", "D"), ncol = 3, font.label = list(size = 14, face = "bold")),
                          nrow = 2, labels = "A", font.label = list(size = 14, face = "bold")))

suppressWarnings(ggsave(paste("summary_plot_",NAMEEE,sep=""), plot = summary_plot, device = "jpeg", width = 18, height = 10.17))

print("========================================================")
print("========================================================")
print("========================================================")
print(paste("Sample : ", NAMEEE, sep = ""))
print(paste("Number of LSTs (> 10Mb) : ", number_LSTs, sep = ""))
print(paste("BRCAness : ", BRCAness, sep = ""))