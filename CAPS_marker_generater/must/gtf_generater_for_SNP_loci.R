
args <- commandArgs(trailingOnly = T)

#Unspecified options set default values
default_args <- c("filtered.citrus_rd_29_minDP10_0.5_samtools.vcf.recode.vcf", "3")   #Default value definition
default_flg <- is.na(args[1:2])
args[default_flg] <- default_args[default_flg]  

# Set file name and parameters
file_name <- args[1]
p <- as.numeric(args[2])

require(vcfR)

dir.create("results")

options(scipen=100)

# read compressed vcf file (*.vcf.gz)
vcf <- read.vcfR(paste(file_name))

# extract genotype data from vcf
gt = extract.gt(vcf)
#gt[1:10,1:2]
#dim(gt)

# extract depyh data from vcf
#gq <- extract.gt(vcf,element = 'GQ', as.numeric = TRUE)
#gq[1:10,1:2]

dp <- extract.gt(vcf,element = 'DP', as.numeric = TRUE)

# extract depyh data from vcf
ref <- getREF(vcf)
#ref[1:10]
alt <- getALT(vcf)
#alt[1:10]

# get marker information (chromosone numbers and positions)
chrom <- getCHROM(vcf)
pos <- getPOS(vcf)

#show the first 10 SNPs of the first 10 lines 
#gt[1:10, 1:2]

# create a matrix of gt scores
gt.score <- matrix(NA, nrow(gt), ncol(gt))
gt.score[gt == "0/0"] <- -1
gt.score[gt == "0/1"] <- 0
gt.score[gt == "1/1"] <- 1
gt.score[gt == "0|0"] <- -1
gt.score[gt == "0|1"] <- 0
gt.score[gt == "1|1"] <- 1

#gt.score[1:10, 1:2]


# name the rows and columns of matrix
rownames(gt.score) <- rownames(gt)
colnames(gt.score) <- colnames(gt)
#gt.score[1:10, 1:2]


#filtering by depth
#set parameter

FilterbyDP = function(a){
  
  filter_dp = a
  gt.score[dp < filter_dp] = NA
  
  all_data = data.frame(chrom,pos,ref,alt,gt.score)
  #all_data[1:10,1:6]
  
  return(all_data)
}

#############Filtered by DP.
all_data = FilterbyDP(p)
#all_data[1:10,1:6]
#dim(all_data)
#str(all_data)

#extracted hoge% typing loci
df = all_data
#i=1

filter_na_rate = function(ex_df){
 1- sum (as.numeric( is.na (  ex_df   ) )) / ( ncol(df) -4 )
}

na.omit_var = function(x){
  y = var(na.omit(x))
  return(y)
  
}

df2 =  df[ apply(df[,5:ncol(df)],1,filter_na_rate) >= 0.5,  ]  ###Specify missing data percentage

df2 = df2[apply(df2[,5:ncol(df)],1,na.omit_var) != 0,]


#No heterogeneity allowed in the case of two systems.
#df2 = df[!is.na(df[,3]) & !is.na(df[,4]) & df[,3] != 0 & df[,4] != 0 & df[,3] != df[,4], ]
#Allow heterogeneity in the case of two systems
#df2 = df[!is.na(df[,3]) & !is.na(df[,4])  & df[,3] != df[,4], ]

#df2[1:10,1:6]
#str(df2)
df3 = df2[order(as.numeric(df2$pos)),]
df3 = df3[order(as.numeric(as.factor(df3$chrom))),]

#dir.create("df")
#write.csv(df2,file="df/df2.csv")

dir.create("df")
write.csv(df3,file="df/df3.csv")
write.csv(df3,file="results/genotype_data.csv")

###########Generate gtf file for CAPS search

seqname   = df2[,1]
source_g  = rep("test",nrow(df2))
feature   = rep("CDS",nrow(df2))
start     = (df2[,2] - 5)
end       = (df2[,2] + 5)
score     = rep(".",nrow(df2))
strand    = rep("+",nrow(df2))
frame     = rep(".",nrow(df2))
attribute = row.names(df2)

gtf = data.frame(seqname,source_g,feature,start,end,score,strand,frame,attribute)  

write.table(gtf,file="results/caps.gtf",sep="\t",col.names=F,row.names=F,quote=F)

