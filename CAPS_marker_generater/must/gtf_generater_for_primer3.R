
args <- commandArgs(trailingOnly = T)

# Unspecified options set default values
default_args <- c("900", "3")   #Default value definition
default_flg <- is.na(args[1:2])
args[default_flg] <- default_args[default_flg]  

#Set file name and parameters
extract_seq_length <- as.numeric(args[1])

options(scipen=100)

df_res4 = read.csv(file="results/candidate_of_snp_amd_enzyme_list.csv",stringsAsFactors = F,row.names = 1)

###########Generate gtf file

seqname   = df_res4[,1]
source_g  = rep("test",nrow(df_res4))
feature   = rep("CDS",nrow(df_res4))
start     = (as.numeric(df_res4[,2]) - (extract_seq_length/2))
end       = (as.numeric(df_res4[,2]) + (extract_seq_length/2))
score     = rep(".",nrow(df_res4))
strand    = rep("+",nrow(df_res4))
frame     = rep(".",nrow(df_res4))
attribute = row.names(df_res4)


start[start<0] = 0


gtf = data.frame(seqname,source_g,feature,start,end,score,strand,frame,attribute,stringsAsFactors = F)  
  
write.table(gtf,file="results/caps_primer.gtf",sep="\t",col.names=F,row.names=F,quote=F)

message(paste("you created gtf file"))


