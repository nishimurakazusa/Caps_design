
args <- commandArgs(trailingOnly = T)

default_args <- c("2000", "3")
default_flg <- is.na(args[1:2])
args[default_flg] <- default_args[default_flg]  

#set parameters
extract_seq_length <- as.numeric(args[1])

options(scipen=100)

dir.create("results")

df_res4 = read.csv(file="caps_list.csv",stringsAsFactors = F)

###########Creating gtf file

seqname   = df_res4[,2]
source_g  = rep("test",nrow(df_res4))
feature   = rep("CDS",nrow(df_res4))
start     = (as.numeric(df_res4[,3]) - (extract_seq_length/2))
end       = (as.numeric(df_res4[,3]) + (extract_seq_length/2))
score     = rep(".",nrow(df_res4))
strand    = rep("+",nrow(df_res4))
frame     = rep(".",nrow(df_res4))
attribute = row.names(df_res4)

start[start<0] = 0

gtf = data.frame(seqname,source_g,feature,start,end,score,strand,frame,attribute,stringsAsFactors = F)  
  
write.table(gtf,file="results/caps_primer.gtf",sep="\t",col.names=F,row.names=F,quote=F)

message(paste("you created gtf file"))


