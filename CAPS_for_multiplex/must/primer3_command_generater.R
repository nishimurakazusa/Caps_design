#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.10")
#BiocManager::install(c("Biostrings"))

require(Biostrings)

args <- commandArgs(trailingOnly = T)

#Unspecified options set default values
default_args <- c("100", "550","2000")   # Default value definition 
default_flg <- is.na(args[1:3])
args[default_flg] <- default_args[default_flg]  

#Set file name, parameters
amprified_length1 <- as.numeric(args[1])
amprified_length2 <- as.numeric(args[2])
ex_length <- as.numeric(args[3])

filename = "results/seq_for_caps_primer.fasta"

fasta = readDNAStringSet(filename)

###This is very important!
fasta_name = fasta@ranges@NAMES
df_fasta_for_p= as.data.frame(fasta) 

##To sort here in the style of the original file.
d1 = as.data.frame(strsplit(fasta_name, "_"))
d1 = t(d1)
head(d1)

d2 = as.data.frame(strsplit(d1[,ncol(d1)],"-"))
d2 = t(d2)
head(d2)

d3 = as.data.frame(strsplit(d2[,2],":"))
d3 = t(d3)

f2 = data.frame(d1[,1],as.numeric(d2[,1]),as.numeric(d3[,1]),fasta_name)
colnames(f2) = c("chrom","start","end","fasta_name_for_primer3")

##Integration once here
df_fasta_for_p = cbind(df_fasta_for_p,f2[,1:4])

df_fasta_for_p = df_fasta_for_p[order(df_fasta_for_p$chrom),]

dir.create("df")

write.csv(file="df/df_fasta_for_p.csv",df_fasta_for_p)
df2_fasta_for_p = read.csv(file="df/df_fasta_for_p.csv",header=T,stringsAsFactors = F)

dir.create("primer3")
dir.create("primer3/input")
dir.create("primer3/output")

primer3_example = read.csv(file="must/primer3_v2.csv",stringsAsFactors = F,header = F)

primer3_example[10,1] =　paste0("PRIMER_PRODUCT_SIZE_RANGE=",amprified_length1,"-",amprified_length2)

primer3_example[10,1] =　paste0("SEQUENCE_TARGET=",ex_length/2 -15,",30")


#i=1
for(i in 1:nrow(df2_fasta_for_p)){
  
  for_massage = paste("File No. ",i,sep="")
  message(for_massage)
    primer3 = matrix("NA",nrow = nrow(primer3_example), ncol = 1)
    primer3[1,1] = paste("SEQUENCE_ID=",df2_fasta_for_p[i,1],sep="")
    primer3[2,1] = paste("SEQUENCE_TEMPLATE=",df2_fasta_for_p[i,2],sep="")
    primer3[3:nrow(primer3_example),1] = primer3_example[3:nrow(primer3_example),1]
  
    cat(primer3,"\n",file = paste("primer3/input/input",i,sep=""), append = F,sep = "\n") 
}

message(paste("you created",nrow(df2_fasta_for_p),"files for primer3 !!"))

