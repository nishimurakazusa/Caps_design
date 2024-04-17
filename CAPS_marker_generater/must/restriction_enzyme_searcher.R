
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.10")
#BiocManager::install(c("Biostrings"))
#BiocManager::install(c("S4Vectors"))
require(Biostrings)

args <- commandArgs(trailingOnly = T)

#Unspecified options set default values
default_args <- c("results/caps.fasta", "3")   #Default value definition
default_flg <- is.na(args[1:2])
args[default_flg] <- default_args[default_flg]  


#Set file name and parameters
file_name <- args[1]

##
f1 = readDNAStringSet(file_name, use.names = T)

names = as.data.frame(names(f1))
df_fasta_seq= as.data.frame(f1)
df_fasta= cbind(names,df_fasta_seq[,1])
#df_fasta[,2] = as.character(df_fasta[,2])
#df_fasta[,2] = as.character(df_fasta[,2])

write.csv(file="df/df_fasta.csv",df_fasta)
df2_fasta = read.csv(file="df/df_fasta.csv",header=T,stringsAsFactors = F,row.names = 1)

Unzip <- function(...) rbind(data.frame(), ...)

d1 = strsplit(df2_fasta[,1], "_")
d1.frame <- do.call(Unzip, d1)

#str(d1.frame)

d2 = strsplit(as.character(d1.frame[,ncol(d1.frame)]),"-")
d2.frame <- do.call(Unzip, d2)

d3 = strsplit(as.character(d2.frame[,2]),":")
d3.frame <- do.call(Unzip, d3)

f2 = data.frame(d1.frame[,1],as.character(d2.frame[,1]),as.character(d3.frame[,1]),df2_fasta[,2])
f2[,2] = as.character(f2[,2])
f2 = cbind(f2,(as.numeric(f2[,2])+5) )
f2[,3] = as.character(f2[,3])
f2[,4] = as.character(f2[,4])
#str(f2)
colnames(f2) = c("chrom","start","end","string","pos")

#f2[,1]=as.numeric(f2[,1])
#str(f2)

f3 = f2[order(as.numeric(f2$pos)),]
f4 = f3[order(as.numeric(as.factor(f3$chrom)) ),]
#head(f4)

df3 = read.csv(file="df/df3.csv",header = T,stringsAsFactors = F)
df_conb = cbind(df3,f4)
#df_conb = cbind(df2,f2)
#head(df_conb)

write.csv(df_conb, file= "df/df_conb.csv")
#error = df_conb[df_conb[,3]-df_conb[12] != 0,]
#error[,3] - error[,12]


ATGC = (df_conb$ref == "A"|
          df_conb$ref == "T"|
          df_conb$ref == "G"|
          df_conb$ref == "C")&
       (df_conb$alt == "A"|
          df_conb$alt == "T"|
          df_conb$alt == "G"|
          df_conb$alt == "C")

near_to_indel = rep(NA,nrow(df_conb))
multiple_snp  = rep(NA,nrow(df_conb))


df_conb = cbind(df_conb,as.numeric(ATGC),near_to_indel,multiple_snp)


##Mark SNPs near indel
#i=4
#for(i in c(1:nrow(df_conb))[ATGC==F]){
#  message(paste0(i))  
#  df_conb$near_to_indel[ df_conb$pos[i]-5 < df_conb$pos &
#                         df_conb$pos[i]+5 > df_conb$pos &
#                         df_conb$chrom == df_conb$chrom[i]] = 2
#}

df_conb = df_conb[ATGC,]


##Create mutant allele file from here
ref_seq11_fasta = as.character(df_conb$string)
ref_seq7_fasta = paste(substr(ref_seq11_fasta, 3, 9))
extracted_alt = df_conb$alt

alt_seq11_fasta = paste(substr(ref_seq11_fasta, 1, 5),extracted_alt,substr(ref_seq11_fasta, 7, 11),sep="")
alt_seq7_fasta = paste(substr(ref_seq11_fasta, 3, 5),extracted_alt,substr(ref_seq11_fasta, 7, 9),sep="")

extracted_chr = as.character(df_conb$chrom)
extracted_chr_fasta = df_conb$chrom

extracted_pos = df_conb$pos

extracted_P1 = df_conb[,4]
extracted_P2 = df_conb[,5]


pos = df_conb$pos
chr = df_conb$chrom

near_to_indel = df_conb$near_to_indel

#i = 26
#j = 1


bp_dif = df_conb$pos - c(df_conb$pos[2:nrow(df_conb)],0 ) 
#as.numeric( 0 > bp_dif & -6 < bp_dif )
multiple_snp = as.numeric( 0 > bp_dif & -6 < bp_dif )


#str(ref_seq_fasta)

##Automatically process multiple restriction enzymes in a repetitive process
restriction_enzyme_list = read.csv(file="must/restriction_enzyme_list.csv", stringsAsFactors = F)

empty_df = as.data.frame(matrix(NA, nrow = nrow(df_conb), ncol = 14),stringsAsFactors = F)

x=1

i=1
for(i in 1:nrow(restriction_enzyme_list)){
     enzyme_name  = restriction_enzyme_list[i,1]
     enzyme_seq   = restriction_enzyme_list[i,2]
     enzyme_seq_length = nchar(enzyme_seq)
  
  if(enzyme_seq_length == 6){
    
    ref_TF  =  grepl(enzyme_seq, ref_seq11_fasta)
    alt_TF  =  grepl(enzyme_seq, alt_seq11_fasta)
    message(paste("restriction endonuclease recognition site length of",enzyme_name,"is 6"))
    
    re_enzyme_site  = as.numeric(ref_TF == alt_TF)
    
    
    ref_cut  = as.numeric(ref_TF)
    alt_cut  = as.numeric(alt_TF)
    
    
    
    df_res1 = data.frame(extracted_chr,extracted_pos,extracted_P1,extracted_P2,
                         ref_seq11_fasta,
                         alt_seq11_fasta,
                         ref_seq7_fasta,
                         alt_seq7_fasta,
                         ref_cut,alt_cut,re_enzyme_site,
                         enzyme_name,near_to_indel,  multiple_snp, stringsAsFactors = F)
    
    df_res2 = df_res1[df_res1$re_enzyme_site  == 0 ,]
    
    #str(df_res2)
    message(paste("you got",nrow(df_res2),"loci !!"))
    
    if( nrow(df_res2) > 0 ){
    empty_df[x:(x+nrow(df_res2)-1),] = df_res2
    }
    
    
    x = x+nrow(df_res2)
    
   }else if (enzyme_seq_length == 4){
    
    ref_TF  =  grepl(enzyme_seq, ref_seq7_fasta)
    alt_TF  =  grepl(enzyme_seq, alt_seq7_fasta)
    message(paste("restriction endonuclease recognition site length of",enzyme_name,"is 4"))
    
    re_enzyme_site  = as.numeric(ref_TF == alt_TF)
    ref_cut  = as.numeric(ref_TF)
    alt_cut  = as.numeric(alt_TF)
     
    df_res1 = data.frame(extracted_chr,extracted_pos,extracted_P1,extracted_P2,
                          ref_seq11_fasta,
                          alt_seq11_fasta,
                          ref_seq7_fasta,
                          alt_seq7_fasta,
                          ref_cut,alt_cut,re_enzyme_site,
                          enzyme_name,
                          near_to_indel,
                          multiple_snp, stringsAsFactors = F)
      

     df_res2 = df_res1[df_res1$re_enzyme_site  == 0 ,]
     
     #str(df_res2)
     message(paste("you got",nrow(df_res2),"loci !!"))
     
     if( nrow(df_res2) > 0 ){
     empty_df[x:(x+nrow(df_res2)-1),] = df_res2
     }
     
     
     x = x+nrow(df_res2)
     
     
     } else {
       message(paste("restriction endonuclease recognition site of",enzyme_name,"is not acceptable"))
     }
 }

##Summarize and save
final_df = empty_df[!is.na(empty_df[,1]),]
colnames(final_df) = colnames(df_res2)

final_df = final_df[order(final_df$extracted_pos),]
final_df = final_df[order(final_df$extracted_chr),]

message(paste("Finally, you got",nrow(final_df),"CAPS !!"))
write.csv(file="results/candidate_of_snp_amd_enzyme_list.csv",final_df)






