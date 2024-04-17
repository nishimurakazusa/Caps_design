#Summarize primer3 output

args <- commandArgs(trailingOnly = T)

# Unspecified options set default values
default_args <- c("900", "3")   # Default value definition
default_flg <- is.na(args[1:2])
args[default_flg] <- default_args[default_flg]  

#Set file name and parameters
extract_seq_length <- as.numeric(args[1])



#df_res4 = read.csv(file="results/candidate_of_snp_and_enzyme_list.csv",stringsAsFactors = F,row.names = 1)
#str(df_res4)
#chr = df_res4[,1]

##Data frame to store information
result_primers_list = matrix(NA,nrow = (nrow(df_res4)*10),ncol=13+14) 

i = 1

primer3_example = read.csv(file="must/primer3_v2.csv",stringsAsFactors = F,header = F)

n_para = nrow(primer3_example)

primer3_res = read.table(file=paste("primer3/output/result",i,sep=""),sep="=",stringsAsFactors = F)
  
s0 = which(primer3_res[,1]  == "PRIMER_PAIR_0_PENALTY")
e0 = which(primer3_res[,1]  == "PRIMER_PAIR_1_PENALTY")-1
s1 = which(primer3_res[,1]  == "PRIMER_PAIR_1_PENALTY")
e1 = which(primer3_res[,1]  == "PRIMER_PAIR_2_PENALTY")-1
s2 = which(primer3_res[,1]  == "PRIMER_PAIR_2_PENALTY")
e2 = which(primer3_res[,1]  == "PRIMER_PAIR_3_PENALTY")-1
s3 = which(primer3_res[,1]  == "PRIMER_PAIR_3_PENALTY")
e3 = which(primer3_res[,1]  == "PRIMER_PAIR_4_PENALTY")-1
s4 = which(primer3_res[,1]  == "PRIMER_PAIR_4_PENALTY")
e4 = which(primer3_res[,1]  == "PRIMER_PAIR_4_PENALTY") + e3 -s3

i=1
for(i in 1:nrow(df_res4)){
  for_massage = paste(i,"回目",sep="")
  message(for_massage)
  
  primer3_res = read.table(file=paste("primer3/output/result",i,sep=""),sep="=",stringsAsFactors = F)
  
  all_primer =  cbind(primer3_res[s0:e0,],
                      primer3_res[s1:e1,2], 
                      primer3_res[s2:e2,2],
                      primer3_res[s3:e3,2],
                      primer3_res[s4:e4,2])
  
  locus = as.character(primer3_res[1,2])
  seq_temp = as.character(primer3_res[2,2])
  
  enzyme_data =  paste0(df_res4[i,]) 
  enzyme_data_names =  colnames(df_res4) 
  
  
  for(j in 0:4){
    primerFname     = paste0(chr[i],"_",df_res4[i,2],"_",j,"F",sep="") 
    primerRname     = paste0(chr[i],"_",df_res4[i,2],"_",j,"R",sep="") 
    
    primerFseq      = as.character(all_primer[4,j+2])
    primerRseq      = as.character(all_primer[5,j+2])
    
    primerFpenalty  = as.character(all_primer[2,j+2])
    primerRpenalty  = as.character(all_primer[3,j+2])
    
    primerFstart      = strsplit(as.character(all_primer[6,j+2]),",")[[1]][1]
    primerRstart      = strsplit(as.character(all_primer[7,j+2]),",")[[1]][1]
    
    primerFTm  = as.character(all_primer[8,j+2])
    primerRTm  = as.character(all_primer[9,j+2])
    
    expected_product_size =  as.character(all_primer[22,j+2])
    
    primerF_digested  = as.character((extract_seq_length/2)- as.numeric(primerFstart)+1)
    primerR_digested  = as.character((as.numeric(primerRstart) - (extract_seq_length/2)))
    
    compl_any_th = as.character(all_primer[20,j+2])
    compl_end_th = as.character(all_primer[21,j+2])
    
    primerF_GC = as.character(all_primer[10,j+2])
    primerR_GC = as.character(all_primer[11,j+2])
  
    
    primerF = c(primerFname, primerFseq, primerFpenalty, compl_any_th,  compl_end_th, primerF_GC,
                primerFstart, primerFTm, primerF_digested, expected_product_size,j,locus,seq_temp, enzyme_data)
    primerR = c(primerRname, primerRseq, primerRpenalty, compl_any_th,  compl_end_th,  primerR_GC,
                primerRstart, primerRTm, primerR_digested, expected_product_size,j,locus,seq_temp, enzyme_data)
    
    
    result_primers_list[((i-1)*10+1+(j*2)),] =  primerF
    result_primers_list[((i-1)*10+2+(j*2)),] =  primerR
    
  }
  
}

colnames(result_primers_list) = c("Primer_names","sequences","penalty", "compl_any_th", "compl_end_th","primer_GC",
                                  "start","Tm","expected_dijested_band",
                                  "expected_product_size","No_of_primer","loci","seq_temp", enzyme_data_names)

result_primers_list2 = result_primers_list[,c(11,14,15,1,2,3,7,8,9,10,22,23,24,25,26,27,4:6,12:13,16:21)]

write.csv(file="result_primer_list.csv",result_primers_list2)

message("completed !!")


