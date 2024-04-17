#!/bin/zsh
#dependencies
#primer3 seqkit R
#Rpackage vcfR Biostrings


extract_length=900
amplicon_length_min=400
amplicon_length_max=600

##Obtain the directories
cd `dirname $0`

#Set the genome file
genome=`echo *.fasta`
vcf=`echo *.vcf`

rm -r primer3
rm -r df
rm -r results

##GTF generate
Rscript must/gtf_generater_for_SNP_loci.R $vcf


##Extract SNP loci
seqkit subseq $genome --gtf results/caps.gtf > results/caps.fasta

##restriction site search
Rscript must/restriction_enzyme_searcher.R

##Generate GTF for primer3
Rscript must/gtf_generater_for_primer3.R $extract_length 3

##Cut fasta for primer3
seqkit subseq  $genome --gtf results/caps_primer.gtf > results/seq_for_caps_primer.fasta

##Generate primer3 command
Rscript must/primer3_command_generater.R amplicon_length_min amplicon_length_max

##Run primer3_core
num=$(ls -1 primer3/input/ | wc -l)

for i in `seq 1 $num`; do
  primer3_core primer3/input/input$i > primer3/output/result$i
done

##Run primer3
Rscript must/summarize_primer_pairs.R
