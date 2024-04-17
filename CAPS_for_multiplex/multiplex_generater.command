#!/bin/zsh
#dependency
#primer3 seqkit R
#Rpackage vcfR Biostrings


min=100
max=1800
ex_length=2000

##Obtaining directories
cd `dirname $0`

#Get vcf and genome file name
genome=`echo *.fasta`
#vcf=`echo *.vcf`

rm -r primer3
rm -r df
rm -r results

##Rscript to generate GTFs to cut out sequences for Primer3
Rscript must/gtf_generater_for_primer3_for_multiplex.R $ex_length

##Extract the sequence to be input to Primer3
seqkit subseq  $genome --gtf results/caps_primer.gtf > results/seq_for_caps_primer.fasta

##Automatically generate commands to throw to primer3_core
##The length of the primer3 target array can be specified as the first-second argument.
##Default is 500-550
Rscript must/primer3_command_generater.R $min $max $ex_length

##primer3_core to make a large number of primers
##The last code says how many files were created, so let's make as many as that number.
num=$(ls -1 primer3/input/ | wc -l)

for i in `seq 1 $num`; do
  primer3_core primer3/input/input$i > primer3/output/result$i
done

##Summarize primer3
Rscript must/summarize_primer_pairs2.R
