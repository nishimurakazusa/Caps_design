# Caps_design

## Script for semi-automatic design of CAPS markers from vcf file  

dependencies
Biostrings  <https://bioconductor.org/packages/release/bioc/html/Biostrings.html>  
vcfR        <https://cran.r-project.org/web/packages/vcfR/index.html>  
Seqkit      <https://bioinf.shenwei.me/seqkit/>  
Primer3     <https://github.com/primer3-org/primer3>  

Vcf file in the CAPS_marker_generater directory for MacOS and, If you place the reference genome file in the CAPS_marker_generater directory and double-click the Shell script 'vcf_caps_generater.command', the resulting file will be automatically output. For Linux OS, just execute the script in vcf_caps_generater,comman and the process will proceed.
The following is an explanation of vcf_caps_generater.command. First, set the length of the sequence to be extracted from the fasta file of the genome and the length of the amplified fragment for primer design.

    extract_length=900  
    amplicon_length_min=400  
    amplicon_length_max=600  


Next, get the path of the vcf and reference genome files and, if running as a Shell script, specify the paths of the corresponding files.

    `genome=`echo *.fasta``  
    `vcf=`echo *.vcf``  

The R script gtf_generater_for_SNP_loci.R is used to create a gtf file to obtain 11 bp sequences before and after the position of polymorphisms from the reference genome.

    `##GTF generate`  
    `Rscript must/gtf_generater_for_SNP_loci.R $vcf`  
    
    `##Extract SNP loci`  
    `seqkit subseq $genome --gtf results/caps.gtf > results/caps.fasta`  


Next, based on the sequence information around the extracted polymorphism, the restriction enzyme cleavage and non-cleavage of the reference and mutant sequences are determined. The list of restriction enzymes is listed in must/restriction_enzyme_list.csv. This is a two-column csv file, where the first column is headed "enzyme" and the second column is headed "seq", indicating the enzyme names and recognition sequence of the restriction enzyme, respectively.
The resulting file is output as 'candidate_of_snp_amd_enzyme_list.csv' in the result folder.


    `##restriction site search`  
    `Rscript must/restriction_enzyme_searcher.R`  

Next, gtf_generater_for_primer3.R creates a gtf file for extracting sequences around loci that can be used to design primer set for CAPS, and the sequences are extracted from the reference genome by seqkit subseq.

    `##Generate GTF for primer3`  
    `Rscript must/gtf_generater_for_primer3.R $extract_length 3`  

    `##Cut fasta for primer3`  
    `seqkit subseq $genome --gtf results/caps_primer.gtf > results/seq_for_caps_primer.fasta`  


Next, the extracted fasta file is used to create the input file for primer 3. At this time, the file must/primer3_v2.csv is used, and the number of lines in this file must not be changed. The input file is created by primer3_command_generater.R in the input of the primer3 directory. Then, the iterative process executes primer3 and the results are saved in primer3/output.


    `##Generate primer3 command`  
    `Rscript must/primer3_command_generater.R amplicon_length_min amplicon_length_max`  
    
    `##Run primer3_core`  
    `num=$(ls -1 primer3/input/ | wc -l)`  
    
    `for i in `seq 1 $num`; do`  
    `  primer3_core primer3/input/input$i > primer3/output/result$i`  
    `done`  


Script to summarize Primer3 results. After this script run, result_primer_list.csv is output.
    `##Run primer3`
    `Rscript must/summarize_primer_pairs.R`

![image](https://github.com/nishimurakazusa/Caps_design/assets/46805695/42fc35a7-80d2-402b-8a87-4773e9b0c545)
