# Haphpipe User Guide

## Basics:

*Made by: Shane Pruitt

*My virus is: HIV

*My links are:*

*Associated paper pararaphs.*
*Paragraph about study:*

*Paragraph abut NGS technique and data:*

## Questions


1. What is a directory structure? Please explain and create a diagram / picture to accompany your explaination.
A directory structure is a heirarchy of files, that give a pathway to your file destination. It allows you to sort files into different "pathways", partitionaing data into more and more specific categories/save locations the further down the path you get. IE: If the pathway is Desktop-->Pictures--->blackandwhitephotos, the desktop contains all files saved to the desktop, the pictures folder contains all the pictures youve saved but only the pictures, and blackandwhitephotos contains only the black and white photos that youve taken. As youve gone down the line, youve specified what you need into more specific and defined categories. This comes in handy when processing data, where you can seperate the files for use. 

2. What is the difference between next-generation sequencing and Sanger sequencing?
Sanger sequencing "creates" a sequence by sequencing one small piece of the "puzzle" at a time. While NGS may be able to decipher the most/large portions of the "puzzle" at a time. NGS is much more efficient, is cheaper, and may detect varients more efficiently than the Sanger method. 


3. What is the difference between the two pipelines, hp01 and hp02? Explain what each step of the pipelines accomplishes and why that step is necessary.

- hp_assemble_01

>This pipeline implements denovo assembly. Reads are first trimmed (trim_reads) and used as input for denovo assembly (assemble_denovo). The denovo assembly stage automatically performs error detection, but trimmed reads are also error corrected (ec_reads). The assembled contigs are used as input for amplicon assembly (assemble_amplicons) along with reference FASTA and GTF files. The assembly is then iteratively refined up to five times (refine_assembly) by mapping corrected reads to the assembled FASTA file and lastly finalized (finalize_assembly), resulting in a FASTA file with final consensus sequences, final VCF, and aligned BAM file.

  

- hp_assemble_02  
>This pipeline implements reference-based mapping assembly. Reads are first trimmed (trim_reads) and error-corrected (ec_reads). The corrected reads are used as input for reference-based mapping assembly (refine_assembly) for up to five iterations. Lastly, the assembly is finalized (finalize_assembly) by mapping reads onto the refined reference sequence. The final output is a FASTA file with final consensus sequences, final VCF, and aligned BAM file.

4. Did you need to learn how to bash script for this? What do you feel like you weren't prepared for after reading the introduction part we provided you with? Yes I needed to learn how to bash script for this. I had previously had no experience whatsoever in using a command line, scripting, or anything of that nature. After going through the material I was able to gain a better understanding of what was going on, but I defenitely needed the help of someone more experienced to get me through the process. For me some of the stuff was literally "bashing" at the keyboard, experimenting to see if it worked or not. It probably took a good amount of time more to figure it out than my counterpart Freddie. 



## Installing HAPHPIPE
These are some of the suite of tools that haphpipe utilizes. Inserting these commands into the terminal "tells" the computer that when running haphpipe, it will be always utilizing these tools, somewhat like installing a program. All of these lines perform a certain function, and by creating an environement it packages all of these functions together. This is more efficient so that you dont need to run each command over and over again when analyzing data. 


conda create -n haphpipe/
 python \
 future \
 pyyaml \
 biopython \
 seqtk \
 bowtie2 \
 bwa \
 flash \
 freebayes \
 mummer \
 picard \
 trimmomatic \
 samtools=1.9 \
 gatk=3.8 \
 spades \
 blast \
 sierrapy
 
## Activating HAPHPIPE
Input the code below into your terminal. This will tell the conda that you want to activate the program Haphpipe. 

	conda activate haphpipe
	
While trying to register the Gatk file, I kept getting the error "Error: Unable to access jarfile ./GenomeAnalysisTK.jar
The version of the jar specified, , does not match the version expected by conda: 3.8". This was because someone else had already downloaded it on the computer a couple years ago. 
this was fixed by changing the directory to downloads and putting the path to the original file download. 

	gatk3-register /User/crandalllab/Downloads/GenomeAnalysisTK-3.8-0-ge9d806836.tar.bz2
	
## HAPHPIPE suite
Please explain what the haphpipe suite is. What is it's purpose how do you use it? What is it good for?

For each of the stages below, please explain the stage. List the files needed for the command, what the command does. The options for the command. An example of how to execute the command.

### hp_reads
        sample_reads
        trim_reads
        join_reads
        ec_reads
### hp_assemble
        assemble_denovo
        assemble_amplicons
        assemble_scaffold
        align_reads
        call_variants
        vcf_to_consensus
        refine_assembly
        finalize_assembly
### hp_haplotype
        predict_haplo
        ph_parser
### hp_annotate
        pairwise_align
        extract_pairwise
        annotate_from_ref
        
## Example usage
Your virus with both pipelines. Document all code and explanation.  

## Helpful resources
List all helpful topics you think or did address here with links / explainations you felt were helpful. Use bullet points.

- introduction to NGS
	- https://learn.gencore.bio.nyu.edu/variant-calling/
- inroduction to bash
- intro to github
	- https://github.com/gwcbi/HPC/blob/master/github.md
- markdown tutorial
	- https://www.markdowntutorial.com
- tips for command line
	- https://github.com/gwcbi/HPC/blob/master/commandline.md
- how to run interactive jobs *(do this for all haphpipe modules and pipelines!)*
	- https://github.com/gwcbi/HPC/blob/master/interactive_jobs.md
	
## FAQ
List all questions you had and provide links / explainations you felt were helpful to answer these questions. Use topics and numbering.
