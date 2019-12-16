# Haphpipe User Guide

## Basics:

*Made by: Shane Pruitt

*My virus is: HIV

*My links are:*

*Associated paper pararaphs.*
The authors are claiming that although NGS may detect more low frequency variants, this may not be useful in standard care. They claim that there was already 90% viral suppression rate in cities like Seattle, so sanger sequencing was already working at a very effective rate in the absence of NGS. According to the study, both tests were equally able to detect high frequency variants.PASeq and HyDRA were just as good at detecting high frequency DRMS as sanger sequencing was, showing their potential worth in analysis pipelines. They claim that although better at detecting low frequency variations, NGS doesn't provide much more utility in the workplace than what is already put into place with Sanger sequencing. They do however note that further studies are needed in areas that have not achieved the 90% threshold. 

## Questions


1. What is a directory structure? Please explain and create a diagram / picture to accompany your explaination.
A directory structure is a heirarchy of files, that give a pathway to your file destination. It allows you to sort files into different "pathways", partitionaing data into more and more specific categories/save locations the further down the path you get. IE: If the pathway is Desktop-->Pictures--->blackandwhitephotos, the desktop contains all files saved to the desktop, the pictures folder contains all the pictures youve saved but only the pictures, and blackandwhitephotos contains only the black and white photos that youve taken. As youve gone down the line, youve specified what you need into more specific and defined categories. This comes in handy when processing data, where you can seperate the files for use. 

2. What is the difference between next-generation sequencing and Sanger sequencing?
Sanger sequencing "creates" a sequence by sequencing one small piece of the "puzzle" at a time. While NGS may be able to decipher the most/large portions of the "puzzle" at a time. NGS is much more efficient, is cheaper, and may detect varients more efficiently than the Sanger method. 


3. What is the difference between the two pipelines, hp01 and hp02? Explain what each step of the pipelines accomplishes and why that step is necessary. It seems that hp01 doesn't have to go through the correction stage because the denovo assembly does that during its process, and then the contigs are oriented to the amplicon reference. hp02 goes through the same first step, but because it does not go through denovo assembly, it goes through a correction phase. After these steps are done, both pipelines go through the same steps/commands. 

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

To me a pipeline is a sequence of commands that perform certain tasks in the order that you give it to them, for the purpose of analyzing data. Conda packages these programs for you into one simpler command so that it's more efficient. The haphpipe suite to me is a set of commands put into a specific order for the purpose of analyzing sequences, creating a package of commands that streamline what would otherwise be a tedious process. 

### hp_reads
Seqtk is a fast and lightweight tool for processing sequences in the FASTA or FASTQ format. It seamlessly parses both FASTA and FASTQ files which can also be optionally compressed 
        sample_reads
        trim_reads
        join_reads
        ec_reads
### hp_assemble
"The denovo assembly stage automatically performs error detection, but trimmed reads are also error corrected". To me, this are programs that take the trimmed reads and clean them up in a way that they can be compiled together into a readable piece. 
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
### I pasted this in here because for me reading these steps actually gave me insight on what was going on. What I did was open the wiki on one half of the screen and then had these steps on the other, and slowly worked my way down to try and figure out what was going on. 
### Stages
Each stage can be run on its own. Stages are grouped into 4 categories: hp_reads, hp_assemble, hp_haplotype, and hp_annotate.

More detailed description of command line options for each stage are available in the [wiki](https://github.com/gwcbi/haphpipe/wiki). 

To view all available stages in HAPHPIPE, run:

```
haphpipe -h

```

  

  

  

### Reads

  

Stages to manipulate reads and perform quality control. Input is reads in FASTQ format, output is modified reads in FASTQ format. This to me states that this program is looking for "anomolies" or formatting issues that shouldnt normally be present. 

  

##### sample_reads

  

Subsample reads using seqtk ([documentation](https://github.com/lh3/seqtk)). Input is reads in FASTQ format. Output is sampled reads in FASTQ format.

Example to execute:

```

haphpipe sample_reads --fq1 read_1.fastq --fq2 read_2.fastq --nreads 1000 --seed 1234

```

  

##### trim_reads

  

Trim reads using Trimmomatic ([documentation](http://www.usadellab.org/cms/?page=trimmomatic)). Input is reads in FASTQ format. Output is trimmed reads in FASTQ format. As the documentation states this just cleans up the data in a way that makes it more manageable for analysis. If the reads are found to be below a certain quality threshold, then they will be trimmed. It allows for a better quality sequence. 

Example to execute:

```

haphpipe trim_reads --fq1 read_1.fastq --fq2 read_2.fastq

```

  

##### join_reads

  

Join reads using FLASH ([paper](https://www.ncbi.nlm.nih.gov/pubmed/21903629)). Input is reads in FASTQ format. Output is joined reads in FASTQ format. For denovo assembly, it can cause problems if there is a large amount of short sequence reads. Joining the reads is a process that will overlap the short reads so that they can be in a long read that better accomodates the assembly programs. As the paper states it has an incredibly high success rate with a very low margin of error, meaning that the joined reads are very accurate. 

Example to execute:

```

haphpipe join_reads --fq1 trimmed_1.fastq --fq2 trimmed_2.fastq

```

  

##### ec_reads

  

Error correction using SPAdes ([documentation](http://cab.spbu.ru/software/spades/)). Input is reads in FASTQ format. Output is error-corrected reads in FASTQ format. For me the documentation provided didn't really help with figuring out what SPAdes does. I also searched on google, but for the purposes of this user guide, there was nothing there that explained things at a basic level that I could find. 

Example to execute:

```

haphpipe ec_reads --fq1 trimmed_1.fastq --fq2 trimmed_2.fastq

```

  

### Assemble

  

Assemble consensus sequence(s). Input reads (in FASTQ format) are assembled

using either denovo assembly or reference-based alignment.

Resulting consensus can be further refined.

Since this is also utilizing SPAdes, my input from the previous section is repeated here. I would provide maybe a personalized explanation that is simpler (if the target audience is not experienced with computers such as myself). 
  

##### assemble_denovo

  

Assemble reads via de novo assembly using SPAdes ([documentation](http://cab.spbu.ru/software/spades/)). Input is reads in FASTQ format. Output is contigs in FNA format.

Example to execute:

```

haphpipe assemble_denovo --fq1 corrected_1.fastq --fq2 corrected_2.fastq --outdir denovo_assembly --no_error_correction TRUE

```

##### assemble_amplicons

  

Assemble contigs from de novo assembly using both a reference sequence and amplicon regions with MUMMER 3+ ([documentation](http://mummer.sourceforge.net/manual/)). Input is contigs and reference sequence in FASTA format and amplicon regions in GTF format. The documentation provided here is quite useful in understanding what the function is. It allows for quick assembly for large amounts of DNA or amino acid sequences. 

Example to execute:

```

haphpipe assemble_amplicons --contigs_fa denovo_contigs.fa --ref_fa refSequence.fasta --ref_gtf refAmplicons.gtf

```

  

##### assemble_scaffold

  

Scaffold contigs against a reference sequence with MUMMER 3+ ([documentation](http://mummer.sourceforge.net/manual/)). Input is contigs in FASTA format and reference sequence in FASTA format. Output is scaffold assembly, alligned scaffold, imputed scaffold, and padded scaffold in FASTA format. Since this section also uses MUMMER 3+ I would repeat what I said last section.

Example to execute:

```

haphpipe assemble_scaffold --contigs_fa denovo_contigs.fa --ref_fa refSequence.fasta

```

  

##### align_reads

  

Map reads to reference sequence (instead of running de novo assembly) using Bowtie2 ([documentation](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)) and Picard ([documentation](https://broadinstitute.github.io/picard/)). Input is reads in FASTQ format and reference sequence in FASTA format. Bowtie2 is an aligning tool that references the data provided by the user to long reference sequences. While Picard is a command line tool that helps you manipulate HTS data (basically that means that you can sequence several different DNA sequences at the same time, meaning that you can sequence large amounts of data concurrently)

Example to execute:

```

haphpipe align_reads --fq1 corrected_1.fastq --fq2 corrected _2.fastq --ref_fa refSequence.fasta

```

  

##### call_variants

  

Variant calling from alignment using GATK ([documentation](https://software.broadinstitute.org/gatk/download/archive)). Input is alignment file in BAM format and reference sequence in FASTA format (either reference from reference-based assembly or consensus final sequence from de novo assembly). Output is a Variant Call File (VCF) format file. At the GATk website there is a very nice user guide that covers almost anything youd like to know about the program. This is also a set of command line tools that allow for manipulation of HTS, just like the previous program. 

Example to execute:

```

haphpipe call_variants --aln_bam alignment.bam --ref_fa refSequence.fasta

```

  

##### vcf_to_consensus

  

Generate a consensus sequence from a VCF file. Input is a VCF file. Output is the consensus sequence in FASTA format.

Example to execute:

```

haphpipe vcf_to_consensus --vcf variants.vcf

```

  

##### refine_assembly

  

Map reads to a denovo assembly or reference alignment. Assembly or alignment is iteratively updated. Input is reads in FASTQ format and reference sequence (assembly or reference alignment) in FASTA format. Output is refined assembly in FASTA format.

Example to execute:

```

haphpipe refine_assembly --fq_1 corrected_1.fastq --fq2 corrected_2.fastq --ref_fa refSequence.fasta

```

  

##### finalize_assembly

  

Finalize consensus, map reads to consensus, and call variants. Input is reads in FASTQ format and reference sequence in FASTA format. Output is finalized reference sequence, alignment, and variants (in FASTA, BAM, and VCF formats, respectively).

```

haphpipe finalize_assembly --fq_1 corrected_1.fastq --fq2 corrected_2.fastq --ref_fa refined.fna

```

  

### Haplotype

  

Haplotype assembly stages. HAPHPIPE implements PredictHaplo ([paper](https://www.ncbi.nlm.nih.gov/pubmed/26355517)), although other haplotype reconstruction programs can be utilized outside of HAPHPIPE using the final output of HAPHPIPE, typically with the final consensus sequence (FASTA) file, reads (raw, trimmed, and/or corrected), and/or final alignment (BAM) file as input.

  

##### predict_haplo

  

Haplotype identification with PredictHaplo. Input is reads in FASTQ format and and reference sequence in FASTA format. Output is the longest global haplotype file and corresponding HTML file. __Note: PredictHaplo must be installed separately before running this stage.__

Example to execute:

```

haphpipe predict_haplo corrected_1.fastq --fq2 corrected_2.fastq --ref_fa final.fna

```

  

##### ph_parser

  

Return PredictHaplo output as a correctly formatted FASTA file. Input is the output file from __predict_haplo__ (longest global .fas file). Output is a correctly formatted FASTA file.

Example to execute:

```

haphpipe ph_parser best.fas

```

  

### Annotate

  

Stages to annotate consensus sequences.

  

##### pairwise_align

  

Apply correct coordinate system to final sequence(s) to facilitate downstream analyses. Input is the final sequence file in FASTA format, a reference sequence in FASTA format, and a reference GFT file. Output is a JSON file to be used in __extract_pairwise__.

Example to execute:

```

haphpipe pairwise_align --amplicons_fa final.fna --ref_fa refSequence.fasta --ref_gtf referenceSeq.gtf

```

  

##### extract_pairwise

  

Extract sequence regions from the pairwise alignment produced in __pairwise_align__. Input is the JSON file from __pairwise_align__. Output is either an unaligned nucleotide FASTA file, an aligned nucleotide FASTA file, an amino acid FASTA file, an amplicon GTF file, or a tab-separated values (TSV) file (default: nucleotide FASTA with regions of interest from GTF file used in __pairwise_align__).




## Helpful resources
List all helpful topics you think or did address here with links / explainations you felt were helpful. Use bullet points.
To be honest, I needed all of these references because I had zero experience with any of these. 
- introduction to NGS
	- https://learn.gencore.bio.nyu.edu/variant-calling/
- inroduction to bash from Github
- intro to github
	- https://github.com/gwcbi/HPC/blob/master/github.md
- markdown tutorial
	- https://www.markdowntutorial.com
- tips for command line
	- https://github.com/gwcbi/HPC/blob/master/commandline.md
- how to run interactive jobs 
	- https://github.com/gwcbi/HPC/blob/master/interactive_jobs.md
-HyDRA's website has a really nice example of what a user friendly user manual might look like. Although it still assumes that you are pretty knowledeable about computers, I think it might be a good place to start. [HyDRA's Website](https://github.com/Findwise/Hydra)
	
## FAQ
List all questions you had and provide links / explainations you felt were helpful to answer these questions. Use topics and numbering.

I was always taught to bring explanations down to the very basic level, and then bring them up as you go in order to prevent any confusion. Athough I dont usually use sites like Quora to me this quora question is a good example of a very simple, easy explanation that can be foundational to learning more. Even if the topic is not very complex, I've always found "explaining it like I'm 5" to be a good way to start. This was a a good resource for me on day 1. 
[Quora Answer](https://www.quora.com/What-are-pipelines-in-Bioinformatics)

My experience with programs/pipelines like this was previously little to none. When I first started I had to look up what a pipeline even was. It was a process that just took time to learn what was going on. Once you understand the basic principle, it makes it much easier to get the hang of whats going on. My groupmate Freddie is CS master student, and luckily he showed me running through his data, and analyzing it so I was able to get a better understanding of what was going on. Although the guides were somewhat helpful, I wasn't very acquainted with bashing, command lines, etc. Going through the process with Freddie to me was much more helpful. 
