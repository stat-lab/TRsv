# TRsv
Detection and Genotyping of Tandem Repeat Expansion/Contraction, Structural Variants (SVs), and Indels using Long Reads  

## Introduction

Tandem repeat variations (TR-CNVs) are copy number variations of repeat units of tandem repeats (TRs). TRsv detects TR-CNVs within TR regions and structural variations (SVs) and short indels outside TR regions separately using long reads. The TR regions are defined using either a pre-built human TR bed file or user-prepared TR bed file. TRsv improves TR-CNV detection by assembling fragmented insertion and deletion alignments within a TR region of one long read alignment and searching for the TR units in the insertion sequence. In addition, TRsv improves SV detection derived from non-HiFi reads (PacBio CLR and Nanopore long reads) through machine learning-based filtering.  

## Requirements

perl 5 or later versions  
R 3.5 or later (required library: [xgboost, Matrix, data.table, jsonlite, and lattice])(https://cran.r-project.org/web/packages/xgboost/)  
[Samtools](https://github.com/samtools/samtools)  
[Minimap2](https://github.com/lh3/minimap2)  
[YASS](https://bioinfo.univ-lille.fr/yass/)  
[TRF](https://github.com/Benson-Genomics-Lab/TRF)  
[MultAlin](http://lipm-bioinfo.toulouse.inrae.fr/download/multalin/)  
The above install directories must be set to PATH or specified with the corresponding tool_path options to enable to execute perl, Rscript, samtools, minimap2, yass, trf, and multalin.   
- Tips for Multalin installation  
  ```
  wget http://lipm-bioinfo.toulouse.inrae.fr/download/multalin/multalin.5.4.1.dynamic.sources.tar.gz
  gzip -dc multalin.5.4.1.dynamic.sources.tar.gz | tar xvf -
  cd multalin.5.4.1
  cd src
  cmake .
  make
  cd ../
  ln -s src/multalin
  chmod -R 755 $PWD
  export PATH=$PWD:$PATH (setenv PATH $PWD:$PATH # for csh shell)
  export MULTALIN=$PWD (or export MULTALIN=$PWD/matrix) (setenv MULTALIN $PWD # for csh shell)
  ```

#### Input file
- Alignment bam/cram file & its index file for PacBio HiFi, PacBio CLR or ONT Nanopore long reads  
- Reference fasta file & its index file  
- Reference TR bed file (Human TR bed files for hg37, hg38, and CNM13-T2T are included in this package)  
  
- gap bed file (gap regions in the input reference genome. Human files for each build are included in this package. optional) 
- transposable element (TE) fasta file (Human files containing ALU, LINE1, and SVA consensus sequences are included in this package. optional but highly recommended) 
- excluded bed file (genome regions to be excluded. Human files for each build are included in this package. optional) 
- Gender list file (optionally required for joint call, 1st column: sample name, 2nd column: M|F, separated by tab)
- GFF annotation file (required for annotation command, available from [Ensembl](http://www.ensembl.org/). Human files for each build are included in this package.)  
  
[Format of reference TR bed file]  
  
The following 8 columns are needed for each TR region.  
1st column (Chr): chromosome name  
2nd column (Start): start position of a TR region  
3rd column (End): end position of a TR region  
4th column (TRID): ID/name of TR  
5th column (UnitSize): repeat unit size of TR  
6th column (CN): repeat unit copy number of TR  
7th column (ME): mobile element, such as ALU and LINE1, that has homology to TR repeat unit ('-' is indicated if no homology, it's OK even if all lines have '-'.)  
8th column (UnitSeq): sequence of TR repeat unit  


## Citation

Unpublished (submitted).  


## Install

Download the latest release from https://github.com/stat-lab/TRsv/releases, and unpack with gzip -dc TRsv_v1.*.tar.gz | tar xvf -  

The Data folder in the TRsv package contains tandem repeat bed files, gap bed files, centromere bed files, and gff annotation files for the human build 37/38/T2T references. The Data folder also contains the training data sets for machine learning-based SV filtering for non-HiFi data. Do not change the name of the files/directories (except config files) and the directory structure within the TRsv package.  

Install the required external tools (samtools, yass, multalin, trf, and Rscript). If the executable names are different from 'trf', 'yass', or 'multalin', please rename the executable names to 'trf', 'yass', or 'multalin' (e.g., trf409.linux64 -> trf). Set the path (the directory path containing the executable) to the $PATH environmental variable (e.g. export PATH=/home/tools/yass-1.15/bin:$PATH). Alternatively, specify the options (samtool_path, yass_path, multalin_path, trf_path, and r_path) on TRsv command or in configure file (e.g., --yass_path /home/tools/yass-1.15/bin).  

Alternatively, a singularity container file (TRsv-v1.1.sif) is available at [http://jenger.riken.jp/static/SVs_bykosugisensei_20220329/TRsv-v1.1.sif]. Run using the singularity TRsv-v1.1.sif file can be done using a run_TRsv_singularity.pl script in the scripts directory.

## <a name="hdata"></a>Human and Non-human Data 

By default, TRsv handles WGS alignment data (bam/cram) generated based on the human build 37 reference (GRCh37 or GRCh37d5). To use the data based on the human build 38 reference or T2T-CHM13.v2.0, run the TRsv command with the ‘--build 38’ or '--build T2T' option. When the --build option is specified, specific data files related to the human reference will be automatically selected from the Data folder. If you need to use files other than those in the Data folder, please specify the relevant options, such as --repeat_bed and --exclude_bed. For non-human species, use the option ‘-nh 1’ and use related options such as --repeat_bed, --gap_bed, and --te_fa to specify external input files such as repeat bed, gap bed, TE fasta file, etc.  

## <a name="gusage"></a>General Usage

A description of the options used in each command is output with the following commands:  
```
TRsv call -h -x hifi (-x ont/clr for non HiFi data)
TRsv joint_call -h
TRsv annotate -h
```

### Call TR-CNVs and SVs/indels

The variant calling consists of four steps and is performed sequentially. In the first step (step-1), variants are collected from the input alignment file, and in the second step (step-2), the output file from step-1 is used to characterize the insertion sequences (e.g., determining the repeat unit content in TR-INS, checking the homology with TE sequences). The third step (step-3) is to annotate the low quality variants in the output file from step-2. The last step removes plausible false positive calls using a machine learning method and outputs the final vcf file. Intermediate files from step-1 to step-3 are output to a ${out_prefix}.temp folder, with ${$out_prefix}.chr*.discov.txt from step-1 and ${$out_prefix}.chr*.discov.out from steps-2/3.  

[When using configure file]  
```
TRsv call -c config.txt
```
(Use template configure files for HiFi and non-HiFi contained in the package.)  
  If you use human data and intend to use its related files contained in the Data folder, you do not need to specify 'repeat_bed', 'repeat_u', 'lowconf_tr', 'te_fasta', and 'gap_bed' in the config file.

[When not using configure file (for human)]  
```
TRsv call -b <bam_file> -r <reference_fasta> -x <platform, hifi|ont|clr> --build <human_ref_build, 37|38|T2T, default: 37> -p <output_prefix> -n <number of threads>
```
  In this case, a reference tandem repeat bed and other related files specific to the human reference build are automatically supplied from the Data folder.  
  To use user-specified files other than these files, specify the relevant options, such as --repeat_bed, --exclude_bed, --te_fa, etc.  
  The platform is either hifi (PacBio HiFi), clr (PacBio CLR) or ont (ONT Nanopore).  
  The human_ref_build is either 37 (GRCh37), 38 (GRCh38) or T2T (CHM13-T2T).  
  The jobs in steps 1-3 are executed in parallel for the N number of chromosome data specified by the -n option.

[When not using configure file (for non-human)]  
```
TRsv call -b <bam_file> -r <reference_fasta> -x <platform, hifi|ont|clr> -p <output_prefix> -n <number_of_threads> -nh 1 -rep <repeat_bed (mandatory)> -gb <gap_bed (optional)> -tf <TE_fasta (optional)> 
```
  The repeat_bed is a bed file describing tandem repeat regions and tandem repeat information (see "Input file" section)
  The gap_bed is a bed file describing tab-separated gap regions (chr start end) in each line.  
  The TE_fasta is a fasta file of transposable elements.  
  
[When using singularity container TRsv.sif] 
```
singularity exec -B <list of bind directories for input file path on the host> <absolute path of TRsv sif file> TRsv call -x <platform, hifi|ont|clr> --build <human_ref_build, 37|38|T2T> -r <absolute path of reference fasta> -b <absolute path of bam file> -p <output_prefix> -n <number_of_threads>
```
Alternative using run_TRsv_singularity.pl script.  
```
run_TRsv_singularity.pl -sif <absolute path of TRsv sif file> -com <TRsv command, call|joint_call|annotate> [other arguments; specified with -x, -b, -r, --build, -p, -n, and other command-specific options]
```
The input files must be specified with the absolute path of the real files (not linked files). When using the run_TRsv_singularity.pl script, the bind directories for the input file path (for singularity -B option) are automatically specified in the script.

Other frequently used options are as follows:  

-ml (min_len) <minimum length of variant> (default: 3 for HiFi reads, 20 for non-HiFi reads)  
-msl (min_tr_len) <minimum length of TR-CNV> (default: 3 for HiFi reads, 20 for non-HiFi reads)  
If this is not specified, the value specified with -ml is used for -msl.  
-xc (exclude_chr) <chromosome name(s) to be excluded, comma-separated> (e.g., -xc Y for male sample)  
-sk (skip) <skip step-1~3> (1: skip step-1, 2: skip step-2, 3: skip step-3, 12: skip steps-1 and -2, 123: skip steps-1, -2, and -3)  
You can use this option if you want to resume from a step in the middle of the process.  
-rp (r_path) path of R (>= v3.5), where xgboost library has been installed if the corresponding R is not set in $PATH (only non-HiFi data)  



#### Output files  

${output_prefix}.discov.vcf (final output vcf file for HiFi data)  
${output_prefix}.discov.filt.vcf (final output vcf file for non-HiFi data after machine learning-based filtering)  
${output_prefix}.discov.INS.fa (fasta file of insertion sequences of TR-INSs inside TR regions and INSs outside TR regions)  
  
In the output vcf file, the TR-CNV line has the 'TR:CNV' string in the fifth field and 'SVTYPE', 'SVLEN', 'CN', 'TRID', 'TREND', 'TRULEN', and other tags in the eighth field. The SVTYPE tag in the TR-CNV line indicates INS (expansion) or DEL (contraction) of the corresponding TR unit. The SVLEN and CN tags represent the INS/DEL length and copy number of the corresponding TR unit, respectively. If there are two non-reference TR-CNV alleles, these alleles are represented with comma-separated values or strings for each SVTYPE, SVLEN, and CN tag. The TRID, TREND, TRULEN tags represent the TR ID, end position, and repeat unit size of the corresponding TR, respectively. A TR-INS with an unrelated INS sequence or low TR unit content (< 50% by default) is assigned as a normal INS with the TRID tag accompanied with a ME or TRUNIT tag in the eighth field. Conversely, when the sequence of a normal INS contains a high content of copies of some tandem repeat unit, the INS line has the TRUNIT tag indicating the TR unit sequence, but no TRID tag.

### Merge vcf files from multiple samples (optional)

[Human data]

```
TRsv joint_call -v <input_vcf_list> -p <output_prefix> -od <output_directory, default: working directory> -g <geneder_list>

```
  (--build 38 for human GRCh38, --build T2T for CHM13-T2T)
  The input_vcf_list is a file in which each line describes tab-separated sample name and the path of the corresponding vcf file.  
  The gender_list is a file containing tab-separated sample name and the corresponding gender (M: male or F: female) in each line.

[Non-humnan data]

```
TRsv joint_call -v <input_vcf_list> -p <output_prefix> -od <output_directory> -nh 1 -r <reference_fasta_index> -gap <gap_bed, optional>

```
  The reference_fasta_index is an index file of the reference fasta used, which can be generated with the samtools faidx command.
  The gap_bed is a bed file describing tab-separated gap regions (chr start end) in each line.  

For large data sets with a large number of samples and/or variants, running the command for each chromosome using the -c option should shorten the run time.

#### Output file

${output_prefix}.All-samples.vcf


### Annotate variants overlapping gene regions (optional)

The annotate command adds gene name/ID and gene region that overlap with the TR regions or SVs/indels to the INFO filed (with SVANN key) of the vcf file. The gene regions include exon/CDS (All-exons if the TR or SV completely overlaps all exons), 5’-/3’-UTR, intron, 5’-/3’-flanking region. Two ranges (5 Kb and 50 Kb) of the flanking regions are specified by default, and these lengths can be changed with the options, -c5, -c3, -f5, and -f3. These annotations are also added to the FORMAT AN subfield for each sample in an additional output vcf file. For human, one of the gff3 gene annotation files (Homo_sapiens.GRCh37.87.gff3.gz, Homo_sapiens.GRCh38.104.gff3.gz, or Homo_sapience.T2T-chm13v2.0.ensemble.gff3.gz), downloaded from Ensembl (ftp://ftp.ensembl.org/pub/grch37/release-87/gff3/homo_sapiens9), is automatically selected by default. For non-human species, a gff3 annotation file obtained from the Ensembl site must be specified with the -r option. Any input SV vcf file with SVTYPE and SVLEN keys in the INFO field can be used. The annotate command can be done as follows:  
<Human data>  

```
TRsv annotate -v <input_vcf> -p <output_prefix>
```
  (--build 38 for human GRCh38, --build T2T for CHM13-T2T)  
  The bp-size specification of gene flanking region from the terminal exons
  (-c5 <5'_proximal_region, default: 5000> -c3 <3'_proximal_region, default: 5000> -f5 <5'_distal_region, default: 50000> -f3 <3'_distal_region, default: 50000>)  

<Non-human data>  

```
TRsv annotate -v <input_vcf> -p <output_prefix> -nh 1 -r <gff3_file>
```

#### Output files  
 
${out_prefix}.annot.vcf  
${out_prefix}.AS.annot.vcf (Annotations for each sample are added to the FORMAT AN subfield.)  


## <a name="quick"></a>Quick start with test data

The test_data folder contains a bam file of NA12878 PacBio HiFi reads aligned to the GRCh38 reference, corresponding to the first 500 Kb of chr20. A fasta file of the GRCh38 chr20 and a config file are also included in this folder.
```
cd TRsv_v1.0/test_data
```
### Variant call
```
../TRsv call -c config_hifi.txt
```
### Variant annotation
```
../TRsv annotate -v test.discov.vcf -p test --build 38
```
