# Dartmouth CQB RNA-seq analysis pipeline

## Introduction 
The pipeline is designed to provide efficient pre-processing and quality control of bulk RNA-sequencing (RNA-seq) data on high performance computing clusters (HPCs) leverging the Torque/PBS scheduler, and has been made available by the *Data Analytics Core (DAC)* of the *Center for Quantitative Biology (CQB)*, located at Dartmouth College. Both single- and paired-end datasets are supported, in addition to both library preparation methods interrogating full-length transcripts as well as 3'-end profiling methods. The pipeline has been built and tested using human and mouse data sets. Required software can be installed using Conda with the enrionment file (environment.yml) located in *Dartmouth-Data-Analytics-Core/DAC-rnaseq-pipeline/*. are managed using a Conda environment currently available on the Dartmouth HPC infrastructure (the Dartmouth Discovery cluster) but will be made more widely accessible in the near future. 

<img src="logo.jpg" width="250" height="140" >

## Pipeline summary:
The major steps implmented in the pipeline include: 

- FASTQ quality control assesment using [*FASTQC*](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- Read trimming for Poly-A tails, specified adapters, and read quality using [*cutadapt*](https://cutadapt.readthedocs.io/en/stable/)
- Alignment using [*STAR*](https://github.com/alexdobin/STAR)
- Quantification with [*HTSeq-count*](https://htseq.readthedocs.io/en/release_0.11.1/count.html) or [*RSEM*](https://deweylab.github.io/RSEM/)

As input, the pipeline takes raw data in FASTQ format, and produces quantified read counts (using *HTSeq-Count* or *RSEM*) as well as a detailed quality control report (including pre- and post-alignment QC metrics) for all processed samples. Quality control reports are aggregated into HTML files using *MultiQC*. 

R-code used to perform downstream exploratory data analysis and gene-level differential expression are also provided, however are currently detatched from the preprocessing and quality control steps and must be run separately. These scripts will be incorporated into the pipeline in the near future. 

## Implementation
The pipeline uses R-scripts to generate and submit jobs to the scheduler, and requires several variables to be defined by the user when running the pipeline: 

**Function arguments**
DAC_RNAseq_process (Lab, FastqRaw, SamNames, SeqMethod, AlignInd, AlignRef, PicardInt, PicardRef, QuantRef, CondaEnv, OutputFolder)   

* **Lab** - The name of the lab, or relevant project (used for file naming).
* **FastqRaw** - The absolute path to raw FASTQ files.
* **SamNames** - a vector of sample names (e.g. SamName <- c("SamName1", "SamName2", etc.)) that make up the prefixes of FASTQ file (e.g. 'SamName1' for 'SamName1_R1_001.fastq.gz'.
* **SeqMethod** - Either "fullLength" for assays profiling full transcripts or "3Prime" for 3'-end profiling assays. 
* **AlignInd** - Absolute path to the STAR index to be used as the reference genome. 
* **AlignRef** - Absolute path to the genome annotation (.gtffile ) to be used during alignment (to determine splice-site coordinates).
This is the reference that you would like to use during the alignment step, please give an absolute path (*.gtf).
* **PicardInt** - Absolute path to coordinates of ribosomal RNA sequences in reference genome, in [interval-list format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531852-Intervals-and-interval-lists).
* **PicardRef** - Absolute path to genome annotation in [RefFlat format](https://gatk.broadinstitute.org/hc/en-us/articles/360040509431-CollectRnaSeqMetrics-Picard-).
* **QuantRef** - Absolute path to genome annotation file (.gtf). 
* **CondaEnv**- This is the environment that includes all of the dependencies needed to run this pipeline, the yml file to create this environment is included in this directory (environment.yml).
* **OutputFolder** - Absolute path to directory for pipeline outputs.

**Argument definitions for example function call below**
* **Lab** <- SomeLab
* **FastqRaw** <- RNAseq_2-13-20/
* **SamNames** <- c("1_S1","2_S3","3_S5","4_S7","5_S")
* **SeqMethod** <- fullSized
* **AlignInd** <- genomic_references/human/STAR/hg38_index
* **AlignRef** <- genomic_references/human/ensembl-annotations/Homo_sapiens.GRCh38.97.gtf
* **PicardInt** <- genomic_references/human/CollectRnaSeqMetrics/Homo_sapiens.GRCh38.97.rRNA.interval_list
* **PicardRef** <- genomic_references/human/CollectRnaSeqMetrics/Homo_sapiens.GRCh38.97.refFlat.txt
* **QunatRef** <- genomic_references/human/RSEM/txdir/RSEMref
* **CondaEnv** <- conda activate conda_envs/rnaseq1
* **OutputFolder** <- ./

**Example function call from within an R terminal**   
source("rnaseq-pipeline.R")

DAC_RNAseq_process("SomeLab", "RNAseq_2-13-20/", c("1_S1","2_S3","3_S5","4_S7","5_S"), "fullSized", "genomic_references/human/STAR/hg38_index", "genomic_references/human/ensembl-annotations/Homo_sapiens.GRCh38.97.gtf", "genomic_references/human/CollectRnaSeqMetrics/Homo_sapiens.GRCh38.97.rRNA.interval_list", "genomic_references/human/CollectRnaSeqMetrics/Homo_sapiens.GRCh38.97.refFlat.txt", "genomic_references/human/RSEM/txdir/RSEMref", "conda activate conda_envs/rnaseq1", "./")


> **Contact & questions:** 
> Please address questions to *DataAnalyticsCore@groups.dartmouth.edu* or generate a issue in the GitHub repository. 

> **This pipeline was created with funds from the COBRE grant **1P20GM130454**. 
> If you use the pipeline in your own work, please acknowledge the pipeline by citing the grant number in your manuscript.**

