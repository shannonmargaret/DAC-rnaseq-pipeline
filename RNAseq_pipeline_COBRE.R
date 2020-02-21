#-------------------------------------#
# [1] RNA-seq pipeline for COBRE #
# From fastq to raw counts           #
#-------------------------------------#

rm(list = ls())
library("glue")
whoseData <- "GaurLab"
myRaw <- "/dartfs-hpc/rc/lab/G/GSR_Active/Labs/Gaur/RNAseq_2-13-20/"
myQC <- "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/rnaseq_pipeline/shannon/fastqc/"
myTri <- "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/rnaseq_pipeline/shannon/trim/"
myAli <- "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/rnaseq_pipeline/shannon/alignment/"
myExp <- "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/rnaseq_pipeline/shannon/rawcounts/"
myTmp <- "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/rnaseq_pipeline/shannon/tmp/"
#--
StarInd <- "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/genomic_references/human/STAR/hg38_index"
StarRef <- "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/genomic_references/human/ensembl-annotations/Homo_sapiens.GRCh38.97.gtf"
PicardRef <- "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/genomic_references/human/CollectRnaSeqMetrics/Homo_sapiens.GRCh38.97.refFlat.txt"
PicardInt <- "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/genomic_references/human/CollectRnaSeqMetrics/Homo_sapiens.GRCh38.97.rRNA.interval_list"
RsemRef <- "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/genomic_references/human/RSEM/txdir/RSEMref"
#--
mytar <- c("1_S1","10_S4","11_S6","12_S8","13_S10","14_S12","15_S14","16_S16","17_S17","18_S18","19_S19","2_S3","20_S20","3_S5","4_S7","5_S9","6_S11","7_S13","8_S15","9_S2","Undetermined_S0")
#-------------------------
for(i in 1:length(mytar)){
	#-------
	myfiles <- dir(myRaw, pattern = mytar[i])
	
	sample_id_1 <- unique(sapply(strsplit(myfiles[1], "_"), "[", 1))
	sample_id_2 <- unique(sapply(strsplit(myfiles[1], "_"), "[", 2))
	sample_id <- paste(sample_id_1, sample_id_2, sep = "_")
	
	print (myfiles)
	#var1<- glue(mytar[i], "_R1_")
	#var2<- glue(mytar[i], "_R2_")
	my_R1 <- myfiles[grep("_R1_", myfiles)]
	my_R2 <- myfiles[grep("_R2", myfiles)]

	myoutf1 <- glue(myTmp, whoseData, "_", mytar[i], ".sp")
	conOut <- file(myoutf1, "w")
	tmp <- "#!/bin/bash"
	writeLines(tmp, conOut)
	writeLines("\n", conOut)
	#-----------------
	# Quality control
	#-----------------
	writeLines("source /optnfs/common/miniconda3/etc/profile.d/conda.sh", conOut)
#	writeLines("conda activate anaconda3", conOut)
	writeLines("conda activate /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/conda_envs/rnaseq1" ,conOut)
	writeLines("echo 'QC step'", conOut)
	writeLines("you made it to line 43")
	comm <- "fastqc myFASTQ --outdir="
	A <- length(myRaw)
	B <- length(my_R1)
	print("length of my_R1")
	print(my_R1)
	print(B)
	#name <- glue(myRaw, my_R1)
	tmp <- gsub("myFASTQ", glue(myRaw, my_R1), comm)
	writeLines("you made it to line 46")
	tmp <- glue(tmp, myQC)
	writeLines("you made it to line 48")
	writeLines(tmp, conOut)
	writeLines("\n", conOut)
	#--
	tmp <- gsub("myFASTQ", glue(myRaw, my_R2), comm)
	tmp <- glue(tmp, myQC)
	writeLines(tmp, conOut)
	writeLines("\n", conOut)
	#------
	# Trim
	#------
	writeLines("echo 'Trim step'", conOut)
	comm_1 <- "cutadapt -A 'A{76}' -g 'T{76}'"
	comm_2 <- "-o R1_fastq_out -p R2_fastq_out R1_fastq_in R2_fastq_in"
	comm_3 <- "-m 5 --nextseq-trim=20 > Sample_ID.report"  #-j 8 need python 3

	my_R1_in <- paste(myRaw, sample_id, "_R1_001.fastq.gz", sep = "")
	my_R2_in <- paste(myRaw, sample_id, "_R2_001.fastq.gz", sep = "")

	my_R1_out <- paste(myTri, sample_id, "_R1_001.fastq.gz", sep = "")
	my_R2_out <- paste(myTri, sample_id, "_R2_001.fastq.gz", sep = "")

	comm_2 <- gsub("R1_fastq_out", my_R1_out, comm_2)
	comm_2 <- gsub("R2_fastq_out", my_R2_out, comm_2)
	comm_2 <- gsub("R1_fastq_in", my_R1_in, comm_2)
	comm_2 <- gsub("R2_fastq_in", my_R2_in, comm_2)

	tag <- paste(myTri, sample_id, sep = "")
	comm_3 <- gsub("Sample_ID", tag, comm_3)

	comm <- paste(comm_1, comm_2, comm_3, sep = " ")
	tmp <- comm
	writeLines(tmp, conOut)
	writeLines("\n", conOut)
	#----------
	# Alignment
	#----------
	writeLines("echo 'Alignment step'", conOut)

	comm_1 <- "STAR --quantMode TranscriptomeSAM --genomeDir myind --sjdbGTFfile mygene"
	comm_2 <- "--runThreadN 4 --twopassMode Basic --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --outFilterMultimapNmax 20 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1"
	comm_3 <- "--readFilesIn trimed_R1_fastq trimed_R2_fastq --readFilesCommand zcat --outFileNamePrefix Sample_ID."

	comm_1 <- gsub("myind", StarInd, comm_1)
	comm_1 <- gsub("mygene", StarRef, comm_1)

	myTrim_R1 <- paste(myTri, sample_id, "_R1_001.fastq.gz", sep = "")
	myTrim_R2 <- paste(myTri, sample_id, "_R2_001.fastq.gz", sep = "")
	tag <- paste(myAli, sample_id, sep = "")

	comm_3 <- gsub("trimed_R1_fastq", myTrim_R1, comm_3)
	comm_3 <- gsub("trimed_R2_fastq", myTrim_R2, comm_3)
	comm_3 <- gsub("Sample_ID", tag, comm_3)

	comm <- paste(comm_1, comm_2, comm_3, sep = " ")
	tmp <- comm
	writeLines(tmp, conOut)
	writeLines("\n", conOut)
	#------------------
	# QC for BAM file
	#------------------
	writeLines("echo 'Alignment QC step'", conOut)
	comm <- "picard CollectRnaSeqMetrics I=myInBam O=myOut REF_FLAT=myTxt STRAND=SECOND_READ_TRANSCRIPTION_STRAND RIBOSOMAL_INTERVALS=interval_list MAX_RECORDS_IN_RAM=1000000"
	myBam_in <- paste(myAli, sample_id, ".Aligned.sortedByCoord.out.bam", sep = "")
	my_out <- paste(myAli, sample_id, ".RNA_Metrics", sep = "")
	comm <- gsub("myInBam", myBam_in, comm)
	comm <- gsub("myOut", my_out, comm)
	comm <- gsub("myTxt", PicardRef, comm)
	comm <- gsub("interval_list", PicardInt, comm)
	tmp <- comm
	writeLines(tmp, conOut)
	writeLines("\n", conOut)
	#--
	comm <- "picard MarkDuplicates I=myInBam O=myOutBam M=myTxt OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 CREATE_INDEX=false"
	myBam_in <- paste(myAli, sample_id, ".Aligned.sortedByCoord.out.bam", sep = "")
	myBam_out <- paste(myAli, sample_id, ".marked_duplicates.bam", sep = "")
	mytxt <- paste(myAli, sample_id, ".marked_dup_metrics.txt", sep = "")
	comm <- gsub("myInBam", myBam_in, comm)
	comm <- gsub("myOutBam", myBam_out, comm)
	comm <- gsub("myTxt", mytxt, comm)
	tmp <- comm
	writeLines(tmp, conOut)
	writeLines("\n", conOut)
	#-----------
	# Expression
	#-----------
	writeLines("echo 'RSEM step'", conOut)
	myIn <- paste(myAli, sample_id, ".Aligned.toTranscriptome.out.bam", sep = "")
	tag <- paste(myExp, sample_id, sep = "")

	comm <- "rsem-calculate-expression --paired-end --alignments --strandedness reverse -p 8 input_bam mygene Sample_ID"

	comm <- gsub("input_bam", myIn, comm)
	comm <- gsub("mygene", RsemRef, comm)
	comm <- gsub("Sample_ID", tag, comm)
	
	tmp <- comm
	writeLines(tmp, conOut)
	writeLines("conda deactivate", conOut)
	close(conOut)
}

mysub <- paste(myTmp, "submit.sp", sep="")
conOut <- file(mysub, "w")
curLine <- paste("mksub -l nodes=1:ppn=20 -l walltime=20:0:0 ", whoseData, "_", mytar, ".sp", sep = "")
writeLines(curLine, conOut)
close(conOut)

comm <- paste("chmod u+x ", mysub, sep="")
system(comm)
"bash submit.sp"


###################
# [2] Multiqc for every step #
###################

rm(list = ls())

whoseData <- "GaurLab"
myQC <- "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/rnaseq_pipeline/shannon/fastqc/"
myTri <- "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/rnaseq_pipeline/shannon/trim/"
myAli <- "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/rnaseq_pipeline/shannon/alignment/"
myExp <- "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/rnaseq_pipeline/shannon/rawcounts/"
myTmp <- "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/rnaseq_pipeline/shannon/tmp/"
myMul <- "/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/rnaseq_pipeline/shannon/multiqc/"
myoutf <- paste(myTmp, whoseData, "_MultiQC.sp", sep = "")
comm <- "multiqc fastqc trim bam -n myfolder"

conOut <- file(myoutf, "w")
tmp <- "#!/bin/bash"
writeLines(tmp, conOut)
writeLines("\n", conOut)
tmp <- "source activate QCpipe"
writeLines(tmp, conOut)
writeLines("\n", conOut)
tmp <- gsub("fastqc", myQC, comm)
tmp <- gsub("trim", myTri, tmp)
tmp <- gsub("bam", myAli, tmp)
tmp <- gsub("myfolder", paste(myMul, whoseData, "_multiqc_report.html", sep = ""), tmp)
writeLines(tmp, conOut)
writeLines("\n", conOut)
tmp <- "source deactivate"
writeLines(tmp, conOut)
writeLines("\n", conOut)
close(conOut)
