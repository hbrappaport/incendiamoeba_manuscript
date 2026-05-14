#Performed read mapping of each rnaseq fastq file against the Incendi Genome using bbmap  - bbmap.sh in=*_R1_001.fastq in2=*_R2_001.fastq out=*_mapped.bam ref=/home/ngilbe02/data/pacbio/incendi/incendi_genomes_metabat2_fna/bins/bin.7.fna covstats=001_covstats.txt 

#tabulated counts with featurecounts using the braker gtf feature table- featureCounts -p -t gene -g Geneid -F GTF -a /home/ngilbe02/data/rnaseq/incendi_pilot_temp_data/bbmap_to_bin7_unzipped/featurecounts/bc2001_braker.gtf -o *_mapped.txt *.bam

#Packages
library( "DESeq2" )
library(ggplot2)
library(dplyr)
library(tidyr)


#import mapped reads
rawcts <- read.table("genome_transcriptome_analysis/RNAseq_experiment/bc2001_bin7_mapped.txt", sep = "\t", header = T)


###fix column headers

colnames(rawcts) <- gsub("X.home.ngilbe02.data.rnaseq.incendi_pilot_temp_data.bbmap_to_bin7_unzipped.", "", colnames(rawcts))
colnames(rawcts) <- gsub("_mapped.bam", "", colnames(rawcts))

#metadata
metadata <- read.csv("genome_transcriptome_analysis/RNAseq_experiment/sample_metadata.csv")


#Construct deseq2dataset object
dds <- DESeqDataSetFromMatrix(countData=rawcts[c(1, 7:14)], 
                              colData=metadata, 
                              design=~temperature, tidy = TRUE)

dds <- DESeq(dds)

res <- results(dds)
head(res)

summary(res)

#Sort list by p-value
  res <- res[order(res$padj),]
  head(res)
  
  results_pilot <- as.data.frame(res)
  
#Get normalized VST transcript values
  vsdata_pilot <- as.data.frame(assay(vsdata))
  vsdata_pilot$Geneid <- rownames(vsdata_pilot)  

#Merge w/deseq log2fc data
  results_pilot$Geneid <- rownames(results_pilot)
  
  results_pilot <- merge(results_pilot, vsdata_pilot, by = "Geneid")

####Perform transcripts per million (TPM) normalization
   # counts: matrix of RAW  read counts (genes x samples)
   # gene_lengths: vector of gene lengths in base pairs
   rownames(rawcts) <- rawcts$Geneid
   
   rpk <- rawcts[c(7:14)] / (rawcts$Length / 1000)
   tpms <- t(t(rpk) / colSums(rpk) * 1e6)
   
  #also calculate averaged tpm for 61 degree C condition 
  tpms <- tpms %>%
    mutate(ave_tpms_61 = ("8123_002"+"8123_004"+"8123_006"+"8123_008")/4)

  write.csv(tpms ,"pilot_tpms.csv")  
  
 