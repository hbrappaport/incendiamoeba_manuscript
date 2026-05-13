Folders:
1. genome_stats
 - genomestats_compare.R - plotting script for Fig. 4A, associated files: busco.csv, genomestats.csv, orthogroup_comparison.csv

2. genome_kegg_functions
- parsekeggbrite.R - 
    assign functional annotations to KEGG assignments for functional enrichment, associated files: keggmatches.csv for each of 5 genomes for comparison, output file: keggbrite_5genomes.csv

3. RNAseq_experiment
- processing_rnaseq_reads.R -
   code for processing fastq rnaseq data from temperature experiment (48C v 61C), performing deseq2 analysis, and calculating transcripts-per-million (TPM) normalized counts.
- bc2001_bin7_mapped.txt : per-feature tabulated raw counts mapped to the Incendi pacbio bin across each temperature replicate
- sample_metadata.csv: metadata associated with sample IDs in the columns of bc2001_bin7_mapped.txt . Used for deseq2 fold-change analysis between 48C and 61C
  
4. orthofinder
- output files from Orthofinder including Orthogroups.GeneCount.csv and Statistics_PerSpecies.csv
