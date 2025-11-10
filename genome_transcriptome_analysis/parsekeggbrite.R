# packages
library(httr)
library(readr)
library(dplyr)
library(stringr)

# Lines 9 to 40 written with assistance from Claude AI

# Flatten BRITE hierarchy to be able to join columns from each genome

url <- "https://rest.kegg.jp/get/br:ko00001"
lines <- content(GET(url), "text")
lines <- strsplit(lines, "\n")[[1]]

level_a <- level_b <- level_c <- NA
records <- list()

for (line in lines) {
  if (grepl("^A ", line)) {
    level_a <- sub("^A\\s+", "", line)
  } else if (grepl("^B ", line)) {
    level_b <- sub("^B\\s+", "", line)
  } else if (grepl("^C ", line)) {
    level_c <- sub("^C\\s+", "", line)
  } else if (grepl("^D ", line)) {
    parts <- strsplit(sub("^D\\s+", "", line), " ", fixed = TRUE)[[1]]
    ko_id <- parts[1]
    desc <- paste(parts[-1], collapse = " ")
    records <- append(records, list(data.frame(
      Level_A = level_a,
      Level_B = level_b,
      Level_C = level_c,
      KO_ID = ko_id,
      Description = desc,
      stringsAsFactors = FALSE
    )))
  }
}

brite_df <- bind_rows(records)

# load kegg annotations for first I. cascadensis genome - prefiltered for genes with kegg hits
icas1_kegg = read.csv("genome_transcriptome_analysis/icas1_keggmatches.csv")
# count number of genes for each KO_ID
icas1_kegg_counts = icas1_kegg %>% group_by(KO_ID) %>% mutate(icas1 = n())
# remove gene ID
icas1_kegg_counts = icas1_kegg_counts %>% subset(select = -c(gene))
# remove duplicate rows
icas1_kegg_counts = icas1_kegg_counts %>% unique()
# join counts from bc2001 with kegg brite annotations
brite_df_5genomes = full_join(brite_df, icas1_kegg_counts)

# do the same for second I. cascadensis genome 
icas2_kegg = read.csv("genome_transcriptome_analysis/icas2_keggmatches.csv")
icas2_kegg_counts = icas2_kegg %>% group_by(KO_ID) %>% mutate(icas2 = n())
icas2_kegg_counts = icas2_kegg_counts %>% subset(select = -c(gene))
icas2_kegg_counts = icas2_kegg_counts %>% unique()
brite_df_5genomes = full_join(brite_df_5genomes, icas2_kegg_counts)

# do the same for V. vermiformis
vverm_kegg = read.csv("genome_transcriptome_analysis/vverm_keggmatches.csv")
vverm_kegg_counts = vverm_kegg %>% group_by(KO_ID) %>% mutate(vverm = n())
vverm_kegg_counts = vverm_kegg_counts %>% subset(select = -c(gene))
vverm_kegg_counts = vverm_kegg_counts %>% unique()
brite_df_5genomes = full_join(brite_df_5genomes, vverm_kegg_counts)

# do the same for A. terricola
aterricola_kegg = read.csv("genome_transcriptome_analysis/aterricola_keggmatches.csv")
aterricola_kegg_counts = aterricola_kegg %>% group_by(KO_ID) %>% mutate(aterricola = n())
aterricola_kegg_counts = aterricola_kegg_counts %>% subset(select = -c(gene))
aterricola_kegg_counts = aterricola_kegg_counts %>% unique()
brite_df_5genomes = full_join(brite_df_5genomes, aterricola_kegg_counts)

# do the same for D. discoideum
dicty_kegg = read.csv("genome_transcriptome_analysis/dicty_keggmatches.csv")
dicty_kegg_counts = dicty_kegg %>% group_by(KO_ID) %>% mutate(dicty = n())
dicty_kegg_counts = dicty_kegg_counts %>% subset(select = -c(gene))
dicty_kegg_counts = dicty_kegg_counts %>% unique()
brite_df_5genomes = full_join(brite_df_5genomes, dicty_kegg_counts)

#now have all 5 genomes for comparison with cou nts to each kegg function
write.csv(brite_df, "genome_transcriptome_analysis/keggbrite_5genomes.csv")

# format for transcriptome analysis - ignoring splice variants
icas1_kegg_nosplice <- read.csv("genome_transcriptome_analysis/icas1_kegg_nosplice.csv")
icas1_counts <- read.csv("genome_transcriptome_analysis/counts_icas1.csv")
icas2_kegg_nosplice <- read.csv("genome_transcriptome_analysis/icas2_kegg_nosplice.csv")
icas2_counts <- read.csv("genome_transcriptome_analysis/counts_icas2.csv")

# collapse genes with multiple duplicated variants
icas1_kegg_nosplice = unique(icas1_kegg_nosplice)
icas2_kegg_nosplice = unique(icas2_kegg_nosplice)

# checking how many genes have variants with different assignments
icas1_variants = icas1_kegg_nosplice %>% group_by(Geneid) %>% mutate(variants = n())
icas2_variants = icas2_kegg_nosplice %>% group_by(Geneid) %>% mutate(variants = n())

# Join brite annotations by KO_ID
icas1_kegg_annot = full_join(brite_df, icas1_kegg_nosplice)
icas2_kegg_annot = full_join(brite_df, icas2_kegg_nosplice)

# Join with transcript counts by Geneid
icas1_annotcounts = full_join(icas1_kegg_annot, icas1_counts)
icas2_annotcounts = full_join(icas2_kegg_annot, icas2_counts)

# Remove columns to remove duplicates
icas1_annottrim = subset(icas1_annotcounts, select = c(KO_ID, Description, Geneid, Chr, Start, End, Strand, Length, 
                                                         mapped_to_Icas1))
icas1_annottrim = unique(icas1_annottrim)
icas2_annottrim = subset(icas2_annotcounts, select = c(KO_ID, Description, Geneid, Chr, Start, End, Strand, Length, 
                                                        mapped_to_Icas2))
icas2_annottrim = unique(icas2_annottrim)

icas1_annottrim["Sample"] <- "Icas1"
icas1_annottrim = icas1_annottrim %>% subset(select = c(KO_ID, Description, Sample, mapped_to_Icas1))
icas1_annottrim = rename(icas1_annottrim, mapped_to_Icas = mapped_to_Icas1)

icas2_annottrim["Sample"] <- "Icas2"
icas2_annottrim = icas2_annottrim %>% subset(select = c(KO_ID, Description, Sample, mapped_to_Icas2))
icas2_annottrim = rename(icas2_annottrim, mapped_to_Icas = mapped_to_Icas2)

combined_annot = full_join(icas1_annottrim, icas2_annottrim)

sorted_combined_annot = combined_annot %>%
  arrange(Description, mapped_to_Icas)

# remove functions with no hits
noNA_combined_annottrim = sorted_combined_annot %>% na.omit() %>%
  unique()

# write individual mapped counts files
write.csv(icas1_annotcounts, "genome_transcriptome_analysis/icas1_annotcounts.csv")                          
write.csv(icas2_annotcounts, "genome_transcriptome_analysis/icas2_annotcounts.csv")
# write combined annotation count file
write.csv(combined_annot, "genome_transcriptome_analysis/combined_rnaannot_full.csv")
