library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(data.table)


#-----------------------------------------------------------------------------
# merge KEGG Brite categories with PICRUSt2 output file and taxonomy from dada2 
#------------------------------------------------------------------------------

#load KEGG Birte categories 
brite_long <- read_csv("brite_database_long_clean.csv", show_col_types = FALSE)

head(brite_long)

#load piecrust results (with contribution) 

ko_contrib_df <- fread(
  "pred_metagenome_contrib.tsv",
  sep = "\t",
  data.table = FALSE)
 
head(ko_contrib_df) 


# Extract unique ASVs
unique_asvs <- unique(ko_contrib_df$taxon)

# Number of unique ASVs represented
length(unique_asvs)

# View them
sort(unique_asvs)



#load ASV taxonomy (ASV names must match the list in the pred_metagenome_contrib.tsv 
tax_df <- read_tsv(
  "ASV-taxonomy.txt",
  show_col_types = FALSE,
  progress = TRUE)
  
head(tax_df)
  
#clean the files 

ko_contrib_tax <- ko_contrib_df %>%
  left_join(tax_df, by = c("taxon" = "ID"))

head(ko_contrib_tax)

brite_map <- brite_long %>%
  distinct(KO, brite_A, brite_B, brite_C)
  
ko_contrib_tax2 <- ko_contrib_tax %>%
  mutate(KO = sub("^ko:", "", `function`))

  
ko_contrib_tax2 <- ko_contrib_tax2 %>%
  left_join(brite_map, by = "KO")
  
   
colnames(ko_contrib_tax2)


#load sample metadata 

metadata <- read_tsv( "metadata.txt",
  col_types = cols(
    samples = col_character(),
    species = col_character()))
	
head(metadata)

ko_contrib_tax2 <- ko_contrib_tax2 %>%
  left_join(metadata, by = c("sample" = "samples"))


write_csv(ko_contrib_tax2, "ko_contrib_tax_metadata.csv")


### -----------------------------------------------------------------
### Calculate percent contribution of of all genera to to the each brite_C category 
### -----------------------------------------------------------------


ko_contrib_tax3 <- fread(
  "ko_contrib_tax_metadata.csv",
  sep = ",",
  data.table = FALSE)
  
head(ko_contrib_tax3)


ko_contrib_tax_metab <- ko_contrib_tax3 %>%
  filter(brite_A == "Metabolism")



sum_briteC_by_genus <- ko_contrib_tax_metab %>%
  mutate(brite_C = if_else(is.na(brite_C) | brite_C == "", "Unclassified (no BRITE C)", brite_C)) %>%
  group_by(brite_C, Genus) %>%
  summarise(contrib = sum(norm_taxon_function_contrib, na.rm = TRUE), .groups = "drop") %>%
  group_by(brite_C) %>%
  mutate(percent = 100 * contrib / sum(contrib))



dim(sum_briteC_by_genus)

briteC_genus_mat_t <- sum_briteC_by_genus %>%
  filter(!is.na(brite_C)) %>%              # important: no empty/NA column names
  select(Genus, brite_C, percent) %>%
  pivot_wider(
    names_from  = brite_C,
    values_from = percent,
    values_fill = 0
  ) %>%
  column_to_rownames("Genus") %>%           # Genus = rownames
  t() %>%                                   # transpose
  as.data.frame() %>%
  rownames_to_column("brite_C")


write.csv(briteC_genus_mat_t,"TableS6-briteC_by_genus_transposed.csv", row.names = FALSE)



--------------------------------------------------------------------------------------------

briteC_KO_counts_dataset <- ko_contrib_tax_metab %>%
  filter(
    !is.na(brite_C),
    brite_C != ""
  ) %>%
  group_by(brite_B, brite_C) %>%
  summarise(
    n_KO_dataset = n_distinct(KO),
    .groups = "drop"
  ) %>%
  arrange(brite_B, brite_C)

briteC_KO_counts_dataset


write.csv(briteC_KO_counts_dataset,"TableS6-additional-column.csv", row.names = FALSE)




###------------ getting a list of all KOs for the top 29 genera with brite_C and brite_D categories,
## and KO names added at the end 


library(readr)
library(dplyr)
library(tibble)
library(tidyr)
library(data.table)

ko_contrib_tax2 <- fread(
  "ko_contrib_tax_metadata.csv",
  sep = ",",
  data.table = FALSE)
  
  
cat(colnames(ko_contrib_tax2), sep = "\n")
  
top29 <- c(
  "Acinetobacter",
  "Alphaproteobacteria UG",
  "Alteromonas",
  "Arenicella",
  "Brevinema",
  "Burkholderiales UG",
  "Cardiobacteriaceae UG",
  "Cetobacterium",
  "Cobetia",
  "Dokdonia",
  "Flavobacteriaceae UG",
  "Gammaproteobacteria UG",
  "Granulosicoccus",
  "Kordiimonas",
  "Lentisphaera",
  "Lewinella",
  "Maribacter",
  "Marivibrio",
  "Paracoccaceae UG",
  "Parendozoicomonas",
  "Phormidesmis",
  "Photobacterium",
  "Propionigenium",
  "Rubritalea",
  "Ruegeria",
  "Schizothrix",
  "Vibrio",
  "Vibrionaceae UG",
  "Woeseia"
)

ko_top29_metabolism <- ko_contrib_tax2 %>%
  filter(
    brite_A == "Metabolism",
    Genus %in% top29
  ) %>%
  select(
    Class,
    Order,
    Genus,
    KO,
    brite_B,
    brite_C
  ) %>%
  distinct() %>%
  arrange(Genus, brite_B, brite_C, KO)
  
  
 write.csv(
  ko_top29_metabolism,
  "Table_S7_top29_genera_unique_KOs_metabolism.csv",
  row.names = FALSE
) 

library(readr)
library(dplyr)
library(stringr)

# Read unformatted KEGG BRITE database
kegg_lines <- read_lines("kegg_brite_ko00001.txt")

# Extract KO number and KO name from D-level lines
ko_names <- tibble(line = kegg_lines) %>%
  filter(str_detect(line, "^D\\s+K\\d{5}")) %>%
  mutate(
    KO = str_extract(line, "K\\d{5}"),
    KO_name = str_replace(line, "^D\\s+K\\d{5}\\s+", "")
  ) %>%
  distinct(KO, KO_name)

# Add KO names to your top 29 genera table
ko_top29_metabolism_with_names <- ko_top29_metabolism %>%
  left_join(ko_names, by = "KO") %>%
  select(
    Class,
    Order,
    Genus,
    KO,
    KO_name,
    brite_B,
    brite_C
  ) %>%
  arrange(Genus, brite_B, brite_C, KO)

# Export
write.csv(
  ko_top29_metabolism_with_names,
  "Table_S7_top29_genera_unique_KOs_metabolism_with_KO_names.csv",
  row.names = FALSE
)



###---------------------------------------------------------------
# List of the 32 BRITE C categories to investigate
###---------------------------------------------------------------

briteC_32 <- c(
  "Betalain biosynthesis",
  "Biosynthesis of various plant secondary metabolites",
  "Flavone and flavonol biosynthesis",
  "Flavonoid biosynthesis",
  "Glucosinolate biosynthesis",
  "Indole alkaloid biosynthesis",
  "Isoflavonoid biosynthesis",
  "Isoquinoline alkaloid biosynthesis",
  "Phenylpropanoid biosynthesis",
  "Stilbenoid, diarylheptanoid and gingerol biosynthesis",
  "Tropane, piperidine and pyridine alkaloid biosynthesis",
  "Glycosaminoglycan biosynthesis - chondroitin sulfate / dermatan sulfate",
  "Glycosaminoglycan biosynthesis - heparan sulfate / heparin",
  "Glycosphingolipid biosynthesis - ganglio series",
  "Glycosphingolipid biosynthesis - globo and isoglobo series",
  "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis",
  "Mucin type O-glycan biosynthesis",
  "N-Glycan biosynthesis",
  "Other types of O-glycan biosynthesis",
  "Various types of N-glycan biosynthesis",
  "alpha-Linolenic acid metabolism",
  "Arachidonic acid metabolism",
  "Cutin, suberine and wax biosynthesis",
  "Linoleic acid metabolism",
  "Primary bile acid biosynthesis",
  "Sphingolipid metabolism",
  "Steroid biosynthesis",
  "Steroid hormone biosynthesis",
  "Retinol metabolism",
  "Insect hormone biosynthesis",
  "Monoterpenoid biosynthesis",
  "Zeatin biosynthesis"
)

KO_32_categories <- ko_contrib_tax2 %>%
  filter(brite_C %in% briteC_32) %>%
  select(
    brite_B,
    brite_C,
    KO
  ) %>%
  distinct() %>%
  arrange(brite_B, brite_C, KO)

KO_32_categories


# Add KO names to your top 29 genera table
KO_32_categories_with_names <- KO_32_categories %>%
  left_join(ko_names, by = "KO") %>%
  select(
    KO,
    KO_name,
    brite_B,
    brite_C
  ) %>%
  arrange(brite_B, brite_C, KO)


write.csv(
  KO_32_categories_with_names,
  "KO_32_categories_with_names.csv",
  row.names = FALSE
) 



### the below is an alternative for Table S6 -contribution to briteB 

library(data.table)
library(dplyr)
library(tidyr)
library(tibble)

### -----------------------------------------------------------------
### Calculate percent contribution of all genera to each BRITE B category
### -----------------------------------------------------------------

ko_contrib_tax <- fread(
  "ko_contrib_tax_metadata.csv",
  sep = ",",
  data.table = FALSE
)

head(ko_contrib_tax)

# Keep only KEGG BRITE A = Metabolism
ko_contrib_tax_metab <- ko_contrib_tax %>%
  filter(brite_A == "Metabolism")

# Sum normalized taxon-function contributions by BRITE B category and genus
sum_briteB_by_genus <- ko_contrib_tax_metab %>%
  mutate(
    brite_B = if_else(
      is.na(brite_B) | brite_B == "",
      "Unclassified (no BRITE B)",
      brite_B
    )
  ) %>%
  group_by(brite_B, Genus) %>%
  summarise(
    contrib = sum(norm_taxon_function_contrib, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(brite_B) %>%
  mutate(
    percent = 100 * contrib / sum(contrib)
  ) %>%
  ungroup()

dim(sum_briteB_by_genus)

# Convert to BRITE B x genus matrix
briteB_genus_mat_t <- sum_briteB_by_genus %>%
  select(Genus, brite_B, percent) %>%
  pivot_wider(
    names_from = brite_B,
    values_from = percent,
    values_fill = 0
  ) %>%
  column_to_rownames("Genus") %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("brite_B")

# Export
write.csv(
  briteB_genus_mat_t,
  "TableS6-briteB_by_genus_transposed.csv",
  row.names = FALSE
)

------------------------------------------------------------------------------------
library(readr)
library(dplyr)
library(tidyr)
library(tibble)

# Complete KEGG BRITE mapping database
brite_long <- read_csv(
  "brite_database_long_clean.csv",
  show_col_types = FALSE
)

# New PICRUSt2 unstratified KO predictions
ko_unstrat <- read_tsv(
  "pred_metagenome_unstrat.tsv",
  show_col_types = FALSE
)

dim(ko_unstrat)
head(ko_unstrat)
colnames(ko_unstrat)
nrow(ko_unstrat)


ko_unstrat2 <- ko_unstrat %>%
  rename(KO = `function`) %>%
  mutate(
    KO = sub("^ko:", "", KO)
  )
  
head(ko_unstrat2$KO)
n_distinct(ko_unstrat2$KO)

brite_map <- brite_long %>%
  distinct(
    KO,
    brite_A,
    brite_B,
    brite_C
  )

ko_unstrat_brite <- ko_unstrat2 %>%
  left_join(
    brite_map,
    by = "KO"
  )
  
ko_unstrat_metabolism <- ko_unstrat_brite %>%
  filter(
    brite_A == "Metabolism",
    !is.na(brite_C),
    brite_C != ""
  )
 

 
tableS6 <- read_csv(
  "TableS6-additional-column.csv",
  show_col_types = FALSE
)

briteC_contrib <- sort(unique(tableS6$brite_C))
briteC_unstrat <- sort(unique(ko_unstrat_metabolism$brite_C))

# Number in each dataset
length(briteC_contrib)
length(briteC_unstrat)

# Categories found only in contribution output
setdiff(briteC_contrib, briteC_unstrat)

# Categories found only in unstratified output
setdiff(briteC_unstrat, briteC_contrib) 

library(dplyr)
library(tibble)

# Start from KO predictions already merged to BRITE hierarchy
# and restricted to brite_A == "Metabolism"

sample_cols <- colnames(ko_unstrat2)[
  colnames(ko_unstrat2) != "KO"
]

briteC_abund_df <- ko_unstrat_metabolism %>%
  group_by(brite_C) %>%
  summarise(
    across(
      all_of(sample_cols),
      ~ sum(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  )

dim(briteC_abund_df)


briteC_rel <- briteC_abund_df

briteC_rel[sample_cols] <- sweep(
  briteC_rel[sample_cols],
  2,
  colSums(briteC_rel[sample_cols]),
  FUN = "/"
) * 100

round(
  colSums(briteC_rel[sample_cols]),
  6
)

briteC_long <- briteC_rel %>%
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "relative_abundance"
  ) %>%
  left_join(
    metadata,
    by = c("sample" = "samples")
  )
  
  
  
table(briteC_long$species)


library(dplyr)
library(tidyr)
library(purrr)

# ------------------------------------------------------------
# Kruskal-Wallis test for each BRITE C category
# ------------------------------------------------------------

kw_results <- briteC_long %>%
  group_by(brite_C) %>%
  summarise(
    H = unname(
      kruskal.test(relative_abundance ~ species)$statistic
    ),
    df = unname(
      kruskal.test(relative_abundance ~ species)$parameter
    ),
    p_value = kruskal.test(
      relative_abundance ~ species
    )$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH")
  ) %>%
  arrange(p_adj_BH)

kw_results

sum(kw_results$p_adj_BH < 0.05)

kw_significant <- kw_results %>%
  filter(p_adj_BH < 0.05)

kw_significant

write.csv(
  kw_results,
  "TableS8_BRITE_C_Kruskal_Wallis_all_results.csv",
  row.names = FALSE
)

library(vegan)

# BRITE C relative-representation matrix:
# currently categories are rows and samples are columns

briteC_sample_mat <- briteC_rel %>%
  column_to_rownames("brite_C") %>%
  t() %>%
  as.data.frame()

dim(briteC_sample_mat)


metadata_func <- metadata %>%
  filter(samples %in% rownames(briteC_sample_mat)) %>%
  slice(match(rownames(briteC_sample_mat), samples))

# Critical checks
all(rownames(briteC_sample_mat) == metadata_func$samples)

table(metadata_func$species)


briteC_bray <- vegdist(
  briteC_sample_mat,
  method = "bray"
)

set.seed(1234)

func_permanova <- adonis2(
  briteC_bray ~ species,
  data = metadata_func,
  permutations = 999
)

func_permanova

func_bd <- betadisper(
  briteC_bray,
  group = metadata_func$species
)

permutest(
  func_bd,
  permutations = 999
)


