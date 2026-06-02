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


#Gammaproteobacteria UF had a typo should be UG 
ko_contrib_tax2 <- ko_contrib_tax2 %>%
  mutate(
    Genus = recode(
      Genus,
      "Gammaproteobacteria UF" = "Gammaproteobacteria UG"))

write_csv(ko_contrib_tax2, "ko_contrib_tax_metadata.csv")

###--------------------------------------------------------------------------------------
### Rename some briteC categories that are not found in bacteria to something more general
###---------------------------------------------------------------------------------------

ko_contrib_tax2 <- fread(
  "ko_contrib_tax_metadata.csv",
  sep = ",",
  data.table = FALSE)


ko_contrib_tax2 %>%
  count(brite_C == "", is.na(brite_C))


#need to rename problematic brite_C (not found in bacteria so renamed to something more general) 

briteC_rename <- fread(
  "briteC-to-rename.txt",
  sep = "\t",
  data.table = FALSE)

# Rename only the categories listed in the map; leave everything else unchanged
ko_contrib_tax2_fixed <- ko_contrib_tax2 %>%
  left_join(briteC_rename, by = "brite_C") %>%
  mutate(
    brite_C = if_else(!is.na(brite_C_new), brite_C_new, brite_C)
  ) %>%
  select(-brite_C_new)


#double check: Which BRITE C values got changed?

present_mapped <- ko_contrib_tax2 %>%
  distinct(brite_C) %>%
  inner_join(briteC_rename, by="brite_C")

present_mapped


fwrite(ko_contrib_tax2_fixed, "ko_contrib_tax_metadata_briteC_renamed.csv")




### -----------------------------------------------------------------
### Calculate percent contribution of of all genera to to the each brite_C category 
### -----------------------------------------------------------------


ko_contrib_tax3 <- fread(
  "ko_contrib_tax_metadata_briteC_renamed.csv",
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


write.csv(briteC_genus_mat_t,"TableS5-briteC_by_genus_transposed.csv", row.names = FALSE)


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
  "Table_S6_top29_genera_unique_KOs_metabolism.csv",
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
  "Table_S6_top29_genera_unique_KOs_metabolism_with_KO_names.csv",
  row.names = FALSE
)