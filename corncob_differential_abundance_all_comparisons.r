library(tidyverse)
library(vegan)
library(ANCOMBC)
library(phyloseq)
library(RColorBrewer)
library(randomcoloR)
library(microViz)
library(microbiome)
library(cowplot)
library(grid)
library(scales)
library(corncob)

#read-in data 

otu <- read.table("final_ASV_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("final_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples <-read.table("metadata.txt",sep="\t",header=TRUE,row.names=1)


taxon<-as.matrix(taxon)
TAX = tax_table(taxon)


ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(TAX))


##-------------------------------------------------------
#testing for differential abundance across enviornments 
#formula is the variable being tested for abundance 
#phi.formula is effect on dispersion 
#formula.null is the effect on abundance 
#phi.formula_null 
##-------------------------------------------------------

#diff_abudance class level 
ps_class <- ps %>% tax_glom("Class")

sample_data(ps_class)$environment <- factor(
  sample_data(ps_class)$environment,
  levels = c("fishskin", "substrate", "seawater")
)
da_class_en <- differentialTest(formula = ~ environment,
                                 phi.formula = ~ environment,
                                 formula_null = ~ 1,
                                 phi.formula_null = ~ environment,
                                 data = ps_class,
                                 test = "Wald", boot = FALSE,
                                 fdr_cutoff = 0.05)					
								 
plot(da_class_en,level="Class")+labs(x = "Coefficient and Confidence Interval", y = "Classes")+ theme(panel.grid = element_blank ())
ggsave("FigureS2-enviro-classes.tiff", units="in", width=10, height=10, dpi=300, compression = 'lzw')
ggsave("FigureS2-enviro-classes.pdf", units="in", width=10, height=10, dpi=300)


# taxonomy table from the class-level phyloseq object
tax_class <- as.data.frame(tax_table(ps_class))

# significant taxa IDs
sig_ids <- da_class_en$significant_taxa

# get corresponding Class names
sig_classes <- tax_class[sig_ids, "Class", drop = FALSE]

sig_classes
sig_class_table <- data.frame(
  taxon_id = sig_ids,
  Class = tax_class[sig_ids, "Class"]
)

da_class_en_FDR <- as.data.frame(da_class_en$p_fdr)

da_class_en_FDR$taxon_id <- rownames(da_class_en_FDR)
da_class_en_FDR$Class <- tax_class[rownames(da_class_en_FDR), "Class"]

# optional: put Class first
da_class_en_FDR <- da_class_en_FDR[, c("taxon_id", "Class",
                                       setdiff(names(da_class_en_FDR),
                                               c("taxon_id", "Class")))]

write.table(da_class_en_FDR, "da_class_en_FDR_with_Class.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)


#addting information on CI and interpretation 


library(dplyr)

# -----------------------------
# 1. Extract coefficients
# -----------------------------

sig_ids <- da_class_en$significant_taxa

coef_sig_class <- do.call(rbind, lapply(seq_along(sig_ids), function(i) {
  
  taxon_id_i <- sig_ids[i]
  mod <- da_class_en$significant_models[[i]]
  
  coefs <- coef(mod)
  
  data.frame(
    taxon_id = taxon_id_i,
    term = rownames(coefs),
    coefficient = coefs[, 1],
    row.names = NULL
  )
}))

# Keep mu coefficients only
coef_sig_class_mu <- coef_sig_class %>%
  filter(grepl("^mu", term))

# -----------------------------
# 2. Add taxonomy and FDR
# -----------------------------

tax_class <- as.data.frame(tax_table(ps_class))

fdr_class <- as.data.frame(da_class_en$p_fdr)
fdr_class$taxon_id <- rownames(fdr_class)
fdr_class$Class <- tax_class[rownames(fdr_class), "Class"]

fdr_class_sig <- fdr_class %>%
  filter(taxon_id %in% sig_ids)

class_table <- coef_sig_class_mu %>%
  left_join(fdr_class_sig, by = "taxon_id")

# -----------------------------
# 3. Clean comparison names
# -----------------------------

class_table <- class_table %>%
  mutate(
    comparison = term,
    comparison = gsub("mu\\.", "", comparison),
    comparison = gsub("environment", "", comparison)
  ) %>%
  filter(comparison != "(Intercept)")

# -----------------------------
# 4. Add enriched_in column
# assumes fishskin is reference
# -----------------------------

class_table <- class_table %>%
  mutate(
    enriched_in = case_when(
      comparison == "seawater" & coefficient > 0 ~ "seawater",
      comparison == "seawater" & coefficient < 0 ~ "fish skin",
      comparison == "substrate" & coefficient > 0 ~ "substrate",
      comparison == "substrate" & coefficient < 0 ~ "fish skin",
      TRUE ~ NA_character_
    )
  )

# -----------------------------
# 5. Extract CI from plot object
# -----------------------------

p_class <- plot(da_class_en, level = "Class")
plot_data_class <- p_class$data

ci_class_clean <- plot_data_class %>%
  mutate(
    Class = taxa,
    comparison = gsub("\nDifferential Abundance", "", variable),
    comparison = gsub("environment", "", comparison),
    coefficient_plot = x,
    CI_excludes_zero = (xmin <= 0 & xmax <= 0) | (xmin >= 0 & xmax >= 0),
    CI_crosses_zero = !CI_excludes_zero
  ) %>%
  select(Class, comparison, coefficient_plot, xmin, xmax,
         CI_excludes_zero, CI_crosses_zero)

# -----------------------------
# 6. Merge CI with class table
# -----------------------------

TableS2_DA_classes <- class_table %>%
  left_join(
    ci_class_clean,
    by = c("Class", "comparison")
  ) %>%
  select(
    Class, comparison, coefficient, xmin, xmax,
    CI_excludes_zero, CI_crosses_zero,
    enriched_in, taxon_id, term, everything()
  ) %>%
  arrange(Class, comparison)

# -----------------------------
# 7. Save Table S2
# -----------------------------

write.table(
  TableS2_DA_classes,
  "TableS2_DA_classes_environment_FDR_mu_CI.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


#----------------------------------------
#differential abundance , order level 
#------------------------------------------										 		
ps_order <- ps %>% tax_glom("Order")

sample_data(ps_order)$environment <- factor(
  sample_data(ps_order)$environment,
  levels = c("fishskin", "substrate", "seawater")
)

da_order_en <- differentialTest(
  formula = ~ environment,
  phi.formula = ~ environment,
  formula_null = ~ 1,
  phi.formula_null = ~ environment,
  data = ps_order,
  test = "Wald",
  boot = FALSE,
  fdr_cutoff = 0.05
)

plot(da_order_en,level="Order")+labs(x = "Coefficient and Confidence Interval", y = "Order")+ theme(panel.grid = element_blank ())
ggsave("FigureS3-enviro-orders.tiff", units="in", width=12, height=12, dpi=300, compression = 'lzw')
ggsave("FigureS3-enviro-orders.pdf", units="in", width=12, height=12, dpi=300)

# Extract taxonomy from order-level phyloseq object
tax_order <- as.data.frame(tax_table(ps_order))

# Extract FDR table
da_order_en_FDR <- as.data.frame(da_order_en$p_fdr)

# Add taxon ID and Order name
da_order_en_FDR$taxon_id <- rownames(da_order_en_FDR)
da_order_en_FDR$Order <- tax_order[rownames(da_order_en_FDR), "Order"]

# Keep only significant taxa
sig_ids <- da_order_en$significant_taxa

da_order_en_FDR_sig <- da_order_en_FDR[
  da_order_en_FDR$taxon_id %in% sig_ids,
]

# Reorder columns
da_order_en_FDR_sig <- da_order_en_FDR_sig[
  , c("taxon_id", "Order",
      setdiff(names(da_order_en_FDR_sig), c("taxon_id", "Order")))
]

#extract coefficient mu from the model 

coef_sig <- do.call(rbind, lapply(seq_along(sig_ids), function(i) {
  
  taxon_id_i <- sig_ids[i]
  mod <- da_order_en$significant_models[[i]]
  
  coefs <- coef(mod)
  
  data.frame(
    taxon_id = taxon_id_i,
    term = rownames(coefs),
    coefficient = coefs[, 1],
    row.names = NULL
  )
}))


coef_sig_mu <- subset(coef_sig, grepl("^mu", term))

tax_order <- as.data.frame(tax_table(ps_order))

fdr_table <- as.data.frame(da_order_en$p_fdr)
fdr_table$taxon_id <- rownames(fdr_table)
fdr_table$Order <- tax_order[rownames(fdr_table), "Order"]

fdr_sig <- fdr_table[fdr_table$taxon_id %in% sig_ids, ]


da_order_en_sig_coef_mu <- merge(
  coef_sig_mu,
  fdr_sig,
  by = "taxon_id",
  all.x = TRUE)


write.table(
  da_order_en_sig_coef_mu,
  "da_order_en_significant_with_mu_FDR.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

##adding an interpretation column to the above table 

da_order_en_sig_coef_mu$comparison <- da_order_en_sig_coef_mu$term
da_order_en_sig_coef_mu$comparison <- gsub("mu\\.", "", da_order_en_sig_coef_mu$comparison)
da_order_en_sig_coef_mu$comparison <- gsub("environment", "", da_order_en_sig_coef_mu$comparison)

da_order_en_sig_coef_mu_env <- da_order_en_sig_coef_mu[
  da_order_en_sig_coef_mu$comparison != "(Intercept)",
]

da_order_en_sig_coef_mu_env$enriched_in <- ifelse(
  da_order_en_sig_coef_mu_env$comparison == "seawater" &
    da_order_en_sig_coef_mu_env$coefficient > 0, "seawater",
  ifelse(
    da_order_en_sig_coef_mu_env$comparison == "seawater" &
      da_order_en_sig_coef_mu_env$coefficient < 0, "fish skin",
    ifelse(
      da_order_en_sig_coef_mu_env$comparison == "substrate" &
        da_order_en_sig_coef_mu_env$coefficient > 0, "substrate",
      ifelse(
        da_order_en_sig_coef_mu_env$comparison == "substrate" &
          da_order_en_sig_coef_mu_env$coefficient < 0, "fish skin",
        NA
      )
    )
  )
)

write.table(
  da_order_en_sig_coef_mu_env,
  "da_order_en_significant_mu_coefficients_enriched_in.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

#adding CI interval to the data_table 
p_order <- plot(da_order_en, level = "Order")

plot_data <- p_order$data
head(plot_data)
names(plot_data)

library(dplyr)

ci_table_clean <- plot_data %>%
  mutate(
    Order = taxa,
    comparison = gsub("\nDifferential Abundance", "", variable),
    comparison = gsub("environment", "", comparison),
    coefficient_plot = x,
    CI_excludes_zero = (xmin <= 0 & xmax <= 0) | (xmin >= 0 & xmax >= 0),
    CI_crosses_zero = !CI_excludes_zero
  ) %>%
  select(Order, comparison, coefficient_plot, xmin, xmax,
         CI_excludes_zero, CI_crosses_zero)

final_table <- da_order_en_sig_coef_mu_env %>%
  left_join(
    ci_table_clean,
    by = c("Order", "comparison")
  )
  
write.table(
  final_table,
  "TableS3_DA_orders_environment_FDR_mu_CI.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)  

#counting order enriched in fish skin 


library(dplyr)


fishskin_vs_sw <- final_table %>%
  filter(comparison == "seawater",
         enriched_in == "fish skin",
		 CI_excludes_zero == TRUE) %>%
  distinct(Order)

nrow(fishskin_vs_sw)
#22

cat(fishskin_vs_sw$Order, sep = "\n")
Fusobacteriales
Gammaproteobacteria UO
Spirochaetales
Verrucomicrobiales
Francisellales
Micrococcales
Alphaproteobacteria UO
Enterobacterales
Arenicellales
Kiloniellales
Deferribacterales
Propionibacteriales
Phormidesmiales
Mycoplasmatales
Kordiimonadales
Balneolales
Burkholderiales
Granulosicoccales
Cardiobacteriales
Cyanobacteriia UO
Staphylococcales
Mycobacteriales

#enriched in skin vs substrate 
fishskin_vs_sub <- final_table %>%
  filter(comparison == "substrate",
         enriched_in == "fish skin",
		 CI_excludes_zero == TRUE) %>%
  distinct(Order)

nrow(fishskin_vs_sub)
#27
cat(fishskin_vs_sub$Order, sep = "\n")

Fusobacteriales
Lysobacterales
Bifidobacteriales
Spirochaetales
Verrucomicrobiales
Francisellales
Peptostreptococcales-Tissierellales
Actinomycetales
Micrococcales
Alphaproteobacteria UO
Enterobacterales
Arenicellales
Lactobacillales
Propionibacteriales
Pseudomonadales
Bacillales
Desulfovibrionales
Mycoplasmatales
Kordiimonadales
Rhodospirillales
Rickettsiales
Burkholderiales
Granulosicoccales
Cardiobacteriales
Erysipelotrichales
Staphylococcales
Mycobacteriales


fishskin_vs_both <- final_table %>%
  filter(
    enriched_in == "fish skin",
    CI_excludes_zero == TRUE
  ) %>%
  group_by(Order) %>%
  summarise(n = n_distinct(comparison)) %>%
  filter(n == 2)   # both substrate AND seawater


fishskin_vs_both$Order

cat(fishskin_vs_both$Order, sep = "\n")


1) Alphaproteobacteria UO
2) Arenicellales
3) Burkholderiales
4) Cardiobacteriales
5) Enterobacterales
6) Francisellales
7) Fusobacteriales
8) Granulosicoccales
9) Kordiimonadales
10) Micrococcales
11) Mycobacteriales
12) Mycoplasmatales
13) Propionibacteriales
14) Spirochaetales
15) Staphylococcales
16) Verrucomicrobiales




### ---------------------------------------------------------------------
### Testing for DA among fish species - skin only 
###----------------------------------------------------------------------

#selecting fishskin only and testing difference across fish host species 
ps_skin <- subset_samples(ps,environment == "fishskin") 

ntaxa(ps_skin)

ps_skin_rm0 <-filter_taxa(ps_skin, function(x) sum(x) >0, TRUE)

ntaxa(ps_skin_rm0) 

ps_order_skin <- ps_skin_rm0 %>% tax_glom("Order")


sample_data(ps_order_skin)$species <- factor(
  sample_data(ps_order_skin)$species,
  levels = c("Alticus arnoldorum", 
             "Blenniella paula", 
             "Praealticus labrovittatus")
)


da_order_skin <- differentialTest(formula = ~ species,
											phi.formula = ~ species,
											formula_null = ~ 1,
											phi.formula_null = ~ species,
											test = "Wald", boot = FALSE,
											data = ps_order_skin,
											fdr_cutoff = 0.05)

plot(da_order_skin,level="Order")+labs(x = "Coefficient and Confidence Interval", y = "Order")+ theme(panel.grid = element_blank ())

ggsave("FigureS4-DA_orders-among-fish-species.tiff", units="in", width=12, height=12, dpi=300, compression = 'lzw')
ggsave("FigureS4-DA_orders-among-fish-species.pdf", units="in", width=12, height=12, dpi=300)


# -----------------------------
# Extract significant taxa
# -----------------------------
sig_ids <- da_order_skin$significant_taxa

# -----------------------------
# Extract coefficients
# -----------------------------
coef_sig <- do.call(rbind, lapply(seq_along(sig_ids), function(i) {
  
  taxon_id_i <- sig_ids[i]
  mod <- da_order_skin$significant_models[[i]]
  
  coefs <- coef(mod)
  
  data.frame(
    taxon_id = taxon_id_i,
    term = rownames(coefs),
    coefficient = coefs[, 1],
    row.names = NULL
  )
}))

# keep mu only
coef_sig_mu <- coef_sig %>%
  filter(grepl("^mu", term))

# -----------------------------
# Add Order + FDR
# -----------------------------
tax_order <- as.data.frame(tax_table(ps_order_skin))

fdr_table <- as.data.frame(da_order_skin$p_fdr)
fdr_table$taxon_id <- rownames(fdr_table)
fdr_table$Order <- tax_order[rownames(fdr_table), "Order"]

fdr_sig <- fdr_table %>%
  filter(taxon_id %in% sig_ids)

table_s4 <- coef_sig_mu %>%
  left_join(fdr_sig, by = "taxon_id")

# -----------------------------
# Clean comparison names
# -----------------------------
table_s4 <- table_s4 %>%
  mutate(
    comparison = term,
    comparison = gsub("mu\\.", "", comparison),
    comparison = gsub("species", "", comparison)
  ) %>%
  filter(comparison != "(Intercept)")

# -----------------------------
# Add enriched_in column
# -----------------------------
# reference = first species in factor

ref_species <- levels(sample_data(ps_order_skin)$species)[1]

table_s4 <- table_s4 %>%
  mutate(
    enriched_in = case_when(
      coefficient > 0 ~ comparison,
      coefficient < 0 ~ ref_species,
      TRUE ~ NA_character_
    )
  )

# -----------------------------
# Extract CI from plot
# -----------------------------
p_skin <- plot(da_order_skin, level = "Order")
plot_data <- p_skin$data

ci_clean <- plot_data %>%
  mutate(
    Order = taxa,
    comparison = gsub("\nDifferential Abundance", "", variable),
    comparison = gsub("species", "", comparison),
    coefficient_plot = x,
    CI_excludes_zero = (xmin <= 0 & xmax <= 0) | (xmin >= 0 & xmax >= 0),
    CI_crosses_zero = !CI_excludes_zero
  ) %>%
  select(Order, comparison, xmin, xmax,
         CI_excludes_zero, CI_crosses_zero)

# -----------------------------
# Merge everything
# -----------------------------
TableS4 <- table_s4 %>%
  left_join(ci_clean, by = c("Order", "comparison")) %>%
  select(
    Order, comparison, coefficient, xmin, xmax,
    CI_excludes_zero, CI_crosses_zero,
    enriched_in, taxon_id, term, everything()
  ) %>%
  arrange(Order, comparison)

# -----------------------------
# Save
# -----------------------------
write.table(
  TableS4,
  "TableS4_DA_orders_skin_species_FDR_mu_CI.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


#---------------------------------------------------
#Differenitally abundant genera between host species 
#---------------------------------------------------


ps_genus_skin <- ps_skin_rm0 %>% tax_glom("Genus")

library(dplyr)

# Important: set species levels BEFORE running differentialTest()
sample_data(ps_genus_skin)$species <- factor(
  sample_data(ps_genus_skin)$species,
  levels = c("Alticus arnoldorum", 
             "Blenniella paula", 
             "Praealticus labrovittatus")
)

# If not already done after setting levels, rerun:
da_genus_skin <- differentialTest(
  formula = ~ species,
  phi.formula = ~ species,
  formula_null = ~ 1,
  phi.formula_null = ~ species,
  data = ps_genus_skin,
  test = "Wald",
  boot = FALSE,
  fdr_cutoff = 0.05
)



plot(da_genus_skin, level="Genus") + labs(x = "Coefficient and Confidence Interval", y = "Genera")+ theme(panel.grid = element_blank())
ggsave("FigureS5-DA-genera-among-fish-species.tiff", units="in", width=12, height=12, dpi=300, compression = 'lzw')
ggsave("FigureS5-DA-genea-among-fish-species.tiff", units="in", width=12, height=12, dpi=300, compression = 'lzw')





# -----------------------------
# Table S5: DA genera among host species
# -----------------------------

sig_ids <- da_genus_skin$significant_taxa

coef_sig <- do.call(rbind, lapply(seq_along(sig_ids), function(i) {
  
  taxon_id_i <- sig_ids[i]
  mod <- da_genus_skin$significant_models[[i]]
  coefs <- coef(mod)
  
  data.frame(
    taxon_id = taxon_id_i,
    term = rownames(coefs),
    coefficient = coefs[, 1],
    row.names = NULL
  )
}))

coef_sig_mu <- coef_sig %>%
  filter(grepl("^mu", term))

tax_genus <- as.data.frame(tax_table(ps_genus_skin))

fdr_table <- as.data.frame(da_genus_skin$p_fdr)
fdr_table$taxon_id <- rownames(fdr_table)
fdr_table$Genus <- tax_genus[rownames(fdr_table), "Genus"]

fdr_sig <- fdr_table %>%
  filter(taxon_id %in% sig_ids)

table_s5 <- coef_sig_mu %>%
  left_join(fdr_sig, by = "taxon_id")

table_s5 <- table_s5 %>%
  mutate(
    comparison = term,
    comparison = gsub("mu\\.", "", comparison),
    comparison = gsub("species", "", comparison)
  ) %>%
  filter(comparison != "(Intercept)")

ref_species <- levels(sample_data(ps_genus_skin)$species)[1]

table_s5 <- table_s5 %>%
  mutate(
    enriched_in = case_when(
      coefficient > 0 ~ comparison,
      coefficient < 0 ~ ref_species,
      TRUE ~ NA_character_
    )
  )

# CI from plot object
p_genus <- plot(da_genus_skin, level = "Genus")
plot_data_genus <- p_genus$data

ci_clean <- plot_data_genus %>%
  mutate(
    Genus = taxa,
    comparison = gsub("\nDifferential Abundance", "", variable),
    comparison = gsub("species", "", comparison),
    coefficient_plot = x,
    CI_excludes_zero = (xmin <= 0 & xmax <= 0) | (xmin >= 0 & xmax >= 0),
    CI_crosses_zero = !CI_excludes_zero
  ) %>%
  select(Genus, comparison, coefficient_plot, xmin, xmax,
         CI_excludes_zero, CI_crosses_zero)

TableS5 <- table_s5 %>%
  left_join(ci_clean, by = c("Genus", "comparison")) %>%
  select(
    Genus, comparison, coefficient, xmin, xmax,
    CI_excludes_zero, CI_crosses_zero,
    enriched_in, taxon_id, term, everything()
  ) %>%
  arrange(Genus, comparison)

write.table(
  TableS5,
  "TableS5_DA_genera_skin_species_FDR_mu_CI.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

TableS5
