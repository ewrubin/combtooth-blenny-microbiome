#ordination plots and permanova 
https://david-barnett.github.io/microViz/

library(phyloseq)
library(RColorBrewer)
library(randomcoloR)
library(microViz)
library(microbiome)
library(cowplot)
library(grid)
library(scales)
library(vegan)
library(pairwiseAdonis)

otu <- read.table("final_ASV_table.txt", sep="\t", header=TRUE, row.names=1)
taxon <- read.table("final_taxa_table.txt", sep="\t", header=TRUE, row.names=1)
samples <- read.table("metadata.txt", sep="\t", header=TRUE, row.names=1)

## set legend order
samples$species <- factor(
  samples$species,
  levels = c(
    "Blenniella paula",
    "Praealticus labrovittatus",
    "Alticus arnoldorum",
    "seawater"
  )
)

taxon <- as.matrix(taxon)
TAX <- tax_table(taxon)

ps <- phyloseq(
  otu_table(otu, taxa_are_rows = FALSE), 
  sample_data(samples), 
  tax_table(TAX)
)

mycolors <- c(
  "darkolivegreen3",
  "darkgoldenrod3",
  "darkmagenta",
  "cornflowerblue"
)

ps %>%
  dist_calc(dist = "bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.6, size = 2, color = "species", shape = "environment") +
  stat_ellipse(aes(color = species)) +
  scale_color_manual(
    values = mycolors,
    breaks = c(
      "Blenniella paula",
      "Praealticus labrovittatus",
      "Alticus arnoldorum",
      "seawater"
    ),
    labels = c(
      expression(italic("Blenniella paula")),
      expression(italic("Praealticus labrovittatus")),
      expression(italic("Alticus arnoldorum")),
      "seawater"
    )
  ) +
  theme_classic(12)

ggsave("Figure1A-PCoA.tiff", units = "in", width = 6, height = 4, dpi = 300, compression = "lzw")
ggsave("Figure1A-PCoA.pdf", units = "in", width = 6, height = 4, dpi = 300)



#------------------------------------------------------------------
# fishskin and substrate, remove seawater for clarity 
#------------------------------------------------------------------



ps_sub_skin <- subset_samples(ps, environment %in% c("fishskin", "substrate"))

ntaxa(ps_sub_skin)

ps_sub_skin_rm0 <- filter_taxa(ps_sub_skin, function(x) sum(x) > 0, TRUE)

ntaxa(ps_sub_skin_rm0) 

ps_sub_skin_rm0 %>%
  dist_calc(dist = "bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.6, size = 2, color = "species", shape = "environment") +
  theme_classic(12) +
  stat_ellipse(aes(color = species)) +
  scale_color_manual(
    values = mycolors,
    breaks = c(
      "Blenniella paula",
      "Praealticus labrovittatus",
      "Alticus arnoldorum"
    ),
    labels = c(
      expression(italic("Blenniella paula")),
      expression(italic("Praealticus labrovittatus")),
      expression(italic("Alticus arnoldorum"))
    )
  ) +
  facet_wrap(~environment, scales = "free")
  
ggsave("Figure1B.tiff", units="in", width=6, height=4, dpi=300, compression = 'lzw') 
ggsave("Figure1B-PCoA-fishskin-substrate.pdf", units = "in", width = 6, height = 4)



#---------------------------------------------------------
#the same plot as above but facet wrap by environment 
#---------------------------------------------------------

ps_sub_skin_rm0 %>% 
  dist_calc(dist = "bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.6, size = 2, color = "environment", shape = "species") +
  theme_classic(12) +
  stat_ellipse(aes(color = environment)) +
  scale_color_manual(values = mycolors) +
  scale_shape_manual(
    breaks = c(
      "Blenniella paula",
      "Praealticus labrovittatus",
      "Alticus arnoldorum"
    ),
    values = c(16, 17, 15),
    labels = c(
      expression(italic("Blenniella paula")),
      expression(italic("Praealticus labrovittatus")),
      expression(italic("Alticus arnoldorum"))
    )
  ) +
  facet_wrap(
    ~species,
    scales = "free",
    labeller = labeller(
      species = as_labeller(c(
        "Blenniella paula" = "italic('Blenniella paula')",
        "Praealticus labrovittatus" = "italic('Praealticus labrovittatus')",
        "Alticus arnoldorum" = "italic('Alticus arnoldorum')"
      ), label_parsed)
    )
  ) +
  theme(
    strip.text = element_text(size = 9)
  )

  
ggsave("FigureS1.tiff", units="in", width=8, height=4, dpi=300, compression = 'lzw') 
ggsave("FigureS1.pdf", units = "in", width = 8, height = 4)

 
##-------------------------------------------------------------------------
#PERMOANOVA - difference among environemnts (fishskin, substrate, seawater) 
##-------------------------------------------------------------------------

# calculate distances
bray_dists <- ps %>% dist_calc("bray")

#PERMANOVA 
bray_perm <- bray_dists %>%
  dist_permanova(
    seed = 1234, # for set.seed to ensure reproducibility of random process
    n_processes = 1, n_perms = 999, # you should use at least 999!
    variables = "environment")

# view the permanova results
perm_get(bray_perm) %>% as.data.frame()


#-------------------------------------------------------
#permanova and betadispersion on substrate samples to see if different by species 
#--------------------------------------------------------------

#select samples on ps object 
ps_sub <- subset_samples(ps,environment == "substrate")

ntaxa(ps_sub)

ps_sub_rm0 <-filter_taxa(ps_sub, function(x) sum(x) >0, TRUE)

ntaxa(ps_sub_rm0) 

#calculate Bray-Curtis distance 

bray_dists_sub <- ps_sub_rm0 %>% dist_calc("bray")

#permanova

bray_perm_sub <- bray_dists_sub %>%
  dist_permanova(
    seed = 1234, # for set.seed to ensure reproducibility of random process
    n_processes = 1, n_perms = 999, # you should use at least 999!
    variables = "species")

#see permanova results 

perm_get(bray_perm_sub) %>% as.data.frame()


library(vegan)

#betadispersion (bd)

dist_mat <- bray_dists_sub %>% dist_get()

bd <- betadisper(
  dist_mat,
  group = sample_data(ps_sub_rm0)$species
)

anova(bd)
permutest(bd, permutations = 999)

# optional: identify which species differ in dispersion
TukeyHSD(bd)

# optional: visualize dispersion
boxplot(bd)

#pairwiseAdonis

metadata_sub <- data.frame(sample_data(ps_sub_rm0))
dist_mat_sub <- bray_dists_sub %>% dist_get()

pairwise_sub <- pairwise.adonis2(
  dist_mat_sub ~ species,
  data = metadata_sub,
  permutations = 999,
  p.adjust.m = "BH"
)

pairwise_sub


#------------------------------------------------------
#permanova and betadispersion - fishskin only to see species separation 
#------------------------------------------------------

# select fishskin samples 


ps_skin <- subset_samples(ps, environment == "fishskin")

ntaxa(ps_skin)

ps_skin_rm0 <- filter_taxa(ps_skin, function(x) sum(x) > 0, TRUE)

ntaxa(ps_skin_rm0)

# Bray-Curtis distance
bray_dists_skin <- ps_skin_rm0 %>% dist_calc("bray")

# PERMANOVA
bray_perm_s <- bray_dists_skin %>%
  dist_permanova(
    seed = 1234,
    n_processes = 1,
    n_perms = 999,
    variables = "species"
  )

perm_get(bray_perm_s) %>% as.data.frame()

# betadisper test
library(vegan)

dist_mat_skin <- bray_dists_skin %>% dist_get()

metadata_skin <- data.frame(sample_data(ps_skin_rm0))

bd_skin <- betadisper(
  dist_mat_skin,
  group = metadata_skin$species
)

anova(bd_skin)
permutest(bd_skin, permutations = 999)

# optional: identify which species differ in dispersion
TukeyHSD(bd_skin)

# optional: visualize dispersion
boxplot(bd_skin)


#pairwiseAdonis

pairwise_skin <- pairwise.adonis2(
  dist_mat_skin ~ species,
  data = metadata_skin,
  permutations = 999,
  p.adjust.m = "BH"
)

pairwise_skin


