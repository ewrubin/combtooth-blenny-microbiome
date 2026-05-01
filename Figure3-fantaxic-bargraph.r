library(ggplot2)
library(phyloseq)
library(dplyr)
library(reshape2)
library(tibble)
library(RColorBrewer)
library(cowplot)
library(grid)
library(scales)
library(fantaxtic)



otu <- read.table("final_ASV_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("final_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples <-read.table("metadata.txt",sep="\t",header=TRUE,row.names=1)


taxon<-as.matrix(taxon)
TAX = tax_table(taxon)


ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(TAX))

ps_class <- ps %>% tax_glom("Class")

meta <- data.frame(sample_data(ps))

meta$group <- dplyr::case_when(
  meta$species == "Blenniella paula" & meta$environment == "fishskin"  ~ "Bp-S",
  meta$species == "Blenniella paula" & meta$environment == "substrate" ~ "Bp-R",
  meta$species == "Praealticus labrovittatus" & meta$environment == "fishskin"  ~ "Pl-S",
  meta$species == "Praealticus labrovittatus" & meta$environment == "substrate" ~ "Pl-R",
  meta$species == "Alticus arnoldorum" & meta$environment == "fishskin"  ~ "Aa-S",
  meta$species == "Alticus arnoldorum" & meta$environment == "substrate" ~ "Aa-R",
  meta$environment == "seawater" ~ "SW",
  TRUE ~ NA_character_
)

meta$group <- factor(
  meta$group,
  levels = c("Bp-S", "Bp-R", "Pl-S", "Pl-R", "Aa-S", "Aa-R", "SW")
)

sample_data(ps) <- sample_data(meta)

top <- top_taxa(
  ps_obj = ps,
  n_taxa = 20,
  tax_level = "Class",
  merged_label = "Other",
  FUN = mean
)

facet_labels <- c(
  "Bp-S" = "italic('Blenniella paula')~'skin'",
  "Bp-R" = "italic('Blenniella paula')~'substrate'",
  
  "Pl-S" = "italic('Praealticus labrovittatus')~'skin'",
  "Pl-R" = "italic('Praealticus labrovittatus')~'substrate'",
  
  "Aa-S" = "italic('Alticus arnoldorum')~'skin'",
  "Aa-R" = "italic('Alticus arnoldorum')~'substrate'",
  
  "SW" = "'Seawater'"
)


fantaxtic.bar=plot_nested_bar(ps_obj = top$ps_obj,
                top_level = "Phylum",
                nested_level = "Class",
                palette = c(Actinomycetota = "lightpink",Bacillota = "darkred", 
				Bacteroidota = "cornsilk3",Balneolota = "antiquewhite1",
				Campylobacterota = "cyan",Cyanobacteriota= "darkmagenta",
				Deinococcota = "mediumpurple",Fusobacteriota = "darkolivegreen",
				Myxococcota = "burlywood2",Planctomycetota = "aquamarine3",
				Pseudomonadota = "cornflowerblue",Spirochaetota = "lightskyblue1",
				Thermodesulfobacteriota = "darkgoldenrod1",Thermoplasmatota = "darkgoldenrod2",
				Thermoproteota = "lightgoldenrod1",Unclassified_Bacteria = "lightgray",
				Verrucomicrobiota = "darkslategray",Other = "azure4"))

fantaxtic.bar +
  facet_wrap(~ group, scales = "free", ncol = 3,labeller = as_labeller(facet_labels, label_parsed)) +
  theme(strip.background = element_rect(color = "black", fill = "white",
                                        linewidth = 0.5, linetype = "solid")) +
  theme(strip.text = element_text(color = "black")) +
  theme(axis.text.x = element_blank()) +
  labs(x = "", y = "Relative Abundance") +
  theme(legend.position = "bottom",
        legend.title = element_blank())
	
ggsave("Figure3-fantaxtic-bargraph-top20-classes.tiff", units="in", width=12, height=10, dpi=300, compression = 'lzw')

ggsave("Figure3-fantaxtic-bargraph-top20-classes.pdf", units="in", width=12, height=10, dpi=300)





fig3 <- fantaxtic.bar +
  facet_wrap(
    ~ group,
    scales = "free",
    ncol = 3,
    labeller = as_labeller(facet_labels, label_parsed)
  ) +
  theme(strip.background = element_rect(color = "black", fill = "white",
                                        linewidth = 0.5, linetype = "solid")) +
  theme(strip.text = element_text(color = "black")) +
  theme(axis.text.x = element_blank()) +
  labs(x = "", y = "Relative Abundance") +
  theme(legend.position = "bottom",
        legend.title = element_blank())

fig3



ggsave("Figure3-fantaxtic-bargraph-top20-classes-2.tiff",
       plot = fig3,
       units = "in", width = 12, height = 10,
       dpi = 300, compression = "lzw")

ggsave("Figure3-fantaxtic-bargraph-top20-classes-2.pdf",
       plot = fig3,
       units = "in", width = 12, height = 10,
       dpi = 300)