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

Fig1A=ps %>%
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
 
 
Fig1A <- ps %>% 
  dist_calc(dist = "bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(
    alpha = 0.6,
    size = 2,
    color = "species",
    shape = "environment"
  ) +
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
  labs(
    x = "PCoA1 [16.3%]",
    y = "PCoA2 [9.9%]"
  ) +
  theme_classic(12) 
 
 
 
 
Fig1A

ggsave("Figure1A-PCoA.tiff", units = "in", width = 6, height = 4, dpi = 300, compression = "lzw")
ggsave("Figure1A-PCoA.pdf", units = "in", width = 6, height = 4, dpi = 300)



#------------------------------------------------------------------
# fishskin and substrate, remove seawater for clarity 
#------------------------------------------------------------------



ps_sub_skin <- subset_samples(ps, environment %in% c("fishskin", "substrate"))

ntaxa(ps_sub_skin)

ps_sub_skin_rm0 <- filter_taxa(ps_sub_skin, function(x) sum(x) > 0, TRUE)

ntaxa(ps_sub_skin_rm0) 



Fig1B <- ps_sub_skin_rm0 %>%
  dist_calc(dist = "bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(
    alpha = 0.6,
    size = 2,
    color = "species",
    shape = "environment"
  ) +
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
  facet_wrap(~environment, scales = "free") +
  labs(
    x = "PCoA1 [18.1%]",
    y = "PCoA2 [9.9%]"
  )

  
Fig1B
  
ggsave("Figure1B.tiff", units="in", width=6, height=4, dpi=300, compression = 'lzw') 
ggsave("Figure1B-PCoA-fishskin-substrate.pdf", units = "in", width = 6, height = 4)

####combined Fig1A and Fig 1B
library(cowplot)

Fig1_combined <- plot_grid(
  Fig1A,
  Fig1B,
  ncol = 1,
  align = "v",
  rel_heights = c(1, 1)
)

Fig1_final <- ggdraw(Fig1_combined) +
  draw_label(
    "A.",
    x = 0.02, y = 0.98,
    hjust = 0,
    size = 12,
    fontface = "plain"
  ) +
  draw_label(
    "B.",
    x = 0.02, y = 0.50,
    hjust = 0,
    size = 12,
    fontface = "plain"
  )

Fig1_final

ggsave(
  "Figure-1-final.pdf",
  Fig1_final,
  width = 8,
  height = 8,
  units = "in",
  device = cairo_pdf
)

ggsave(
  "Figure-1-final.tiff",
  Fig1_final,
  width = 8,
  height = 8,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
#---------------------------------------------------------
#the same plot as above but facet wrap by environment 
#---------------------------------------------------------
FigureS1 <- ps_sub_skin_rm0 %>% 
  dist_calc(dist = "bray") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(
    alpha = 0.6,
    size = 2,
    color = "environment",
    shape = "species"
  ) +
  stat_ellipse(aes(color = environment)) +
  scale_color_manual(
    values = mycolors,
    labels = c(
      "fishskin" = "fish skin",
      "substrate" = "substrate"
    )
  ) +
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
      species = as_labeller(
        c(
          "Blenniella paula" =
            "italic('Blenniella paula')",
          "Praealticus labrovittatus" =
            "italic('Praealticus labrovittatus')",
          "Alticus arnoldorum" =
            "italic('Alticus arnoldorum')"
        ),
        label_parsed
      )
    )
  ) +
  labs(
    x = "PCoA1 [18.1%]",
    y = "PCoA2 [9.9%]"
  ) +
  theme_classic(12) +
  theme(
    strip.text = element_text(size = 9)
  )

FigureS1

ggsave(
  "FigureS1.tiff",
  FigureS1,
  units = "in",
  width = 8,
  height = 4,
  dpi = 300,
  compression = "lzw"
)

ggsave(
  "FigureS1.pdf",
  FigureS1,
  units = "in",
  width = 8,
  height = 4
)
  
ggsave("FigureS1.tiff", units="in", width=8, height=4, dpi=300, compression = 'lzw') 
ggsave("FigureS1.pdf", units = "in", width = 8, height = 4)
