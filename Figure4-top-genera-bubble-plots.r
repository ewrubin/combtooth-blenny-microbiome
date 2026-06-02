library(ggplot2)
library(reshape2)
library(scales)
library(tidyverse) 



data1A <- read.table("data1A_v2.txt", sep = "\t", header = TRUE)
data2A <- read.table("data2A_v2.txt", sep = "\t", header = TRUE)
data3 <- read.table("data3.txt", sep = "\t", header = TRUE)

data1A$Genus <- stringr::str_squish(data1A$Genus)
data2A$Genus <- stringr::str_squish(data2A$Genus)

#Figure 5A 
data.long = melt(data1A, id = c("Genus"),variable.name = "sample", value.name = "proportion")
data.long.taxon <- merge(data.long, data2A, by.x = "Genus", by.y = "Genus")
data.long.taxon.metadata <- merge(data.long.taxon, data3, by.x = "sample", by.y = "sample")

data.long.taxon.metadata$Class <- factor(data.long.taxon.metadata$Class,levels=c("Bacteroidia",
																				"Cyanobacteriia",
																				"Fusobacteriia",
																				"Alphaproteobacteria",
																				"Gammaproteobacteria",
																				"Verrucomicrobiia"))
																				
																				
data.long.taxon.metadata$Genus <- factor(data.long.taxon.metadata$Genus,levels=c("Dokdonia",
																				"Lewinella",
																				"Flavobacteriaceae UG",
																				"Maribacter",
																				"Phormidesmis",
																				"Schizothrix",
																				"Cetobacterium",
																				"Kordiimonas",
																				"Ruegeria",
																				"Acinetobacter",
																				"Arenicella",
																				"Burkholderiales UG",
																				"Cardiobacteriaceae UG",
																				"Gammaproteobacteria UG",
																				"Granulosicoccus",
																				"Parendozoicomonas",
																				"Woeseia",
																				"Rubritalea"))

colors = c("cornsilk3",
			"maroon3",
			"darkolivegreen",
			"#B8D0F7",
			"#6495DE",
			"darkslategray")


data.long.taxon.metadata$species <- factor(data.long.taxon.metadata$species,levels=c("Blenniella paula","Praealticus labrovittatus","Alticus arnoldorum"))



show_col(colors)

max(data.long.taxon.metadata$proportion)
min(data.long.taxon.metadata$proportion)

Fig5A = ggplot(data.long.taxon.metadata, aes(x = sample, y = Genus, fill = Class)) + 
  geom_point(aes(size = proportion), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0,0.50), range = c(0,15), breaks = c(0,0.05,0.1,0.2,0.3,0.4,0.5)) + 
  labs(size = "Proportion", fill = "Class") +
  scale_y_discrete(limits = rev(levels(data.long.taxon.metadata$Genus)))+
  scale_x_discrete(position = "top")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position="right",axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_fill_manual(values = colors)+
  facet_grid(~species, scales="free")+theme(strip.text = element_text(size = 8, color = "black"))+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(size = 8, face = "italic", color = "black"))+
  theme(axis.text.y = element_text(size = 10, face = "italic"))+
  guides(
  fill = guide_legend(
    title.theme = element_text(face = "italic"),
    label.theme = element_text(face = "italic")
  ),
  size = guide_legend(
    title.theme = element_text(face = "plain"),
    label.theme = element_text(face = "plain")
  )
)

Fig5A

ggsave("Figure-5A-top-genera-SDA-8x6.tiff", units="in", width=8, height=6, dpi=300, compression = 'lzw')
ggsave("Figure-5A-top-genera-SDA-8x6.pdf", units="in", width=8, height=6, dpi=300)



#Figure5B plotting 

data1B <- read.table("data1B_v2.txt", sep = "\t", header = TRUE)
data2B <- read.table("data2B_v2.txt", sep = "\t", header = TRUE)
data3 <- read.table("data3.txt", sep = "\t", header = TRUE)
data.long = melt(data1B, id = c("Genus"),variable.name = "sample", value.name = "proportion")

data.long.taxon <- merge(data.long, data2B, by.x = "Genus", by.y = "Genus")

data.long.taxon.metadata <- merge(data.long.taxon, data3, by.x = "sample", by.y = "sample")

data.long.taxon.metadata$Class <- factor(data.long.taxon.metadata$Class,levels=c("Brevinematia",
																				"Fusobacteriia",
																				"Alphaproteobacteria",
																				"Gammaproteobacteria",
																				"Lentisphaeria"))
																					
data.long.taxon.metadata$Genus <- factor(data.long.taxon.metadata$Genus,levels=c("Brevinema",
																				"Propionigenium",
																				"Alphaproteobacteria UG",
																				"Marivibrio",
																				"Paracoccaceae UG",
																				"Alteromonas",
																				"Cobetia",
																				"Photobacterium",
																				"Vibrio",
																				"Vibrionaceae UG",
																				"Woeseia",
																				"Lentisphaera"))																					
																					
colors2 = c("coral2",
			"darkolivegreen",
			"#B8D0F7",
			"#6495DE",
			"darkgoldenrod2")

data.long.taxon.metadata$species <- factor(data.long.taxon.metadata$species,levels=c("Blenniella paula","Praealticus labrovittatus","Alticus arnoldorum"))


show_col(colors2)

max(data.long.taxon.metadata$proportion)
min(data.long.taxon.metadata$proportion)																					

Fig5B = ggplot(data.long.taxon.metadata, aes(x = sample, y = Genus, fill = Class)) + 
  geom_point(aes(size = proportion), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0,0.50), range = c(0,15), breaks = c(0,0.05,0.1,0.2,0.3,0.4)) + 
  labs(size = "Proportion", fill = "Class") +
  scale_y_discrete(limits = rev(levels(data.long.taxon.metadata$Genus)))+
  scale_x_discrete(position = "top")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),legend.position="right",axis.text.x=element_blank(),axis.ticks.x=element_blank())+
  scale_fill_manual(values = colors2)+
  facet_grid(~species, scales="free")+theme(strip.text = element_text(size = 8, color = "black"))+
  theme(strip.background =element_rect(fill="white"))+
  theme(strip.text = element_text(size = 8, face = "italic", color = "black"))+
  theme(axis.text.y = element_text(size = 10, face = "italic"))+
  guides(
  fill = guide_legend(
    title.theme = element_text(face = "italic"),
    label.theme = element_text(face = "italic")
  ),
  size = guide_legend(
    title.theme = element_text(face = "plain"),
    label.theme = element_text(face = "plain")
  )
)

Fig5B


ggsave("Figure-5B-top-genera-notDA-8x5.tiff", units="in", width=8, height=5, dpi=300, compression = 'lzw')
ggsave("Figure-5B-top-genera-notDA-8x5.pdf", units="in", width=8, height=5, dpi=300)

library(cowplot)

Fig5_combined <- plot_grid(
  Fig5A,
  Fig5B,
  ncol = 1,
  align = "v",
  rel_heights = c(1.25, 1)
)

Fig5_final <- ggdraw(Fig5_combined) +
  draw_label(
    "A. SDA top genera",
    x = 0.02, y = 0.985,
    hjust = 0,
    size = 12,
    fontface = "plain"
  ) +
  draw_label(
    "B. Not SDA top genera",
    x = 0.005, y = 0.44,
    hjust = 0,
    size = 12,
    fontface = "plain"
  )

Fig5_final

ggsave(
  "Figure-5-final.pdf",
  Fig5_final,
  width = 8,
  height = 10,
  units = "in",
  device = cairo_pdf
)

ggsave(
  "Figure-5-final.tiff",
  Fig5_final,
  width = 8,
  height = 10,
  units = "in",
  dpi = 600,
  compression = "lzw"
)


