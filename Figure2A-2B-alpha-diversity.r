library(phyloseq)
library(RColorBrewer)
library(randomcoloR)
library(microViz)
library(microbiome)
library(cowplot)
library(grid)
library(scales)
library(mia)
library(ggpubr)
library(knitr)
library(dplyr)
library(paletteer)

otu <- read.table("final_ASV_table.txt",sep="\t",header=TRUE, row.names=1)
taxon <- read.table("final_taxa_table.txt",sep="\t",header=TRUE,row.names=1)
samples <-read.table("metadata.txt",sep="\t",header=TRUE,row.names=1)


taxon<-as.matrix(taxon)
TAX = tax_table(taxon)


ps <- phyloseq(otu_table(otu, taxa_are_rows=FALSE), 
               sample_data(samples), 
               tax_table(TAX))
ps

#The function below creates a table with selected (or all) diversity indicators.

alpha_div <- microbiome::alpha(ps, index = "all")
alpha_div_table <- kable(alpha_div)
write.table(alpha_div_table, file = "alpha_diversity_results.txt")
# Get the metadata from the `phyloseq` object

ps.meta <- meta(ps)
#Add the diversity table to metadata

ps.meta$Shannon <- alpha_div$diversity_shannon 
ps.meta$InverseSimpson <-  alpha_div$diversity_inverse_simpson
ps.meta$Richness <-  alpha_div$observed


#Plotting shannon diveristy boxplot 

ps.meta$cat <- factor(
  ps.meta$cat,
  levels = c("Bp-S","Bp-R","Pl-S","Pl-R","Aa-S","Aa-R","SW"))
  
  
pShannon = ggplot(ps.meta,aes(x = cat, y = Shannon, color = cat)) +
geom_boxplot(outlier.shape = NA)+
labs(x = "", y = "Shannon Diversity") +
theme_bw()+
theme(panel.grid = element_blank())+
theme(axis.text.x=element_text(size=10))+
theme(legend.text = element_text(size=10))+
theme(legend.title=element_text(size=12))+
scale_color_paletteer_d(
  "MetBrewer::Austria",
  name = "Sample-type",
  breaks = c("Bp-S","Bp-R","Pl-S","Pl-R","Aa-S","Aa-R","SW"),
  labels = c(
    expression("Bp-S:"~italic("Blenniella paula")~" skin"),
    expression("Bp-R:"~italic("Blenniella paula")~" substrate"),
    expression("Pl-S:"~italic("Praealticus labrovittatus")~" skin"),
    expression("Pl-R:"~italic("Praealticus labrovittatus")~" substrate"),
    expression("Aa-S:"~italic("Alticus arnoldorum")~" skin"),
    expression("Aa-R:"~italic("Alticus arnoldorum")~" substrate"),
    "SW: seawater"
  )
)

pShannon
ggsave("Figure2A-shannon_diversity_boxplot.tiff", units="in", width=8, height=4, dpi=300, compression = 'lzw')
ggsave("Figure2A-shannon_diversity_boxplot.pdf", units="in", width=8, height=4, dpi=300)



##plotting ASV Richness 

pRichness = ggplot(ps.meta,aes(x = cat, y = Richness, color = cat)) +
geom_boxplot(outlier.shape = NA)+
labs(x = "", y = "ASV Richness") +
theme_bw()+
theme(panel.grid = element_blank())+
theme(axis.text.x=element_text(size=10))+
theme(legend.text = element_text(size=10))+
theme(legend.title=element_text(size=12))+
scale_color_paletteer_d(
  "MetBrewer::Austria",
  name = "Sample-type",
  breaks = c("Bp-S","Bp-R","Pl-S","Pl-R","Aa-S","Aa-R","SW"),
  labels = c(
    expression("Bp-S:"~italic("Blenniella paula")~" skin"),
    expression("Bp-R:"~italic("Blenniella paula")~" substrate"),
    expression("Pl-S:"~italic("Praealticus labrovittatus")~" skin"),
    expression("Pl-R:"~italic("Praealticus labrovittatus")~" substrate"),
    expression("Aa-S:"~italic("Alticus arnoldorum")~" skin"),
    expression("Aa-R:"~italic("Alticus arnoldorum")~" substrate"),
    "SW: seawater"
  )
)

pRichness
ggsave("Figure2B-ASV-Richiness_boxplot.tiff", units="in", width=8, height=4, dpi=300, compression = 'lzw')
ggsave("Figure2B-ASV-Richiness_boxplot.pdf", units="in", width=8, height=4, dpi=300)



#plotting InverseSimpson 

pInverseSimpson = ggplot(ps.meta,aes(x = cat, y = InverseSimpson, color = cat)) +
geom_boxplot(outlier.shape = NA)+
labs(x = "", y = "Inverse Simpson Diversity") +
theme_bw()+
theme(panel.grid = element_blank())+
theme(axis.text.x=element_text(size=10))+
theme(legend.text = element_text(size=10))+
theme(legend.title=element_text(size=12))+
scale_color_paletteer_d(
  "MetBrewer::Austria",
  name = "Sample-type",
  breaks = c("Bp-S","Bp-R","Pl-S","Pl-R","Aa-S","Aa-R","SW"),
  labels = c(
    expression("Bp-S:"~italic("Blenniella paula")~" skin"),
    expression("Bp-R:"~italic("Blenniella paula")~" substrate"),
    expression("Pl-S:"~italic("Praealticus labrovittatus")~" skin"),
    expression("Pl-R:"~italic("Praealticus labrovittatus")~" substrate"),
    expression("Aa-S:"~italic("Alticus arnoldorum")~" skin"),
    expression("Aa-R:"~italic("Alticus arnoldorum")~" substrate"),
    "SW: seawater"
  )
)

pInverseSimpson
ggsave("Figure2C-InverseSimpson_diversity_boxplot.tiff", units="in", width=8, height=4, dpi=300, compression = 'lzw')
ggsave("Figure2C-InverseSimpson_diversity_boxplot.pdf", units="in", width=8, height=4, dpi=300)




library(dplyr)

# Shannon
ps.meta %>%
  group_by(cat) %>%
  summarise(p_value = shapiro.test(Shannon)$p.value)

# Inverse Simpson
ps.meta %>%
  group_by(cat) %>%
  summarise(p_value = shapiro.test(InverseSimpson)$p.value)

# Richness
ps.meta %>%
  group_by(cat) %>%
  summarise(p_value = shapiro.test(Richness)$p.value)
  
  
library(car)

leveneTest(Shannon ~ cat, data = ps.meta)
leveneTest(InverseSimpson ~ cat, data = ps.meta)
leveneTest(Richness ~ cat, data = ps.meta)



🔬 Summary of results
✅ Shannon
Normality: all groups p > 0.05 → OK
Homogeneity: p = 0.19 → OK

👉 Meets assumptions

⚠️ Inverse Simpson
Normality: violated in:
Aa-S (p = 0.025)
Aa-R (p = 0.043)
Homogeneity: strongly violated (p = 2.3e-06)

👉 Fails assumptions

⚠️ Richness
Normality: violated in Pl-S (p = 0.044)
Homogeneity: OK (p = 0.91)

👉 borderline / mixed


#Shannon paramteric test 

anova(lm(Shannon ~ cat, data = ps.meta))


#Richness non paramteric 

kruskal.test(Richness ~ cat, data = ps.meta)

#InverseSimpson - non paramteric 

kruskal.test(InverseSimpson ~ cat, data = ps.meta)

#Shannon non parametric test 

kruskal.test(Shannon ~ cat, data = ps.meta)


library(FSA)
dunnTest(InverseSimpson ~ cat, data = ps.meta, method = "bh")