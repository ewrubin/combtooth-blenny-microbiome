library(readxl)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(grid)

# -----------------------------
# 1. Read data
# -----------------------------

fig6 <- read_excel("Fig6A-data.xlsx")
genus_class <- read_excel("data2A_v2.xlsx")

# Check names
names(fig6)
names(genus_class)


# -----------------------------
# 2. Class colors matching Figure 5A + class-specific shading
# -----------------------------

# Base colors for the top annotation and Class legend
class_cols <- c(
  "Bacteroidia"         = "cornsilk4",
  "Cyanobacteriia"      = "maroon4",
  "Fusobacteriia"       = "darkolivegreen",
  "Alphaproteobacteria" = "cornflowerblue",
  "Gammaproteobacteria" = "darkslateblue",
  "Verrucomicrobiia"    = "darkslategray"
)

# Create 10 shades for each Class.
# Shade [1] is white, so for the heatmap we use [2] for 0% and [10] for 40%.
cornsilk4_shades <- colorRampPalette(c("white", "cornsilk4"))(10)
maroon4_shades <- colorRampPalette(c("white", "maroon4"))(10)
darkolivegreen_shades <- colorRampPalette(c("white", "darkolivegreen"))(10)
blue_shades <- colorRampPalette(c("white", "cornflowerblue"))(10)
blue_shades2 <- colorRampPalette(c("white", "darkslateblue"))(10)
grayshades <- colorRampPalette(c("white", "darkslategray"))(10)

max_value <- 40

class_color_funs <- list(
  "Bacteroidia" = colorRamp2(
    c(0, max_value),
    c(cornsilk4_shades[2], cornsilk4_shades[10])
  ),
  "Cyanobacteriia" = colorRamp2(
    c(0, max_value),
    c(maroon4_shades[2], maroon4_shades[10])
  ),
  "Fusobacteriia" = colorRamp2(
    c(0, max_value),
    c(darkolivegreen_shades[2], darkolivegreen_shades[10])
  ),
  "Alphaproteobacteria" = colorRamp2(
    c(0, max_value),
    c(blue_shades[2], blue_shades[10])
  ),
  "Gammaproteobacteria" = colorRamp2(
    c(0, max_value),
    c(blue_shades2[2], blue_shades2[10])
  ),
  "Verrucomicrobiia" = colorRamp2(
    c(0, max_value),
    c(grayshades[2], grayshades[10])
  )
)

# Function used inside cell_fun. Values above 40 are capped at the darkest shade.
shade_fun <- function(value, class_name) {
  value <- pmin(value, max_value)
  class_color_funs[[class_name]](value)
}

# Optional: inspect the exact shades
# scales::show_col(coral3_shades)
# scales::show_col(maroon4_shades)
# scales::show_col(darkolivegreen_shades)
# scales::show_col(blue_shades)
# scales::show_col(blue_shades3)
# scales::show_col(grayshades)

# -----------------------------
# 3. Put Fig6 data into long format
# -----------------------------

fig6_long <- fig6 %>%
  pivot_longer(
    cols = -c(short_briteB, short_briteC),
    names_to = "Genus",
    values_to = "percent"
  ) %>%
  mutate(percent = as.numeric(percent))

# -----------------------------
# 4. Add Class information
# -----------------------------

fig6_long <- fig6_long %>%
  left_join(
    genus_class %>% distinct(Genus, Class),
    by = "Genus"
  )

# Check if any genera failed to match
fig6_long %>%
  filter(is.na(Class)) %>%
  distinct(Genus)

# -----------------------------
# 5. Make heatmap matrix
# -----------------------------

mat_df <- fig6_long %>%
  select(short_briteB, short_briteC, Genus, percent) %>%
  pivot_wider(
    names_from = Genus,
    values_from = percent,
    values_fill = 0
  )

row_info <- mat_df %>%
  select(short_briteB, short_briteC)

mat <- mat_df %>%
  select(-short_briteB, -short_briteC) %>%
  as.data.frame()

rownames(mat) <- mat_df$short_briteC
mat <- as.matrix(mat)

# -----------------------------
# 6. Optional: set genus order
# -----------------------------

genus_order <- c("Dokdonia","Lewinella","Flavobacteriaceae UG","Maribacter","Phormidesmis","Schizothrix",
                 "Cetobacterium","Kordiimonas","Ruegeria","Acinetobacter","Arenicella","Burkholderiales UG",
                 "Cardiobacteriaceae UG","Gammaproteobacteria UG","Granulosicoccus","Parendozoicomonas","Woeseia","Rubritalea")

genus_order <- genus_order[genus_order %in% colnames(mat)]
mat <- mat[, genus_order, drop = FALSE]

# Match class info to column order
col_class <- genus_class %>%
  distinct(Genus, Class) %>%
  filter(Genus %in% colnames(mat)) %>%
  arrange(match(Genus, colnames(mat)))

stopifnot(all(col_class$Genus == colnames(mat)))

# -----------------------------
# 7. Function to shade by value within class color
# -----------------------------

# Already defined above as shade_fun(), using shade[2] for 0% and shade[10] for 40%.

# -----------------------------
# 8. Column annotation: Class color bar
# -----------------------------

top_ha <- HeatmapAnnotation(
  Class = col_class$Class,
  col = list(Class = class_cols),
  annotation_name_side = "left",
  show_annotation_name = TRUE,
  show_legend = FALSE)


# -----------------------------
# 9. Row split by KEGG BRITE B category
# -----------------------------

row_split <- row_info$short_briteB
names(row_split) <- row_info$short_briteC
row_split <- row_split[rownames(mat)]

row_split <- dplyr::recode(
  row_split,

  "Energy Metab" = "Energy\nMetabolism",
  "Lipid Metab" = "Lipid\nMetabolism",
  "Terpenoids -Polyketides" = "Terpenoids\n& Polyketides",
  "Secondary Metabolites" = "Secondary\nMetabolites",
  "Glycan Biosyn Metab" = "Glycan\nBiosynthesis\n& Metabolism",
  "Cabrohydrate Metab" = "Carbohydrate\nMetabolism",
  "Cofactors-Vitamins" = "Cofactors\n& Vitamins",
  "Amino acid" = "Amino Acid\nMetabolism",
  "Nucleotide Metab" = "Nucleotide\nMetabolism",
  "Xenobiotics" = "Xenobiotics\nBiodegradation\n& Metabolism",

  .default = row_split
)



# -----------------------------
# 10. Draw heatmap
# -----------------------------

ht <- Heatmap(
  mat,
  name = "Percent\ncontribution",
  
  top_annotation = top_ha,
  row_split = row_split,
  
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  show_row_dend = TRUE,
  show_column_dend = FALSE,
  
  row_names_side = "right",
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8,fontface = "italic"),
  column_names_rot = 45,
  row_title_gp = gpar(fontsize = 6, fontface = "bold"),
  show_heatmap_legend = FALSE,
  
  rect_gp = gpar(col = "grey85", lwd = 0.4),
  
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    value <- mat[i, j]
    class_name <- col_class$Class[j]
    cell_col <- shade_fun(value, class_name)
    
    grid.rect(
      x = x, y = y,
      width = width,
      height = height,
      gp = gpar(fill = cell_col, col = "grey85", lwd = 0.4)
    )
    
    if (!is.na(value) && value >= 3) {
      grid.text(
        round(value),
        x = x, y = y,
        gp = gpar(fontsize = 5, col = "black")
      )
    }
  }
)

# Class legend
class_legend <- Legend(
  title = "Class",
  labels = gt_render(paste0("*", names(class_cols), "*")),
  legend_gp = gpar(fill = class_cols))

# Percent legend, generic blue scale for intensity
percent_legend <- Legend(
  title = "Percent\ncontribution",
  col_fun = colorRamp2(
    c(0, 10, 20, 30, 40),
    c("grey95", "grey82", "grey65", "grey45", "grey25")
  ),
  at = c(0, 10, 20, 30, 40)
)

draw(
  ht,
  annotation_legend_list = list(class_legend, percent_legend),
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)

# -----------------------------
# 11. Save figure
# -----------------------------

pdf("Figure5_v_5_20_2026.pdf", width = 8, height = 12)
draw(
  ht,
  annotation_legend_list = list(class_legend, percent_legend),
  heatmap_legend_side = "right",
  annotation_legend_side = "right"
)
dev.off()

