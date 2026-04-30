library(dplyr)
library(tibble)
library(stringr)
library(dplyr)
library(stringr)
library(tidyverse)

options(timeout = 300)

download.file(
  url      = "https://rest.kegg.jp/get/br:ko00001",
  destfile = "kegg_brite_ko00001.txt",
  quiet    = TRUE
)

brite_raw <- readLines(tmp_brite)
head(brite_raw, 20)

brite_raw <- readLines("kegg_brite_ko00001.txt")


brite_A_lines <- brite_raw[str_starts(brite_raw, "A")]
brite_A_lines

"A09100 Metabolism"                           
"A09120 Genetic Information Processing"       
"A09130 Environmental Information Processing"
"A09140 Cellular Processes"                   
"A09150 Organismal Systems"                   
"A09160 Human Diseases"                      
"A09180 Brite Hierarchies"                    
"A09190 Not Included in Pathway or Brite"    


keep_A <- c(
  "Metabolism",
  "Genetic Information Processing",
  "Environmental Information Processing",
  "Cellular Processes",
  "Brite Hierarchies",                
  "Not Included in Pathway or Brite")

parse_brite_with_A <- function(brite_lines) {
  lvlA <- lvlB <- lvlC <- NA_character_
  out <- list()

  for (ln in brite_lines) {
    if (str_starts(ln, "A")) {
      lvlA <- str_trim(sub("^A\\d*\\s*", "", ln))
    } else if (str_starts(ln, "B")) {
      lvlB <- str_trim(sub("^B\\d*\\s*", "", ln))
    } else if (str_starts(ln, "C")) {
      lvlC <- str_trim(sub("^C\\d*\\s*", "", ln))
    } else if (str_starts(ln, "D")) {
      ko <- str_extract(ln, "K\\d{5}")
      if (!is.na(ko)) {
        out[[length(out) + 1]] <- tibble(
          KO = ko,
          brite_A = lvlA,
          brite_B = lvlB,
          brite_C = lvlC
        )
      }
    }
  }
  bind_rows(out)
}

brite_ko_full <- parse_brite_with_A(brite_raw)


head(brite_ko_full)
nrow(brite_ko_full)


distinct(brite_ko_full$brite_A)


brite_ko_prok <- brite_ko_full %>%
  filter(brite_A %in% keep_A)


table(brite_ko_prok$brite_A)



write.csv(
  brite_ko_prok,
  "brite-database-2.csv",
  row.names = FALSE)

library(readr)
library(dplyr)
library(stringr)


# 1) Load BRITE database as a table (NOT readLines)
brite_db2 <- read_csv("brite-database.csv", show_col_types = FALSE)

# 2) Clean category strings

brite_clean <- brite_db2 %>%
  mutate(
    KO = sub("^ko:", "", KO),

    brite_A = str_trim(brite_A),

    brite_B = brite_B %>%
      str_remove("^\\d+\\s+") %>%          # drop leading digits + space
      str_remove("\\s*\\[BR:.*\\]$") %>%   # drop trailing [BR:...]
      str_trim(),

    brite_C = brite_C %>%
      str_remove("^\\d+\\s+") %>%              # drop leading digits + space
      str_remove("\\s*\\[PATH:.*\\]$") %>%     # drop trailing [PATH:...]
      str_remove("\\s*\\[BR:.*\\]$") %>%       # drop trailing [BR:...]
      str_trim()
  )

write_csv(brite_clean, "brite_database_long_clean.csv")



# 3) Collapse to one row per KO with comma-separated unique categories
brite_collapsed <- brite_clean %>%
  group_by(KO) %>%
  summarise(
    brite_A_all = paste(sort(unique(na.omit(brite_A))), collapse = ", "),
    brite_B_all = paste(sort(unique(na.omit(brite_B))), collapse = ", "),
    brite_C_all = paste(sort(unique(na.omit(brite_C))), collapse = ", "),
    n_brite_B = n_distinct(brite_B, na.rm = TRUE),
    n_brite_C = n_distinct(brite_C, na.rm = TRUE),
    .groups = "drop"
  )

# 4) Quick sanity check
head(brite_collapsed, 10)

# 5) Save for future use
write_csv(brite_collapsed, "brite_database_short_clean.csv")


---------------------------------------------------------------------------
#summary of KOs is each brite_C category 

# Agglomerate KEGG BRITE by brite_C and count how many unique KOs are in each category

# Packages
library(readr)
library(dplyr)
library(stringr)

# ---- Input ----
infile <- "brite_database_long_clean.csv"  # <-- change path if needed

brite <- read_csv(infile, show_col_types = FALSE) %>%
  # basic cleanup (optional but helpful)
  mutate(
    KO = str_trim(KO),
    brite_A = str_trim(brite_A),
    brite_B = str_trim(brite_B),
    brite_C = str_trim(brite_C)
  ) %>%
  filter(!is.na(KO), KO != "", !is.na(brite_C), brite_C != "")

# ---- 1) Count KOs per brite_C (unique KOs) ----
kos_per_briteC <- brite %>%
  group_by(brite_C) %>%
  summarise(
    n_KOs = n_distinct(KO),
    .groups = "drop"
  ) %>%
  arrange(desc(n_KOs), brite_C)

print(kos_per_briteC, n = 50)

# Save results
write_csv(kos_per_briteC, "KO_counts_by_brite_C.csv")

# ---- 2) (Optional) Keep hierarchy context: brite_A / brite_B / brite_C ----
kos_per_briteC_with_context <- brite %>%
  group_by(brite_A, brite_B, brite_C) %>%
  summarise(
    n_KOs = n_distinct(KO),
    .groups = "drop"
  ) %>%
  arrange(desc(n_KOs), brite_A, brite_B, brite_C)

write_csv(kos_per_briteC_with_context, "KO_counts_by_brite_C_with_context.csv")

# ---- 3) (Optional) “Agglomerate” mapping table: each KO -> brite_C list ----
# Useful if a KO appears in multiple brite_C categories.
ko_to_briteC <- brite %>%
  distinct(KO, brite_C) %>%
  group_by(KO) %>%
  summarise(
    brite_C_list = paste(sort(unique(brite_C)), collapse = " | "),
    n_briteC = n_distinct(brite_C),
    .groups = "drop"
  ) %>%
  arrange(desc(n_briteC), KO)

write_csv(ko_to_briteC, "KO_to_brite_C_mapping.csv")

