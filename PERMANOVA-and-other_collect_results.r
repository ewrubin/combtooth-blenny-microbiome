
#------------------------------------------------------------
# Environment: PERMANOVA, pairwise PERMANOVA, betadisper, Tukey
# fishskin vs substrate vs seawater
# Export all results into one table
#------------------------------------------------------------

library(vegan)
library(pairwiseAdonis)
library(dplyr)
library(tibble)

# Metadata and Bray-Curtis distance matrix
metadata_env <- data.frame(sample_data(ps))

bray_dists <- ps %>% dist_calc("bray")
dist_mat_env <- bray_dists %>% dist_get()

#------------------------------------------------------------
# 1. Overall PERMANOVA
#------------------------------------------------------------

bray_perm_env <- bray_dists %>%
  dist_permanova(
    seed = 1234,
    n_processes = 1,
    n_perms = 999,
    variables = "environment"
  )

permanova_env <- perm_get(bray_perm_env) %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  mutate(
    test = "PERMANOVA",
    comparison = "Overall: environment",
    p_adjust_method = NA
  ) %>%
  rename(
    df = Df,
    sum_of_squares = SumOfSqs,
    r2 = R2,
    f_value = F,
    p_value = `Pr(>F)`
  ) %>%
  select(test, comparison, term, df, sum_of_squares, r2, f_value, p_value, p_adjust_method)

#------------------------------------------------------------
# 2. Pairwise PERMANOVA
#------------------------------------------------------------

pairwise_env <- pairwise.adonis2(
  dist_mat_env ~ environment,
  data = metadata_env,
  permutations = 999,
  p.adjust.m = "BH"
)

pairwise_env_results <- lapply(names(pairwise_env)[-1], function(x) {
  
  res <- pairwise_env[[x]] %>%
    as.data.frame() %>%
    rownames_to_column("term")
  
  res %>%
    filter(term == "Model") %>%
    mutate(
      test = "Pairwise PERMANOVA",
      comparison = gsub("_vs_", " vs ", x),
      p_adjust_method = "BH"
    ) %>%
    rename(
      df = Df,
      sum_of_squares = SumOfSqs,
      r2 = R2,
      f_value = F,
      p_value = `Pr(>F)`
    ) %>%
    select(test, comparison, term, df, sum_of_squares, r2, f_value, p_value, p_adjust_method)
  
}) %>%
  bind_rows()

#------------------------------------------------------------
# 3. Betadisper
#------------------------------------------------------------

bd_env <- betadisper(
  dist_mat_env,
  group = metadata_env$environment
)

bd_anova_env <- anova(bd_env)

betadisper_env <- bd_anova_env %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  filter(term == "Groups") %>%
  mutate(
    test = "Betadisper",
    comparison = "Overall: environment dispersion",
    r2 = NA,
    p_adjust_method = NA
  ) %>%
  rename(
    df = Df,
    sum_of_squares = `Sum Sq`,
    f_value = `F value`,
    p_value = `Pr(>F)`
  ) %>%
  select(test, comparison, term, df, sum_of_squares, r2, f_value, p_value, p_adjust_method)

#------------------------------------------------------------
# 4. Tukey test on betadisper
#------------------------------------------------------------

tukey_env <- TukeyHSD(bd_env)

tukey_env_results <- tukey_env$group %>%
  as.data.frame() %>%
  rownames_to_column("comparison") %>%
  mutate(
    test = "TukeyHSD on betadisper",
    term = "group",
    df = NA,
    sum_of_squares = NA,
    r2 = NA,
    f_value = NA,
    p_adjust_method = "TukeyHSD"
  ) %>%
  rename(
    difference = diff,
    lower_ci = lwr,
    upper_ci = upr,
    p_value = `p adj`
  ) %>%
  select(
    test, comparison, term, df, sum_of_squares, r2,
    f_value, difference, lower_ci, upper_ci,
    p_value, p_adjust_method
  )

#------------------------------------------------------------
# 5. Combine all results
#------------------------------------------------------------

all_env_results <- bind_rows(
  permanova_env %>%
    mutate(difference = NA, lower_ci = NA, upper_ci = NA),
  
  pairwise_env_results %>%
    mutate(difference = NA, lower_ci = NA, upper_ci = NA),
  
  betadisper_env %>%
    mutate(difference = NA, lower_ci = NA, upper_ci = NA),
  
  tukey_env_results
)

# Round numeric columns
all_env_results <- all_env_results %>%
  mutate(across(where(is.numeric), ~round(.x, 4)))

# View table
all_env_results

# Export
write.csv(
  all_env_results,
  "environment_PERMANOVA_pairwise_betadisper_results.csv",
  row.names = FALSE
)



#------------------------------------------------------------
# Fish skin: PERMANOVA, pairwise PERMANOVA, betadisper, Tukey
# Export all results into one table
#------------------------------------------------------------

library(vegan)
library(pairwiseAdonis)
library(dplyr)
library(tibble)

# Metadata and distance matrix
metadata_skin <- data.frame(sample_data(ps_skin_rm0))
dist_mat_skin <- bray_dists_skin %>% dist_get()

#------------------------------------------------------------
# 1. Overall PERMANOVA
#------------------------------------------------------------

bray_perm_s <- bray_dists_skin %>%
  dist_permanova(
    seed = 1234,
    n_processes = 1,
    n_perms = 999,
    variables = "species"
  )

permanova_skin <- perm_get(bray_perm_s) %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  mutate(
    test = "PERMANOVA",
    comparison = "Overall: species",
    p_adjust_method = NA
  ) %>%
  rename(
    df = Df,
    sum_of_squares = SumOfSqs,
    r2 = R2,
    f_value = F,
    p_value = `Pr(>F)`
  ) %>%
  select(test, comparison, term, df, sum_of_squares, r2, f_value, p_value, p_adjust_method)

#------------------------------------------------------------
# 2. Pairwise PERMANOVA
#------------------------------------------------------------

pairwise_skin <- pairwise.adonis2(
  dist_mat_skin ~ species,
  data = metadata_skin,
  permutations = 999,
  p.adjust.m = "BH"
)

pairwise_skin_results <- lapply(names(pairwise_skin)[-1], function(x) {
  
  res <- pairwise_skin[[x]] %>%
    as.data.frame() %>%
    rownames_to_column("term")
  
  res %>%
    filter(term == "Model") %>%
    mutate(
      test = "Pairwise PERMANOVA",
      comparison = gsub("_vs_", " vs ", x),
      p_adjust_method = "BH"
    ) %>%
    rename(
      df = Df,
      sum_of_squares = SumOfSqs,
      r2 = R2,
      f_value = F,
      p_value = `Pr(>F)`
    ) %>%
    select(test, comparison, term, df, sum_of_squares, r2, f_value, p_value, p_adjust_method)
  
}) %>%
  bind_rows()

#------------------------------------------------------------
# 3. Betadisper
#------------------------------------------------------------

bd_skin <- betadisper(
  dist_mat_skin,
  group = metadata_skin$species
)

bd_anova <- anova(bd_skin)

betadisper_skin <- bd_anova %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  filter(term == "Groups") %>%
  mutate(
    test = "Betadisper",
    comparison = "Overall: species dispersion",
    r2 = NA,
    p_adjust_method = NA
  ) %>%
  rename(
    df = Df,
    sum_of_squares = `Sum Sq`,
    f_value = `F value`,
    p_value = `Pr(>F)`
  ) %>%
  select(test, comparison, term, df, sum_of_squares, r2, f_value, p_value, p_adjust_method)

#------------------------------------------------------------
# 4. Tukey test on betadisper
#------------------------------------------------------------

tukey_skin <- TukeyHSD(bd_skin)

tukey_skin_results <- tukey_skin$group %>%
  as.data.frame() %>%
  rownames_to_column("comparison") %>%
  mutate(
    test = "TukeyHSD on betadisper",
    term = "group",
    df = NA,
    sum_of_squares = NA,
    r2 = NA,
    f_value = NA,
    p_adjust_method = "TukeyHSD"
  ) %>%
  rename(
    difference = diff,
    lower_ci = lwr,
    upper_ci = upr,
    p_value = `p adj`
  ) %>%
  select(
    test, comparison, term, df, sum_of_squares, r2,
    f_value, difference, lower_ci, upper_ci,
    p_value, p_adjust_method
  )

#------------------------------------------------------------
# 5. Combine all results
#------------------------------------------------------------

all_skin_results <- bind_rows(
  permanova_skin %>%
    mutate(difference = NA, lower_ci = NA, upper_ci = NA),
  
  pairwise_skin_results %>%
    mutate(difference = NA, lower_ci = NA, upper_ci = NA),
  
  betadisper_skin %>%
    mutate(difference = NA, lower_ci = NA, upper_ci = NA),
  
  tukey_skin_results
)

# optional: round numeric columns
all_skin_results <- all_skin_results %>%
  mutate(across(where(is.numeric), ~round(.x, 4)))

# view table
all_skin_results

# export
write.csv(
  all_skin_results,
  "fishskin_PERMANOVA_pairwise_betadisper_results.csv",
  row.names = FALSE
)


#------------------------------------------------------------
# Substrate: PERMANOVA, pairwise PERMANOVA, betadisper, Tukey
# Export all results into one table
#------------------------------------------------------------

library(vegan)
library(pairwiseAdonis)
library(dplyr)
library(tibble)

# Substrate only
ps_sub <- subset_samples(ps, environment == "substrate")

ntaxa(ps_sub)

ps_sub_rm0 <- filter_taxa(ps_sub, function(x) sum(x) > 0, TRUE)

ntaxa(ps_sub_rm0)

# Metadata and distance matrix
metadata_sub <- data.frame(sample_data(ps_sub_rm0))
bray_dists_sub <- ps_sub_rm0 %>% dist_calc("bray")
dist_mat_sub <- bray_dists_sub %>% dist_get()

#------------------------------------------------------------
# 1. Overall PERMANOVA
#------------------------------------------------------------

bray_perm_sub <- bray_dists_sub %>%
  dist_permanova(
    seed = 1234,
    n_processes = 1,
    n_perms = 999,
    variables = "species"
  )

permanova_sub <- perm_get(bray_perm_sub) %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  mutate(
    test = "PERMANOVA",
    comparison = "Overall: species",
    p_adjust_method = NA
  ) %>%
  rename(
    df = Df,
    sum_of_squares = SumOfSqs,
    r2 = R2,
    f_value = F,
    p_value = `Pr(>F)`
  ) %>%
  select(test, comparison, term, df, sum_of_squares, r2, f_value, p_value, p_adjust_method)

#------------------------------------------------------------
# 2. Pairwise PERMANOVA
#------------------------------------------------------------

pairwise_sub <- pairwise.adonis2(
  dist_mat_sub ~ species,
  data = metadata_sub,
  permutations = 999,
  p.adjust.m = "BH"
)

pairwise_sub_results <- lapply(names(pairwise_sub)[-1], function(x) {
  
  res <- pairwise_sub[[x]] %>%
    as.data.frame() %>%
    rownames_to_column("term")
  
  res %>%
    filter(term == "Model") %>%
    mutate(
      test = "Pairwise PERMANOVA",
      comparison = gsub("_vs_", " vs ", x),
      p_adjust_method = "BH"
    ) %>%
    rename(
      df = Df,
      sum_of_squares = SumOfSqs,
      r2 = R2,
      f_value = F,
      p_value = `Pr(>F)`
    ) %>%
    select(test, comparison, term, df, sum_of_squares, r2, f_value, p_value, p_adjust_method)
  
}) %>%
  bind_rows()

#------------------------------------------------------------
# 3. Betadisper
#------------------------------------------------------------

bd_sub <- betadisper(
  dist_mat_sub,
  group = metadata_sub$species
)

bd_anova_sub <- anova(bd_sub)

betadisper_sub <- bd_anova_sub %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  filter(term == "Groups") %>%
  mutate(
    test = "Betadisper",
    comparison = "Overall: species dispersion",
    r2 = NA,
    p_adjust_method = NA
  ) %>%
  rename(
    df = Df,
    sum_of_squares = `Sum Sq`,
    f_value = `F value`,
    p_value = `Pr(>F)`
  ) %>%
  select(test, comparison, term, df, sum_of_squares, r2, f_value, p_value, p_adjust_method)

#------------------------------------------------------------
# 4. Tukey test on betadisper
#------------------------------------------------------------

tukey_sub <- TukeyHSD(bd_sub)

tukey_sub_results <- tukey_sub$group %>%
  as.data.frame() %>%
  rownames_to_column("comparison") %>%
  mutate(
    test = "TukeyHSD on betadisper",
    term = "group",
    df = NA,
    sum_of_squares = NA,
    r2 = NA,
    f_value = NA,
    p_adjust_method = "TukeyHSD"
  ) %>%
  rename(
    difference = diff,
    lower_ci = lwr,
    upper_ci = upr,
    p_value = `p adj`
  ) %>%
  select(
    test, comparison, term, df, sum_of_squares, r2,
    f_value, difference, lower_ci, upper_ci,
    p_value, p_adjust_method
  )

#------------------------------------------------------------
# 5. Combine all results
#------------------------------------------------------------

all_sub_results <- bind_rows(
  permanova_sub %>%
    mutate(difference = NA, lower_ci = NA, upper_ci = NA),
  
  pairwise_sub_results %>%
    mutate(difference = NA, lower_ci = NA, upper_ci = NA),
  
  betadisper_sub %>%
    mutate(difference = NA, lower_ci = NA, upper_ci = NA),
  
  tukey_sub_results
)

# optional: round numeric columns
all_sub_results <- all_sub_results %>%
  mutate(across(where(is.numeric), ~round(.x, 4)))

# view table
all_sub_results

# export
write.csv(
  all_sub_results,
  "substrate_PERMANOVA_pairwise_betadisper_results.csv",
  row.names = FALSE
)
