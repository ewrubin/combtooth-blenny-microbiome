# PERMANOVA and beta-dispersion analyses reported in the manuscript
# Bray-Curtis dissimilarity; 999 permutations; seed = 1234 for PERMANOVA

library(phyloseq)
library(microViz)
library(vegan)

# -------------------------------------------------------------------------
# 1. Overall comparison among environments
#    Model: Bray-Curtis dissimilarity ~ environment
#    Groups: fish skin, substrate, seawater
# -------------------------------------------------------------------------

bray_dists_env <- ps %>%
  dist_calc("bray")

# PERMANOVA
bray_perm_env <- bray_dists_env %>%
  dist_permanova(
    seed = 1234,
    n_processes = 1,
    n_perms = 999,
    variables = "environment"
  )

perm_get(bray_perm_env) %>%
  as.data.frame()

# Homogeneity of multivariate dispersion
dist_mat_env <- bray_dists_env %>%
  dist_get()

bd_env <- betadisper(
  dist_mat_env,
  group = sample_data(ps)$environment
)

permutest(
  bd_env,
  permutations = 999
)


# -------------------------------------------------------------------------
# 2. Fish-skin communities among host species
#    Model: Bray-Curtis dissimilarity ~ species
# -------------------------------------------------------------------------

ps_skin <- subset_samples(
  ps,
  environment == "fishskin"
)

ps_skin_rm0 <- prune_taxa(
  taxa_sums(ps_skin) > 0,
  ps_skin
)

bray_dists_skin <- ps_skin_rm0 %>%
  dist_calc("bray")

# PERMANOVA
bray_perm_skin <- bray_dists_skin %>%
  dist_permanova(
    seed = 1234,
    n_processes = 1,
    n_perms = 999,
    variables = "species"
  )

perm_get(bray_perm_skin) %>%
  as.data.frame()

# Homogeneity of multivariate dispersion
dist_mat_skin <- bray_dists_skin %>%
  dist_get()

bd_skin <- betadisper(
  dist_mat_skin,
  group = sample_data(ps_skin_rm0)$species
)

permutest(
  bd_skin,
  permutations = 999
)


# -------------------------------------------------------------------------
# 3. Substrate communities among host-associated habitats
#    Model: Bray-Curtis dissimilarity ~ species
# -------------------------------------------------------------------------

ps_sub <- subset_samples(
  ps,
  environment == "substrate"
)

ps_sub_rm0 <- prune_taxa(
  taxa_sums(ps_sub) > 0,
  ps_sub
)

bray_dists_sub <- ps_sub_rm0 %>%
  dist_calc("bray")

# PERMANOVA
bray_perm_sub <- bray_dists_sub %>%
  dist_permanova(
    seed = 1234,
    n_processes = 1,
    n_perms = 999,
    variables = "species"
  )

perm_get(bray_perm_sub) %>%
  as.data.frame()

# Homogeneity of multivariate dispersion
dist_mat_sub <- bray_dists_sub %>%
  dist_get()

bd_sub <- betadisper(
  dist_mat_sub,
  group = sample_data(ps_sub_rm0)$species
)

permutest(
  bd_sub,
  permutations = 999
)
