setwd("C:/Users/karlp/OneDrive/Desktop/OTUsequencing2025")

# Read your CSV file
df <- read.csv("OTU_table_ECM_and_nonECM.csv")

# Replace all non-zero values with 1
df[df != 0] <- 1

# Save the modified table
write.csv(df, "OTU_table_ECM_and_nonECM_converted.csv", row.names = FALSE)


## ---- Setup ----
# Install (once) if needed:
# install.packages(c("lme4", "MuMIn", "DHARMa", "performance"))

library(lme4)        # glmer.nb
library(MuMIn)       # model selection (dredge, AICc, r.squaredGLMM)
library(DHARMa)      # diagnostics
library(performance) # extra checks (optional)
library(car)

## ---- Load data ----
dat <- read.csv("ECMASV_GLMM_OTU.csv")

# Ensure expected columns exist
stopifnot(all(c("alpha","pH","clay","site") %in% names(dat)))

# Coerce types
dat$site <- as.factor(dat$site)

# Quick sanity checks
summary(dat[c("alpha","pH","clay","site")])

## ---- Full model (NB-GLMM) ----
m_full <- glmer.nb(alpha ~ pH + clay + sand + (1 | site), data = dat)
cat("\n--- Full model summary ---\n")
print(summary(m_full))

##colinearity
model <- glmer.nb(alpha ~ pH + silt + clay + sand + (1 | site), data = dat)

# Compute VIFs
vif(model)

## ---- Model selection (AICc, all subsets of fixed effects) ----
# Required for dredge
op_old <- options(na.action = "na.fail")
on.exit(options(op_old), add = TRUE)

ms <- dredge(m_full)                  # ranks by AICc
cat("\n--- Model selection table (AICc) ---\n")
print(ms)

best_model <- get.models(ms, 1)[[1]]  # best AICc model
cat("\n--- Best model formula ---\n")
print(formula(best_model))

## ---- Best model summary (coefficients + p-values) ----
cat("\n--- Best model summary ---\n")
sbest <- summary(best_model)
print(sbest)

# Neat coefficient table (Estimate, SE, z, p)
cat("\n--- Coefficients (fixed effects) ---\n")
coef_tab <- as.data.frame(coef(summary(best_model)))
names(coef_tab) <- c("Estimate","Std.Error","z.value","Pr(>|z|)")
print(coef_tab)

## ---- Goodness-of-fit: R2 ----
cat("\n--- R² (marginal = fixed only; conditional = fixed + random) ---\n")
print(r.squaredGLMM(best_model))

## ---- Diagnostics (DHARMa) ----
cat("\n--- DHARMa diagnostics ---\n")
sim_res <- simulateResiduals(best_model, n = 1000)

# Plots (opens a series of DHARMa diagnostic plots)
plot(sim_res)

# Uniformity (KS test), dispersion, zero-inflation, outliers
cat("\nUniformity test:\n");         print(testUniformity(sim_res))
cat("\nDispersion test:\n");         print(testDispersion(sim_res))
cat("\nZero-inflation test:\n");     print(testZeroInflation(sim_res))
cat("\nOutliers test:\n");           print(testOutliers(sim_res))

## ---- Extra optional checks ----
# Check singularity (random-effect identifiability) & overall diagnostics
cat("\n--- Optional performance checks ---\n")
print(check_singularity(best_model))
print(check_overdispersion(best_model))  # should be OK with NB, but good to see

##---Model Average---##
avg_model <- model.avg(ms)
summary(avg_model)

#model average for those models with AICc less than two
top.models <- subset(ms, delta<2)
top.models

new.avg.model <- model.avg(top.models)
summary(new.avg.model)
#Results show no difference.

library(vegan)
library(tidyverse)

# --- Prepare data in R ---
# Load ASV table
asv <- read.csv("OTU_table_ECM_and_nonECM_CCA.csv") %>%
  rename(sample = X)

# Load environmental metadata
env <- read.csv("ECMASV_GLMM_OTU.csv")

# Merge by sample ID
merged <- inner_join(asv, env, by = "sample")

# Extract ASV matrix and env variables
asv_matrix <- merged %>% select(contains("OTU"))
env_data   <- merged %>% select(sample, pH, clay, sand, site)

# Remove samples with zero total abundance
keep <- rowSums(asv_matrix) > 0
asv_matrix <- asv_matrix[keep, ]
env_data   <- env_data[keep, ]

# --- Full partial CCA model ---
cca_full <- cca(asv_matrix ~ pH + clay + sand + Condition(site), data = env_data)

# Overall test
anova(cca_full)

# Test each term
anova(cca_full, by = "term")

# Test axes
anova(cca_full, by = "axis")

# --- Model selection (adjusted R²) ---
cca_null <- cca(asv_matrix ~ Condition(site), data = env_data)
cca_step <- ordiR2step(cca_null, scope = formula(cca_full),
                       direction = "both", data = env_data)

summary(cca_step)

# --- Plot selected model ---
ordiplot(cca_step, type = "n", main = "Partial CCA (selected model)")
points(cca_step, display = "sites", col = "blue", pch = 19)
text(cca_step, display = "species", cex = 0.7, col = "darkgreen")

# Fit env vectors (only pH, clay) and add significant ones
envfit_res <- envfit(cca_step, env_data[, c("pH", "clay")], permutations = 999)
plot(envfit_res, col = "red", p.max = 0.05)

# --- Save high-resolution plots ---
png("Partial_CCA_selected.png", width = 2000, height = 1600, res = 300)
ordiplot(cca_step, type = "n", main = "Partial CCA (selected model)")
points(cca_step, display = "sites", col = "blue", pch = 19)
text(cca_step, display = "species", cex = 0.7, col = "darkgreen")
plot(envfit_res, col = "red", p.max = 0.05)
dev.off()

pdf("Partial_CCA_selected.pdf", width = 8, height = 6)
ordiplot(cca_step, type = "n", main = "Partial CCA (selected model)")
points(cca_step, display = "sites", col = "blue", pch = 19)
text(cca_step, display = "species", cex = 0.7, col = "darkgreen")
plot(envfit_res, col = "red", p.max = 0.05)
dev.off()

library(tidyverse)

# --- Load data ---
env <- read.csv("ECMASV_GLMM_OTU_alphahisto.csv")

# --- Calculate mean, sd, and se of alpha per site ---
alpha_stats <- env %>%
  group_by(site) %>%
  summarise(
    mean_alpha = mean(alpha, na.rm = TRUE),
    sd_alpha   = sd(alpha, na.rm = TRUE),
    n          = n(),
    se_alpha   = sd_alpha / sqrt(n)
  )

print(alpha_stats)

ggplot(alpha_stats, aes(x = site, y = mean_alpha, fill = site)) +
  geom_col() +
  labs(title = "Average Number of ECM Fungal OTUs per Location",
       x = "Park Location", y = "Number of ECM Fungal OTUs") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")   # removes redundant legend since site is on x-axis

library(tidyverse)
library(car)       # Levene's test for homogeneity of variance
library(FSA)       # for Dunn's test (non-parametric post-hoc)

# Load data
env <- read.csv("ECMASV_GLMM_OTU.csv")

# --- Check assumptions ---
# Shapiro-Wilk test of normality per site
by(env$alpha, env$site, shapiro.test)

# Levene's test for homogeneity of variance
leveneTest(alpha ~ site, data = env)

# --- ANOVA (if assumptions OK) ---
anova_model <- aov(alpha ~ site, data = env)
summary(anova_model)

# Tukey post-hoc
TukeyHSD(anova_model)

# --- Kruskal-Wallis test (non-parametric) ---
kruskal.test(alpha ~ site, data = env)

# Post-hoc Dunn test with Bonferroni correction
dunnTest(alpha ~ site, data = env, method = "bonferroni")

# Load libraries
library(vegan)
library(tidyverse)

# --- Load data ---
asv <- read.csv("OTU_table_ECM_and_nonECM_CCA.csv") %>%
  rename(sample = X)

env <- read.csv("ECMASV_GLMM_OTU.csv")

# --- Merge datasets by sample ID ---
merged <- inner_join(asv, env, by = "sample")

# Extract ASV abundance data
asv_data <- merged %>%
  select(contains("OTU"))

# Ensure numeric and replace NA with 0
asv_data <- asv_data %>%
  mutate(across(everything(), ~ as.numeric(.))) %>%
  replace(is.na(.), 0)

# Remove samples with zero total abundance
row_sums <- rowSums(asv_data)
merged <- merged[row_sums > 0, ]
asv_data <- asv_data[row_sums > 0, ]

# Extract grouping factor
site <- merged$site

# --- ANOSIM ---
bray_dist <- vegdist(asv_data, method = "bray")

anosim_result <- anosim(bray_dist, site)

# Output
print(anosim_result)
plot(anosim_result)

