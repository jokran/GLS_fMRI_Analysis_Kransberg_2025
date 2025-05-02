# grid_analysis.R
# Analysis script for Kransberg et al. (2025) - GLS fMRI Study
# Author: Jonas Kransberg
# Date last modified: 02.05.2025
# Author contact information: jonas.kransberg@psykologi.uio.no

# 0. Load required libraries
library(tidyverse)
library(lmerTest)

# 1. Generate dummy data
set.seed(123)
# 1.1 Create participant-level info for manual segmentation
n_manual <- 110
manual_ids <- factor(1:n_manual)
participant_manual <- tibble(
  subject_id             = manual_ids,
  visit_age              = runif(n_manual, 16, 80),
  subject_sex            = factor(sample(c("Male","Female"), n_manual, replace = TRUE)),
  total_correct_combined = sample(1:7, n_manual, replace = TRUE)
) %>%
  mutate(
    age_group = factor(ifelse(visit_age < 40, "Under40", "Above40"),
                       levels = c("Under40","Above40"))
  )

# 1.2 Expand to row-level for manual segmentation
manual <- expand_grid(
  participant_manual,
  segmentation = "manual_segmentation",
  xfoldsym     = c(1,5,6,7),
  roi          = factor(c("rh_erc","lh_erc")),
  smoothing    = c(0,4)
) %>%
  mutate(
    beta_gridcode_mean_combined = rnorm(n(), 0, 0.1),
    rayleigh_z        = runif(n(), 0, 1),        # spatial stability metric
    temporally_stable = runif(n(), 0, 100)       # temporal stability metric (percent)
  )

# 1.3 Create participant-level info for automatic segmentation
n_auto <- 207
auto_ids <- factor(1:n_auto)
participant_auto <- tibble(
  subject_id             = auto_ids,
  visit_age              = runif(n_auto, 18.2, 78.6),
  subject_sex            = factor(sample(c("Male","Female"), n_auto, replace = TRUE)),
  total_correct_combined = sample(1:7, n_auto, replace = TRUE)
) %>%
  mutate(
    age_group = factor(ifelse(visit_age < 40, "Under40", "Above40"),
                       levels = c("Under40","Above40"))
  )

automatic <- expand_grid(
  participant_auto,
  segmentation = "automatic_segmentation",
  xfoldsym     = c(1,5,6,7),
  roi          = factor(c("rh_erc","lh_erc")),
  smoothing    = c(0,4)
) %>%
  mutate(
    beta_gridcode_mean_combined = rnorm(n(), 0, 0.1)
  )

# Combine datasets
combined <- bind_rows(manual, automatic)

# 2. Grid-Like Signal (GLS) Analyses
## 2.1 6-Fold GLS vs Zero
for(region in c("rh_erc","lh_erc")){
  cat("6-fold GLS vs zero in", region, "\n")
  print(t.test(filter(manual, xfoldsym==6, roi==region, smoothing == 0)$beta_gridcode_mean_combined,
               mu=0, alternative = "greater", na.rm = TRUE))
}

## 2.2 Control Symmetry Checks (5- and 7-fold)
for(sym in c(5,7)){
  for(region in c("rh_erc","lh_erc")){
    cat(sym, "-fold GLS vs zero in", region, "\n")
    print(t.test(filter(manual, xfoldsym==sym, roi==region, smoothing == 0)$beta_gridcode_mean_combined,
                 mu=0, alternative = "greater", na.rm = TRUE))
  }
}

## 2.3 Age-Group Stratification (<40 vs ≥40)
### Younger group (<40)
for(sym in c(6,5,7)){
  for(region in c("rh_erc","lh_erc")){
    cat("Younger (<40):", sym, "-fold GLS vs zero in", region, "\n")
    data_y <- filter(manual, age_group=="Under40", xfoldsym==sym, roi==region, smoothing == 0)
    print(t.test(data_y$beta_gridcode_mean_combined, mu=0, alternative = "greater", na.rm = TRUE))
  }
}

### Older group (≥40)
for(sym in c(6,5,7)){
  for(region in c("rh_erc","lh_erc")){
    cat("Older (≥40):", sym, "-fold GLS vs zero in", region, "\n")
    data_o <- filter(manual, age_group=="Above40", xfoldsym==sym, roi==region, smoothing == 0)
    print(t.test(data_o$beta_gridcode_mean_combined, mu=0, alternative = "greater", na.rm = TRUE))
  }
}

## 2.4 Between-Group Comparison for 6-Fold
for(region in c("rh_erc","lh_erc")){
  cat("Between-group Welch t-test for 6-fold in", region, "\n")
  dt <- filter(manual, xfoldsym==6, roi==region, smoothing == 0)
  print(t.test(beta_gridcode_mean_combined ~ age_group, data=dt, alternative = "greater"))
}

## 2.5 Continuous Age Effects (LME)
model_6f <- lmer(beta_gridcode_mean_combined ~ visit_age * roi + subject_sex +
                   (1|subject_id), data = filter(manual, xfoldsym==6, smoothing == 0))
summary(model_6f)

# 3.  Temporal and Spatial Stability Analyses
temp <- manual %>%
  filter(xfoldsym == 6, smoothing == 0) 
## 3.1 Temporal Stability Analyses
for(region in c("rh_erc","lh_erc")){
  df <- filter(temp, roi == region)
  cat("Temporal stability one-sample t-test in", region, "(H0: mean = 50%)\n")
  print(t.test(df$temporally_stable, mu = 50))
  cat("Temporal stability age-group comparison in", region, "\n")
  print(t.test(temporally_stable ~ age_group, data = df))
}

## 3.2 Spatial Stability Analyses
# Using Rayleigh's Z, test age-group differences
for(region in c("rh_erc","lh_erc")){
  df <- filter(temp, roi == region)
  cat("Spatial stability age-group t-test for rayleigh_z in", region, "\n")
  print(t.test(rayleigh_z ~ age_group, data = df))
}


# 4 Control Analyses: Smoothing
## 4.1 6-Fold GLS vs Zero
for(region in c("rh_erc","lh_erc")){
  cat("6-fold GLS vs zero in", region, "\n")
  print(t.test(filter(manual, xfoldsym==6, roi==region, smoothing == 4)$beta_gridcode_mean_combined,
               mu=0))
}

## 4.2 Control Symmetry Checks (5- and 7-fold)
for(sym in c(5,7)){
  for(region in c("rh_erc","lh_erc")){
    cat(sym, "-fold GLS vs zero in", region, "\n")
    print(t.test(filter(manual, xfoldsym==sym, roi==region, smoothing == 4)$beta_gridcode_mean_combined,
                 mu=0))
  }
}

## 4.3 Age-Group Stratification (<40 vs ≥40)
### Younger group (<40)
for(sym in c(6,5,7)){
  for(region in c("rh_erc","lh_erc")){
    cat("Younger (<40):", sym, "-fold GLS vs zero in", region, "\n")
    data_y <- filter(manual, age_group=="Under40", xfoldsym==sym, roi==region, smoothing == 4)
    print(t.test(data_y$beta_gridcode_mean_combined, mu=0))
  }
}

### Older group (≥40)
for(sym in c(6,5,7)){
  for(region in c("rh_erc","lh_erc")){
    cat("Older (≥40):", sym, "-fold GLS vs zero in", region, "\n")
    data_o <- filter(manual, age_group=="Above40", xfoldsym==sym, roi==region, smoothing == 4)
    print(t.test(data_o$beta_gridcode_mean_combined, mu=0))
  }
}

## 4.4 Between-Group Comparison for 6-Fold
for(region in c("rh_erc","lh_erc")){
  cat("Between-group Welch t-test for 6-fold in", region, "\n")
  dt <- filter(manual, xfoldsym==6, roi==region, smoothing == 4)
  print(t.test(beta_gridcode_mean_combined ~ age_group, data=dt))
}

## 4.5 Correlate GLS by smoothing type
manualsmoothing_rh <- manual %>%
  filter(xfoldsym == 6, roi == "rh_erc") %>% 
  select(subject_id, smoothing, beta_gridcode_mean_combined) %>%
  pivot_wider(names_from = smoothing, 
              values_from = beta_gridcode_mean_combined,
              names_prefix = "smoothing_")


cor.test(manualsmoothing_rh$smoothing_0, manualsmoothing_rh$smoothing_4)


manualsmoothing_lh <- manual %>% 
filter(xfoldsym == 6, roi == "lh_erc") %>% 
  select(subject_id, smoothing, beta_gridcode_mean_combined) %>%
  pivot_wider(names_from = smoothing, 
              values_from = beta_gridcode_mean_combined,
              names_prefix = "smoothing_")

cor.test(manualsmoothing_lh$smoothing_0, manualsmoothing_lh$smoothing_4)


# 5 Control Analyses: Segmentation
#Re-do primary analyses with automatically segmented EC masks. 
## 5.1 6-Fold GLS vs Zero
for(region in c("rh_erc","lh_erc")){
  cat("6-fold GLS vs zero in", region, "\n")
  print(t.test(filter(automatic, xfoldsym==6, roi==region, smoothing == 0)$beta_gridcode_mean_combined,
               mu=0))
}

## 5.2 Control Symmetry Checks (5- and 7-fold)
for(sym in c(5,7)){
  for(region in c("rh_erc","lh_erc")){
    cat(sym, "-fold GLS vs zero in", region, "\n")
    print(t.test(filter(automatic, xfoldsym==sym, roi==region, smoothing == 0)$beta_gridcode_mean_combined,
                 mu=0))
  }
}

## 5.3 Age-Group Stratification (<40 vs ≥40)
### Younger group (<40)
for(sym in c(6,5,7)){
  for(region in c("rh_erc","lh_erc")){
    cat("Younger (<40):", sym, "-fold GLS vs zero in", region, "\n")
    data_y <- filter(automatic, age_group=="Under40", xfoldsym==sym, roi==region, smoothing == 0)
    print(t.test(data_y$beta_gridcode_mean_combined, mu=0))
  }
}

### Older group (≥40)
for(sym in c(6,5,7)){
  for(region in c("rh_erc","lh_erc")){
    cat("Older (≥40):", sym, "-fold GLS vs zero in", region, "\n")
    data_o <- filter(automatic, age_group=="Above40", xfoldsym==sym, roi==region, smoothing == 0)
    print(t.test(data_o$beta_gridcode_mean_combined, mu=0))
  }
}

## 5.4 Between-Group Comparison for 6-Fold
for(region in c("rh_erc","lh_erc")){
  cat("Between-group Welch t-test for 6-fold in", region, "\n")
  dt <- filter(automatic, xfoldsym==6, roi==region, smoothing == 0)
  print(t.test(beta_gridcode_mean_combined ~ age_group, data=dt))
}

## 5.5 Correlate GLS by smoothing type
# Right hemisphere
combinedcorrh <- combined %>%
  filter(xfoldsym == 6, roi == "rh_erc", smoothing == 0) %>% 
  select(subject_id, segmentation, beta_gridcode_mean_combined) %>%
  pivot_wider(names_from = segmentation, 
              values_from = beta_gridcode_mean_combined,
              names_prefix = "seg_")


cor.test(combinedcorrh$seg_manual_segmentation, combinedcorrh$seg_automatic_segmentation)

# Left hemisphere
combinedcorlh <- combined %>%
  filter(xfoldsym == 6, roi == "lh_erc", smoothing == 0) %>% 
  select(subject_id, segmentation, beta_gridcode_mean_combined) %>%
  pivot_wider(names_from = segmentation, 
              values_from = beta_gridcode_mean_combined,
              names_prefix = "seg_")


cor.test(combinedcorlh$seg_manual_segmentation, combinedcorlh$seg_automatic_segmentation)


# 6 High-Performance Participants Analysis
## 6. 1 Behavioral Performance Ceiling Effects
behavrh <- manual %>% filter(xfoldsym == 6, roi == "rh_erc", smoothing == 0)
behavlh <- manual %>% filter(xfoldsym == 6, roi == "lh_erc", smoothing == 0)
table(behavrh$age_group, behavrh$total_correct_combined >= 6)

## 6.2 Age-Group Difference in Task Scores
print(t.test(total_correct_combined ~ age_group, data=behavrh))


## 6.3 GLS After Excluding Low Scorers (<=6)
high_perf <- filter(manual, total_correct_combined >= 6)

## 6.4 6-Fold GLS vs Zero
for(region in c("rh_erc","lh_erc")){
  cat("6-fold GLS vs zero in", region, "\n")
  print(t.test(filter(high_perf, xfoldsym==6, roi==region, smoothing == 0)$beta_gridcode_mean_combined,
               mu=0))
}

## 6.5 Control Symmetry Checks (5- and 7-fold)
for(sym in c(5,7)){
  for(region in c("rh_erc","lh_erc")){
    cat(sym, "-fold GLS vs zero in", region, "\n")
    print(t.test(filter(high_perf, xfoldsym==sym, roi==region, smoothing == 0)$beta_gridcode_mean_combined,
                 mu=0))
  }
}

## 6.6 Age-Group Stratification (<40 vs ≥40)
### Younger group (<40)
for(sym in c(6,5,7)){
  for(region in c("rh_erc","lh_erc")){
    cat("Younger (<40):", sym, "-fold GLS vs zero in", region, "\n")
    data_y <- filter(high_perf, age_group=="Under40", xfoldsym==sym, roi==region, smoothing == 0)
    print(t.test(data_y$beta_gridcode_mean_combined, mu=0))
  }
}

### Older group (≥40)
for(sym in c(6,5,7)){
  for(region in c("rh_erc","lh_erc")){
    cat("Older (≥40):", sym, "-fold GLS vs zero in", region, "\n")
    data_o <- filter(high_perf, age_group=="Above40", xfoldsym==sym, roi==region, smoothing == 0)
    print(t.test(data_o$beta_gridcode_mean_combined, mu=0))
  }
}

## 6.7 Between-Group Comparison for 6-Fold
for(region in c("rh_erc","lh_erc")){
  cat("Between-group Welch t-test for 6-fold in", region, "\n")
  dt <- filter(high_perf, xfoldsym==6, roi==region, smoothing == 0)
  print(t.test(beta_gridcode_mean_combined ~ age_group, data=dt))
}


## 6.8 Correlation Analyses
# Behavioural task scores and Age
print(cor.test(behavrh$visit_age, behavrh$total_correct_combined))

# Behavioural task scores and GLS magnitude
# Right hemisphere
print(cor.test(behavrh$beta_gridcode_mean_combined, behavrh$total_correct_combined))

# Left hemisphere
print(cor.test(behavlh$beta_gridcode_mean_combined, behavlh$total_correct_combined))


# 7. Unimodal (1-Fold) Signal Analyses
## 7.1 1-Fold Signal vs Behavior
filtered <- manual %>% filter(smoothing == 0)
#rh_cor
print(cor.test(filter(filtered, xfoldsym==1, roi == "rh_erc")$beta_gridcode_mean_combined,
               filter(filtered, xfoldsym==1, roi == "rh_erc")$total_correct_combined))

#lh_cor
print(cor.test(filter(filtered, xfoldsym==1, roi == "lh_erc")$beta_gridcode_mean_combined,
               filter(filtered, xfoldsym==1, roi == "lh_erc")$total_correct_combined))

## 7.2 1-Fold Signal
lme_1f <- lmer(beta_gridcode_mean_combined ~ visit_age * roi + subject_sex +
                 (1|subject_id), data = filter(filtered, xfoldsym==1))
summary(lme_1f)

## 7.3 1-Fold Predicting 7-Fold & 6-Fold
wide <- manual %>%
  filter(smoothing == 0) %>%
  select(subject_id, roi, xfoldsym, beta_gridcode_mean_combined, visit_age) %>%
  pivot_wider(
    names_from = xfoldsym,
    names_prefix = "fold_",
    values_from = beta_gridcode_mean_combined
  )

lme_1v7 <- lmer(fold_7 ~
                  fold_1 * roi + visit_age +
                  (1|subject_id), data=wide)
summary(lme_1v7)

lme_1v6 <- lmer(fold_6 ~
                  fold_1 * roi + visit_age +
                  (1|subject_id), data=wide)
summary(lme_1v6)
