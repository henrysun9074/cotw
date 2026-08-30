library(tidyverse)
library(knitr)
library(readxl)
library(data.table)
library(vegan)
library(scales)
library(cowplot)
library(RColorBrewer)
library(ggpubr)
library(forcats)
library(ggvenn)
library(ggrepel)
library(ggforce)
library(here)

df <- read.csv(here("Cleaned data CSVs", "qc_data_PQN.csv"))
df$X <- NULL

df <- df %>%
  mutate(
    bleaching = case_when(
      bleaching == "B"  ~ "Bleached",
      bleaching == "NB" ~ "Non-Bleached",
      is.na(bleaching)  ~ "Not Applicable",
      TRUE              ~ as.character(bleaching)
    ),
    bleaching = factor(bleaching, levels = c("Bleached", "Non-Bleached", "Not Applicable")),
    
    scleractinia = if_else(host_order == "Scleractinia", "1", "0"),
    scleractinia = factor(scleractinia, levels = c("1", "0")),
    location = factor(location),
    symbiont.potential = factor(symbiont.potential),
    host_order = fct_relevel(factor(host_order), "Scleractinia"),
    host_family = factor(host_family),
    host_phylum = factor(host_phylum)
  )

## defining compound names
met_df <- read.csv(here("Cleaned data CSVs", "merged_met_plot_df.csv"))
met_df <- met_df %>%
  mutate(
    compound_class = recode(
      compound_class,
      "Carotenoids (C40, Î²-Î²)" = "Carotenoids",
      "Oxidized glycerophospholipids" = "OxPL",
      "Glycerophosphoethanolamines" = "GPEtn",
      "Neutral glycosphingolipids" = "Neutral GSL",
      "Triacylglycerols" = "TAG",
      "Diacylglycerols" = "DAG",
      "Prenyl quinone meroterpenoids" = "TQ/THQs"
    )
  )

### take top 20 compound_class by count, with Unknown forced to end
target_classes <- met_df %>%
  count(compound_class, sort = TRUE) %>%
  slice_head(n = 20) %>%
  pull(compound_class) %>%
  trimws()

target_classes <- c(
  setdiff(target_classes, "Unknown"),
  intersect(target_classes, "Unknown")
)

### colors
provided_hex <- c(
  "#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF",
  "#9467BDFF", "#8C564BFF", "#E377C2FF", "deepskyblue4", "#BCBD22FF",
  "#17BECFFF", "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF",
  "#C5B0D5FF", "#C49C94FF", "#F7B6D2FF", "#9EDAE5FF", "#DBDB8DFF",
  "#C7C7C7FF"
)

spec_colors <- setNames(provided_hex[seq_along(target_classes)], target_classes)

final_palette <- c(spec_colors, "Other" = "gray30")
ordered_levels <- c(target_classes, "Other")

origin_shapes <- c("Host" = 16, "Symbiont" = 3, "Both" = 17, "Unknown" = 8)

# read in nina glycerolipid annotation dataframe
glycero_df <- read.csv(here("Cleaned data CSVs", "glycerolipids_fa_TyCOTW.csv"))

################################################################################
## datasets to compare

### define 100% df 
met_df_scler_all <- met_df %>%
  filter(scler_ubiquity == 100)

# define 95 %ile df 
met_df_scler_95 <- met_df %>%
  filter(scler_ubiquity >= 95)

# define 90 %ile df 
met_df_scler_90 <- met_df %>%
  filter(scler_ubiquity >= 90)

# define df of non scler
met_df_nonscler <- met_df %>%
  filter(non_scler_ubiquity >= 95)

# define ML dfs
xgb_df <- met_df %>%
  filter(XGBoost_Importance > 0) %>%
  arrange(desc(XGBoost_Importance)) %>%
  mutate(metabolite = fct_reorder(metabolite, XGBoost_Importance, .desc = TRUE))

rf_df <- met_df %>%
  filter(RandomForest_Importance > 0.001) %>%
  arrange(desc(RandomForest_Importance)) %>%
  mutate(metabolite = fct_reorder(metabolite, RandomForest_Importance, .desc = TRUE))

# define core
df_scler <- df %>% filter(host_order == "Scleractinia", !is.na(host_family))

group_summary_core <- df_scler %>%
  pivot_longer(cols = starts_with("x"), names_to = "metabolite", values_to = "val") %>%
  group_by(host_family, metabolite) %>%
  summarise(present = any(val > 0, na.rm = TRUE), .groups = "drop")

total_families <- length(unique(df_scler$host_family))

core_metabolite_ids <- group_summary_core %>%
  group_by(metabolite) %>%
  summarise(n_families = sum(present)) %>%
  filter(n_families == total_families) %>%
  pull(metabolite)

core_df <- met_df %>%
  filter(metabolite %in% core_metabolite_ids)

################################################################################

# generate contingency table for different compound classes

fisher <- function(target_df, background_df, class_col = "compound_class") {
  all_classes <- unique(background_df[[class_col]])
  
  results <- lapply(all_classes, function(current_class) {
    in_target_is_class <- sum(target_df[[class_col]] == current_class, na.rm = TRUE)
    in_target_not_class <- nrow(target_df) - in_target_is_class
    remainder_df <- background_df %>% filter(!(metabolite %in% target_df$metabolite))
    in_rem_is_class <- sum(remainder_df[[class_col]] == current_class, na.rm = TRUE)
    in_rem_not_class <- nrow(remainder_df) - in_rem_is_class
    contingency_matrix <- matrix(c(in_target_is_class, in_target_not_class, 
                                   in_rem_is_class, in_rem_not_class), 
                                 nrow = 2, byrow = TRUE)
    test <- fisher.test(contingency_matrix)
    data.frame(
      compound_class = current_class,
      count_in_target = in_target_is_class,
      count_in_total = in_rem_is_class,
      not_in_target = in_target_not_class,
      not_in_total = in_rem_not_class,
      odds_ratio = test$estimate,
      p_value = test$p.value,
      stringsAsFactors = FALSE
    )
  })
  
  bind_rows(results) %>%
    mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
    arrange(p_value)
}

enrichment_results <- list(
  xgb   = fisher(xgb_df, met_df),
  xgb_origin = fisher(xgb_df, met_df, "refined_origin"),
  rf    = fisher(rf_df, met_df),
  rf_origin = fisher(rf_df, met_df, "refined_origin"),
  core  = fisher(core_df, met_df),
  scl90 = fisher(met_df_scler_90, met_df)
)

# compound class enrichment in different datasets
significant_xgb_origin <- enrichment_results$xgb_origin %>% 
  filter(p_value < 0.05, odds_ratio > 1)
significant_rf_origin <- enrichment_results$rf_origin %>% 
  filter(p_value < 0.05, odds_ratio > 1)

significant_xgb <- enrichment_results$xgb %>% 
  filter(p_value < 0.05, odds_ratio > 1)
significant_xgb$dataset <- "xgboost"
significant_rf <- enrichment_results$rf %>% 
  filter(p_value < 0.05, odds_ratio > 1)
significant_rf$dataset <- "rf"

significant_core <- enrichment_results$core %>% 
  filter(p_value < 0.05, odds_ratio > 1)
significant_core$dataset <- "host_family"

significant_scler90 <- enrichment_results$scl90 %>% 
  filter(p_value < 0.05, odds_ratio > 1)
significant_scler90$dataset <- "scler_90"


################################################################################

## check for overlaps between dfs -- use 90%ile scleractinia

list_of_metabolites_xg <- list(
  Scler_90  = met_df_scler_90$metabolite,
  XGBoost   = xgb_df$metabolite,
  Core_Family = core_df$metabolite
)

common_metabolites_xg <- Reduce(intersect, list_of_metabolites_xg)

super_core_table_xg <- met_df %>%
  filter(metabolite %in% common_metabolites_xg) %>%
  select(metabolite, compound_class, refined_origin, XGBoost_Importance, RandomForest_Importance,
         scler_ubiquity, non_scler_ubiquity, total_ubiquity) %>%
  mutate(ubiquity_diff = scler_ubiquity - non_scler_ubiquity) %>%
  arrange(desc(XGBoost_Importance))

###########################################

list_of_metabolites_rf <- list(
  Scler_90  = met_df_scler_90$metabolite,
  RandomForest = rf_df$metabolite,
  Core_Family = core_df$metabolite
)

common_metabolites_rf <- Reduce(intersect, list_of_metabolites_rf)

# Create the summary table with metadata
super_core_table_rf <- met_df %>%
  filter(metabolite %in% common_metabolites_rf) %>%
  select(metabolite, compound_class, refined_origin, XGBoost_Importance, RandomForest_Importance,
         scler_ubiquity, non_scler_ubiquity, total_ubiquity) %>%
  mutate(ubiquity_diff = scler_ubiquity - non_scler_ubiquity) %>%
  arrange(desc(XGBoost_Importance))

# converged core from both RF and XGB
super_core_table <- rbind(super_core_table_rf, super_core_table_xg)
super_core_table <- unique(super_core_table)

class_counts <- super_core_table %>%
  summarise(
    DAG = sum(compound_class == "DAG"),
    MADAG = sum(compound_class == "MADAG"),
    TAG = sum(compound_class == "TAG"),
    Ceramide            = sum(compound_class == "Ceramides"),
    Unknown             = sum(compound_class == "Unknown")
  )
class_counts

significant_supercore <- fisher(super_core_table, met_df )%>% 
  filter(p_value < 0.05, odds_ratio > 1)
significant_supercore$dataset <- "supercore"

# combined dataset of compound classes enriched in Scleractinia across different comparisons
fishers_csv <- rbind(significant_core, significant_rf, significant_xgb, 
                     significant_scler90, significant_supercore)

#########################################
# TAG/DAG/MADAG of interest 
tag_core <- super_core_table %>%
  filter(compound_class == "TAG" | compound_class == "DAG" | compound_class == "MADAG")

mean(tag_core$ubiquity_diff)
# on average these compounds 70.3% more ubiquitous in Scleractinia than outgroups
tag_core_vec <- tag_core$metabolite

# ceramides of interest 
ceramide_core <- super_core_table %>%
  filter(compound_class == "Ceramides")

mean(ceramide_core$ubiquity_diff)
# on average these compounds 37.91% more ubiquitous in Scleractinia than outgroups
ceramide_core_vec <- ceramide_core$metabolite

################################################################################

# glycerolipid df lookup
glycero_df_core <- glycero_df %>%
  filter(glycero_df$metabolite %in% tag_core_vec)

glycero_df_core %>%
  group_by(metabolite) %>%
  group_walk(~ {
    cat("\n", paste(rep("-", 30), collapse = ""), "\n")
    cat("Metabolite:", .y$metabolite, "\n")
    cat(paste(rep("-", 30), collapse = ""), "\n")
    
    sub_table <- .x %>%
      select(fatty_acid, compound_name, molecular_formula)
    
    print(as.data.frame(sub_table))
  })

glycero_df_summary <- glycero_df_core %>%
  left_join(
    met_df %>% select(metabolite, scler_ubiquity, non_scler_ubiquity), 
    by = "metabolite"
  ) %>%
  group_by(metabolite) %>%
  summarise(
    FA_Composition = paste(fatty_acid, collapse = " | "),
    Compound_Names = paste(unique(compound_name), collapse = " | "),
    Formulas       = paste(unique(molecular_formula), collapse = " | "),
      Scleractinian_Ubiquity     = paste(unique(scler_ubiquity), collapse = " | "),
    Non_Scleractinian_Ubiquity = paste(unique(non_scler_ubiquity), collapse = " | "),
    .groups = "drop"
  )
print(glycero_df_summary)

### updated table 1 for paper 
met_df_ranked <- met_df %>%
  mutate(
    XGBoost_Rank      = min_rank(desc(XGBoost_Importance)),
    RandomForest_Rank = min_rank(desc(RandomForest_Importance))
  )

glycero_core_paper_summary <- glycero_df_core %>%
  left_join(
    met_df_ranked %>% 
      select(
        metabolite, 
        XGBoost_Rank, 
        RandomForest_Rank, 
        refined_origin
      ), 
    by = "metabolite"
  ) %>%
  group_by(metabolite) %>%
  summarise(
    Compound_Names    = paste(unique(compound_name), collapse = " | "),
    FA_Composition    = paste(fatty_acid, collapse = " | "),
    XGBoost_Rank      = paste(unique(XGBoost_Rank), collapse = " | "),
    RandomForest_Rank = paste(unique(RandomForest_Rank), collapse = " | "),
    Refined_Origin    = paste(unique(refined_origin), collapse = " | "),
    .groups = "drop"
  )
print(glycero_core_paper_summary)

################################################################################
# fisher's exact tests to compare fatty acid chain lengths

test_df <- glycero_df %>%
  mutate(is_core = if_else(metabolite %in% tag_core_vec, "Core", "Non-Core"))

# extract FA chain lengths
extract_lengths <- function(x) {
  unique(as.numeric(unlist(str_extract_all(x, "\\d+(?=:)"))))
}

# create a long-format dataframe where every metabolite-FA combo is a row
fa_analysis <- test_df %>%
  mutate(lengths = map(fatty_acid, extract_lengths)) %>%
  unnest(lengths) %>%
  filter(lengths %in% c(12, 14, 16, 18, 20, 22, 24, 26, 28)) %>%
  select(metabolite, is_core, lengths) %>%
  distinct()

target_lengths <- c(12, 14, 16, 18, 20, 22, 24, 26, 28)

fisher_results <- map_df(target_lengths, function(len) {
  
  # Create contingency table
  tab <- matrix(c(
    sum(fa_analysis$is_core == "Core" & fa_analysis$lengths == len),
    sum(fa_analysis$is_core == "Non-Core" & fa_analysis$lengths == len),
    length(unique(test_df$metabolite[test_df$is_core == "Core"])) - 
      sum(fa_analysis$is_core == "Core" & fa_analysis$lengths == len),
    length(unique(test_df$metabolite[test_df$is_core == "Non-Core"])) - 
      sum(fa_analysis$is_core == "Non-Core" & fa_analysis$lengths == len)
  ), nrow = 2, byrow = TRUE)
  
  test <- fisher.test(tab)
  tibble(
    FA_Length = as.character(len),
    p_value = test$p.value,
    odds_ratio = test$estimate,
    conf_low = test$conf.int[1],
    conf_high = test$conf.int[2]
  )
})

fisher_results <- fisher_results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)
print(fisher_results)

#########################################

# function to generate contingency table and run Fisher's test for Fatty acids
test_fa_length <- function(target_length, core_vec, full_df) {
  
  all_metabolites <- unique(full_df$metabolite)
  relevant_core <- intersect(core_vec, all_metabolites)
  relevant_non_core <- setdiff(all_metabolites, relevant_core)
  
  has_target <- full_df %>%
    mutate(lengths = map(fatty_acid, extract_lengths)) %>%
    unnest(lengths) %>%
    filter(lengths == target_length) %>%
    pull(metabolite) %>%
    unique()
  
  tab <- matrix(c(
    sum(relevant_core %in% has_target),      # Core with FA
    sum(relevant_non_core %in% has_target),  # Non-Core with FA
    sum(!(relevant_core %in% has_target)),   # Core without FA
    sum(!(relevant_non_core %in% has_target)) # Non-Core without FA
  ), nrow = 2, byrow = TRUE)
  
  rownames(tab) <- c(paste0("C", target_length, "_Present"), 
                     paste0("C", target_length, "_Absent"))
  colnames(tab) <- c("Core", "Non-Core")
  
  cat("\n--- Contingency Table for FA Length", target_length, "---\n")
  print(tab)
  
  res <- fisher.test(tab)
  cat("\np-value:", res$p.value, "\n")
  return(res)
}
test_fa_length(22, tag_core_vec, glycero_df)
