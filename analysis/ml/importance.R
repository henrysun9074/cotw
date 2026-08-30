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
library(caret)
library(tibble)
library(stringr)
library(RColorBrewer)
library(ggrepel)
library(here)

# read in data
df <- read.csv(here("Cleaned data CSVs", "qc_data_PQN.csv"))
df$X <- NULL

met_df <- read.csv(here("Cleaned data CSVs", "merged_met_plot_df.csv"))

present_metabolites <- df %>% 
  select(starts_with("x")) %>% 
  colnames()

met_df <- met_df %>%
  filter(met_df$metabolite %in% present_metabolites)

cols_bleaching <- c(
  "Bleached" = "#FF847CFF", 
  "Non-Bleached" = "#019875FF", 
  "Not Applicable" = "#D3D3D3")

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

# color palettes
cols_location <-c("#002594FF", "#E0B2CDFF", "#54C4E3FF", "#F3AA4FFF")
cols_symbiont  <- c("#D84D16FF", "#FFF800FF", "#8FDA04FF")
cols_phylum <- c("#24492EFF", "#015B58FF", "#2C6184FF", "#59629BFF", "#89689DFF", "#BA7999FF", "#E69B99FF")
cols_sclero    <- c("1" = "#DE7862FF", "0" = "#D8AF39FF")

##################################################

# for compound class - custom spectral library

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

target_classes <- met_df %>%
  count(compound_class, sort = TRUE) %>%
  slice_head(n = 20) %>%
  pull(compound_class) %>%
  trimws()

target_classes <- c(
  setdiff(target_classes, "Unknown"),
  intersect(target_classes, "Unknown")
)

provided_hex <- c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF", "#D62728FF", 
                  "#9467BDFF", "#8C564BFF", "#E377C2FF", "deepskyblue4", "#BCBD22FF", 
                  "#17BECFFF", "#AEC7E8FF", "#FFBB78FF", "#98DF8AFF", "#FF9896FF", 
                  "#C5B0D5FF", "#C49C94FF", "#F7B6D2FF", "#9EDAE5FF", "#DBDB8DFF", 
                  "#C7C7C7FF")

spec_colors <- setNames(provided_hex, target_classes)

final_palette <- c(spec_colors, "Other" = "gray30")
ordered_levels <- c(target_classes, "Other")

process_importance_data <- function(df) {
  df %>%
    mutate(compound_class = trimws(as.character(compound_class))) %>%
    mutate(display_class = if_else(compound_class %in% names(final_palette), 
                                   compound_class, 
                                   "Other")) %>%
    mutate(display_class = fct_relevel(factor(display_class), "Other", after = Inf)) 
}

met_plot_df <- process_importance_data(met_df)

ordered_levels <- c(target_classes, "Other")
met_plot_df$display_class <- factor(met_plot_df$display_class, levels = ordered_levels)


##################################################

origin_shapes <- c("Host" = 16, "Symbiont" = 3, "Both" = 17, "Unknown" = 8)

#################################################################################

# total ubiquity 
sum(met_plot_df$total_ubiquity == 100, na.rm = TRUE)
# 10 ubiquitous compounds found in all samples
sum(met_plot_df$total_ubiquity > 90, na.rm = TRUE)
# 413 compounds found in 90% or more of all samples

# Scleractinia summary stats
sum(met_plot_df$scler_ubiquity == 100, na.rm = TRUE)
# 59 ubiquitous compounds detected across all Scleractinian samples
sum(met_plot_df$scler_ubiquity > 90, na.rm = TRUE)
# 854 compounds found in over 90% of Scleractinian samples
sum(met_plot_df$scler_ubiquity > 80, na.rm = TRUE)
# 2157 compounds found in over 80% of Scleractinian samples
sum(met_plot_df$scler_ubiquity == 0, na.rm = TRUE)
# 841 compounds never found in Scleractinians

# Non-Scleractinia summary stats
sum(met_plot_df$non_scler_ubiquity == 100, na.rm = TRUE)
# 146 compounds found across all non Scleractinian samples
sum(met_plot_df$non_scler_ubiquity > 90, na.rm = TRUE)
# 425 compounds found in over 90% of non Scleractinian samples
sum(met_plot_df$non_scler_ubiquity > 80, na.rm = TRUE)
# 754 compounds found in over 90% of non Scleractinian samples
sum(met_plot_df$non_scler_ubiquity == 0, na.rm = TRUE)
# 2105 compounds never found in non Scleractinian samples


#################################################################################

# for figure S4

met_presence_long <- df %>%
  pivot_longer(
    cols = starts_with("x"),
    names_to = "metabolite",
    values_to = "value"
  ) %>%
  mutate(present = value > 0)

ubiquity_overall <- met_presence_long %>%
  group_by(metabolite) %>%
  summarise(ubiquity_all = mean(present) * 100, .groups = "drop")

ubiquity_corals <- met_presence_long %>%
  filter(host_order == "Scleractinia") %>%
  group_by(metabolite) %>%
  summarise(ubiquity_coral = mean(present) * 100, .groups = "drop")

met_summary <- ubiquity_overall %>%
  left_join(ubiquity_corals, by = "metabolite") %>%
  mutate(ubiquity_coral = ifelse(is.na(ubiquity_coral), 0, ubiquity_coral))

coral_present <- met_presence_long %>%
  filter(host_order == "Scleractinia", present) %>%
  distinct(metabolite)

noncoral_present <- met_presence_long %>%
  filter(host_order != "Scleractinia", present) %>%
  distinct(metabolite)

coral_only <- setdiff(coral_present$metabolite, noncoral_present$metabolite)

met_summary <- met_summary %>%
  left_join(
    met_presence_long %>%
      group_by(metabolite) %>%
      summarise(avg_abundance = mean(value, na.rm = TRUE), .groups = "drop"),
    by = "metabolite"
  ) %>%
  mutate(category = ifelse(metabolite %in% coral_only, "Coral-only", "Other"))
met_plot_df$display_class <- factor(met_plot_df$display_class, levels = ordered_levels)

x_vline_pos <- met_summary %>%
  filter(category == "Coral-only") %>%
  pull(ubiquity_all) %>%
  { if(length(.) == 0) NA_real_ else max(., na.rm = TRUE) }

class_mapping <- met_plot_df %>%
  select(metabolite, display_class) %>%
  distinct()

met_plot_df$display_class <- factor(met_plot_df$display_class, levels = ordered_levels)
class_mapping <- met_plot_df %>%
  select(metabolite, display_class, compound_class) %>%
  distinct()

met_summary_classed <- met_summary %>%
  left_join(class_mapping, by = "metabolite") %>%
  mutate(display_class = replace_na(display_class, "Other"))

met_summary_classed <- met_summary_classed %>%
  left_join(met_df %>% select(metabolite, refined_origin), by = "metabolite") 

met_summary_classed$refined_origin <- factor(
  met_summary_classed$refined_origin, 
  levels = c("Host", "Symbiont", "Both", "Unknown")
)

### Figure S4A
pa <- ggplot(met_summary_classed, aes(x = ubiquity_all, y = avg_abundance)) +
  geom_point(
    aes(color = display_class, shape = refined_origin),
    size = 3,      
    stroke = 0.8,  
    alpha = 0.85
  ) +
  scale_color_manual(values = final_palette) +
  scale_shape_manual(values = origin_shapes) +
  scale_y_continuous(
    labels = label_number(scale_cut = cut_short_scale())
  ) + 
  scale_x_continuous(
    limits = c(0, 100), 
    breaks = seq(0, 100, by = 20)
  ) +
  labs(
    x = "Ubiquity",
    y = "Average Abundance",
    color = "Compound Class", # Updated label
    shape = "Metabolite Origin"
  ) + 
  theme_pubr() + 
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    axis.line         = element_line(color = "black"),
    legend.position   = "right",
    legend.text       = element_text(size = 9),
    plot.title        = element_text(hjust = 0.5, face = "bold", size = 16)
  )
pa

#################################################################################

### Figure S4B

## code to make this plot using only Scleractinian samples and metabolites (regardless of origin)
## that are only found in Scleractinia
host_only_mets <- met_df %>%
  filter(scler_ubiquity != 0) %>%
  pull(metabolite)

met_presence_scler_host <- df %>%
  filter(scleractinia == 1) %>% # Only Scleractinian samples
  select(sample, all_of(intersect(names(.), host_only_mets))) %>% # Only Host metabolites
  pivot_longer(
    cols = -sample,
    names_to = "metabolite",
    values_to = "value"
  ) %>%
  mutate(present = value > 0)

met_summary_scler <- met_presence_scler_host %>%
  group_by(metabolite) %>%
  summarise(
    ubiquity_all = mean(present) * 100,
    avg_abundance = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

met_summary_classed_scler <- met_summary_scler %>%
  left_join(class_mapping, by = "metabolite") %>%
  left_join(met_df %>% select(metabolite, refined_origin), by = "metabolite") %>%
  mutate(
    display_class = trimws(as.character(display_class)),
    display_class = if_else(
      is.na(display_class) | !(display_class %in% names(final_palette)),
      "Other",
      display_class
    ),
    display_class = factor(display_class, levels = names(final_palette)),
    refined_origin = factor(
      refined_origin,
      levels = c("Host", "Symbiont", "Both", "Unknown")
    )
  )
pb <- ggplot(met_summary_classed_scler, aes(x = ubiquity_all, y = avg_abundance)) +
  geom_point(
    aes(color = display_class, shape = refined_origin),
    size = 3,
    stroke = 0.4,
    alpha = 0.85
  ) +
  scale_color_manual(values = final_palette) +
  scale_shape_manual(values = origin_shapes) +
  scale_y_continuous(
    labels = label_number(scale_cut = cut_short_scale())
  ) + 
  scale_x_continuous(
    limits = c(0, 100), 
    breaks = seq(0, 100, by = 20)
  ) +
  labs(
    x = "Scleractinian Ubiquity",
    y = "Average Abundance",
    color = "Compound Class",
    shape = "Metabolite Origin",
  ) + 
  theme_pubr() + 
  theme(
    plot.background   = element_rect(fill = "white", color = NA),
    panel.background  = element_rect(fill = "white", color = NA),
    axis.line         = element_line(color = "black"),
    legend.position   = "right",
    legend.text       = element_text(size = 9),
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 4, shape = 21))
  )
pb

#################################################################################

xgb_plot_data <- met_plot_df %>%
  filter(XGBoost_Importance > 0) %>%
  arrange(desc(XGBoost_Importance)) %>%
  mutate(metabolite = fct_reorder(metabolite, XGBoost_Importance, .desc = TRUE))
# 61 metabolites with feature importance > 0, all over 0.0005
# 57 with feature importance > 0.001

rf_plot_data <- met_plot_df %>%
  filter(RandomForest_Importance > 0.001) %>%
  arrange(desc(RandomForest_Importance)) %>%
  mutate(metabolite = fct_reorder(metabolite, RandomForest_Importance, .desc = TRUE))
# 2765 metabolites with feature importance > 0
# 267 with feature importance > 0.001


xgb_origin_summary <- xgb_plot_data %>%
  count(refined_origin, name = "n") %>%
  mutate(percentage = round(n / sum(n) * 100, 1))

rf_origin_summary <- rf_plot_data %>%
  count(refined_origin, name = "n") %>%
  mutate(percentage = round(n / sum(n) * 100, 1))


#################### make feature importance plots ###########################

# Figure 5A

p1 <- ggbarplot(xgb_plot_data, x = "metabolite", y = "XGBoost_Importance",
                fill = "display_class", color = "transparent",
                xlab = "Metabolite", ylab = "XGBoost Feature Importance") +
  theme_pubr() +
  scale_fill_manual(values = final_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
  scale_x_discrete(expand = expansion(add = c(1, 0.5))) + 
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.position = "none")
p1

# Figure 5B

p2 <- ggbarplot(rf_plot_data, x = "metabolite", y = "RandomForest_Importance",
                fill = "display_class", color = "transparent",
                xlab = "Metabolite", ylab = "RF Feature Importance") +
  theme_pubr() +
  scale_fill_manual(values = final_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_x_discrete(expand = expansion(add = c(1, 0.5))) +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        legend.position = "none")
p2


# Figure 5C

top_labels <- met_plot_df %>%
  mutate(dist = sqrt(XGBoost_Importance^2 + RandomForest_Importance^2)) %>%
  arrange(desc(dist)) %>%
  slice_head(n = 20) %>%
  pull(metabolite)

origin_shapes <- c("Host" = 16, "Symbiont" = 3, "Both" = 17, "Unknown" = 8)

fit <- lm(RandomForest_Importance ~ XGBoost_Importance, data = met_plot_df)
r2 <- summary(fit)$r.squared
pval <- summary(fit)$coefficients[2, 4]
label_text <- paste0("R² = ", round(r2, 3),
                     "\n p = ", signif(pval, 3))

p3 <- ggscatter(
  met_plot_df, 
  x = "XGBoost_Importance", 
  y = "RandomForest_Importance",
  color = "display_class", 
  shape = "refined_origin", 
  palette = final_palette,
  size = 5,
  font.x = c(22, "bold"),
  font.y = c(16, "bold"),
  repel = FALSE,  
  add = "reg.line",
  add.params = list(color = "steelblue4", fill = "lightgray", linetype = "dashed", linewidth = 0.8),
  cor.coeff = TRUE, 
  cor.method = "pearson",
  xlab = "XGBoost Feature Importance", 
  ylab = "RF Feature Importance"
) +
  geom_text_repel(
    data = met_plot_df %>% filter(metabolite %in% top_labels),
    aes(label = coraldb_compound_name, color = display_class),
    size = 3,
    fontface = "italic",
    force = 20,            
    box.padding = 2,  
    point.padding = 1,
    max.overlaps = Inf
  ) +
  scale_shape_manual(values = origin_shapes) +
  theme_pubr() +
  theme(
    legend.position = "right",
    axis.title = element_text(size = 16)
  ) +
  guides(
    color = "none",
    fill = "none",
    shape = guide_legend(
      title = "Metabolite Origin",
      override.aes = list(size = 5)
    )
  )
p3

#################################################################################

# take metabolites to plot ubiquity in Scler vs non Scler
top_xgb_mets <- met_plot_df %>% arrange(desc(XGBoost_Importance)) %>% head(61) %>% pull(metabolite)
top_rf_mets  <- met_plot_df %>% arrange(desc(RandomForest_Importance)) %>% head(267) %>% pull(metabolite)

########################################
met_summary_xgb <- met_plot_df %>% filter(metabolite %in% top_xgb_mets)
met_summary_rf  <- met_plot_df %>% filter(metabolite %in% top_rf_mets)

sum(met_summary_xgb$scler_ubiquity == 100, na.rm = TRUE)
sum(met_summary_xgb$scler_ubiquity > 90, na.rm = TRUE)
sum(met_summary_xgb$scler_ubiquity > 80, na.rm = TRUE)
sum(met_summary_xgb$non_scler_ubiquity > 90, na.rm = TRUE)
sum(met_summary_xgb$non_scler_ubiquity > 80, na.rm = TRUE)
# 4/61 Scleractinian ubiquitous compounds, 25 ubiquity > 90%, 32 ubiquity > 80%
# 11/61 compounds non-Scleractinian ubiquity > 90%, 15 > 80%

sum(met_summary_rf$scler_ubiquity == 100, na.rm = TRUE)
sum(met_summary_rf$scler_ubiquity > 90, na.rm = TRUE)
sum(met_summary_rf$scler_ubiquity > 80, na.rm = TRUE)
sum(met_summary_rf$non_scler_ubiquity > 90, na.rm = TRUE)
sum(met_summary_rf$non_scler_ubiquity > 80, na.rm = TRUE)
# 6/267 Scleractinian ubiquitous compounds, 69 ubiquity > 90%, 134 ubiquity > 80%
# 60/267 compounds non-Scleractinian ubiquity > 90%, 68 > 80%

########################################
sort(table(met_summary_xgb$compound_class), decreasing = TRUE)
# 37 annotations, 24 unknown
# 12 TAG, 7 ceramides, 5 glycerophosphocolines, everything else 2 or fewer
sort(table(met_summary_rf$compound_class), decreasing = TRUE)
# 119 annotations, 148 unknown
# 48 TAG, 18 ceramides, 12 glycerophosphocholines, 7 MADAG, everything else 5 or fewer

# Figure 5D
p_xgb_ubiq_scatter <- ggplot(met_summary_xgb, aes(x = non_scler_ubiquity, y = scler_ubiquity)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray60") +
  geom_point(
    aes(color = display_class, shape = refined_origin), 
    size = 4, stroke = 1, alpha = 0.9
  ) +
  scale_color_manual(values = final_palette) +
  scale_shape_manual(values = origin_shapes) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(
    x = "Non-Scleractinian Ubiquity (%)", 
    y = "Scleractinian Ubiquity (%)", 
    color = "Compound Class",
    shape = "Metabolite Origin"
  ) +
  theme_pubr() +
  theme(
    legend.position = "none", 
    plot.title = element_text(size = 12, face = "bold")
  )
p_xgb_ubiq_scatter


# Figure 5E
p_rf_ubiq_scatter <- ggplot(met_summary_rf, aes(x = non_scler_ubiquity, y = scler_ubiquity)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray60") +
  geom_point(
    aes(color = display_class, shape = refined_origin), 
    size = 4, stroke = 1, alpha = 0.9
  ) +
  scale_color_manual(values = final_palette) +
  scale_shape_manual(values = origin_shapes) + 
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(
    x = "Non-Scleractinian Ubiquity (%)", 
    y = "Scleractinian Ubiquity (%)", 
    color = "Compound Class",
    shape = "Metabolite Origin"
  ) +
  theme_pubr() +
  theme(
    legend.position = "none", 
    plot.title = element_text(size = 12, face = "bold")
  )
p_rf_ubiq_scatter

################################################################################

# volcano plot for Figure S15

comm_long <- df %>%
  pivot_longer(
    cols = starts_with("x"), 
    names_to  = "metabolite",
    values_to = "abundance"
  )

comm_long <- comm_long %>%
  mutate(abundance = as.numeric(as.character(abundance)))

stats_data <- comm_long %>%
  select(-scleractinia) %>%                    
  left_join(df %>% select(sample, scleractinia), by = "sample") %>%
  filter(!is.na(scleractinia)) %>%
  mutate(group = if_else(as.character(scleractinia) == "1", "Scleractinia", "Other"))

# compute L2FC and p-values per metabolite
volcano_results <- stats_data %>%
  group_by(metabolite) %>%
  summarise(
    mean_scler = mean(abundance[group == "Scleractinia"], na.rm = TRUE),
    mean_other = mean(abundance[group == "Other"], na.rm = TRUE),
    p_val_raw = wilcox.test(abundance ~ group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_val_raw, method = "bonferroni"),
    log2FC = log2((mean_scler + 1) / (mean_other + 1)),
    # Use adjusted p-value for the y-axis
    neg_log_p_adj = -log10(p_adj)
  )

class_order <- c(target_classes, "Other")

plot_data_volcano <- volcano_results %>%
  inner_join(
    met_df %>% select(metabolite, compound_class, refined_origin),
    by = "metabolite"
  ) %>%
  mutate(
    compound_class = trimws(as.character(compound_class)),
    refined_origin = as.character(refined_origin),
    display_class = if_else(
      is.na(compound_class) | !(compound_class %in% names(final_palette)),
      "Other",
      compound_class
    ),
    refined_origin = if_else(is.na(refined_origin), "Unknown", refined_origin),
    display_class = factor(display_class, levels = class_order)
  )

classes <- levels(droplevels(plot_data_volcano$display_class))
class_colors <- final_palette[classes]
m <- nrow(volcano_results)   
sig_threshold <- -log10(0.05 / m)

## plot
make_volcano <- function(df, title = NULL) {
  ggplot(df, aes(x = log2FC, y = neg_log_p_adj)) +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey70") +
    geom_hline(yintercept = sig_threshold, linetype = "dashed", color = "grey70", linewidth = 0.8) +
    
    geom_point(aes(color = display_class), alpha = 0.75, size = 2.5) +
    
    scale_color_manual(
      name = "Compound Class",
      values = class_colors,
      breaks = classes,
      na.value = "gray60"
    ) +
    ylim(0, 100) +
    xlim(-25, 25) +
    labs(
      x = "log2 Fold Change",
      y = "-log10(adj. p-value)",
      title = title
    ) +
    
    theme_pubr() +
    theme(
      legend.position = "right",   # show legend now
      axis.title = element_text(size = 15),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

plot_data_xgb <- plot_data_volcano %>%
  filter(metabolite %in% met_summary_xgb$metabolite)

plot_data_rf <- plot_data_volcano %>%
  filter(metabolite %in% met_summary_rf$metabolite)

overlap_count <- length(intersect(plot_data_xgb$metabolite, plot_data_rf$metabolite)) #31

##### summmary stats xgb
dems <- plot_data_xgb %>%
  filter(p_adj < 0.05 & log2FC > 2)
class_counts <- dems %>%
  count(compound_class, sort = TRUE)
class_counts

total_class_counts <- plot_data_xgb %>%
  count(compound_class, name = "total_in_dataset")

class_representation <- total_class_counts %>%
  left_join(class_counts, by = "compound_class") %>%
  rename(count_in_dems = n) %>%
  mutate(count_in_dems = replace_na(count_in_dems, 0)) %>%
  mutate(percent_is_dem = (count_in_dems / total_in_dataset) * 100) %>%
  arrange(desc(percent_is_dem))
class_representation

##### summmary stats rf
dems <- plot_data_rf %>%
  filter(p_adj < 0.05 & log2FC > 2)
class_counts <- dems %>%
  count(compound_class, sort = TRUE)
class_counts

total_class_counts <- plot_data_rf %>%
  count(compound_class, name = "total_in_dataset")

class_representation <- total_class_counts %>%
  left_join(class_counts, by = "compound_class") %>%
  rename(count_in_dems = n) %>%
  mutate(count_in_dems = replace_na(count_in_dems, 0)) %>%
  mutate(percent_is_dem = (count_in_dems / total_in_dataset) * 100) %>%
  arrange(desc(percent_is_dem))
class_representation

p_volcano_xgb <- make_volcano(plot_data_xgb)
p_volcano_rf  <- make_volcano(plot_data_rf)

legend <- get_legend(
  p_volcano_rf +
    theme(
      legend.position = "top",
      legend.title = element_text(size = 16, face = "bold"),
      legend.text  = element_text(size = 14),
      legend.key.size = unit(0.8, "cm")
    )
)
p_combined <- plot_grid(
  legend,
  plot_grid(
    p_volcano_xgb + theme(legend.position = "none"),
    p_volcano_rf + theme(legend.position = "none"),
    labels = c("A", "B"),
    label_size = 18,
    ncol = 2,
    align = "hv"
  ),
  ncol = 1,
  rel_heights = c(0.15, 1)
)
ggsave(
  here("misc", "figs/pqn", "volcano_ml.jpg"),
  p_combined,
  width = 16,
  height = 9,
  dpi = 300
)

#################################################################################

# combine 5A-E into single figure 5

levels(met_plot_df$display_class) <- str_wrap(levels(met_plot_df$display_class), width = 20)
names(final_palette) <- str_wrap(names(final_palette), width = 20)

legend_dummy <- ggplot(met_plot_df, aes(x = XGBoost_Importance, y = RandomForest_Importance)) +
  geom_point(aes(fill = display_class, shape = refined_origin)) +
  scale_fill_manual(
    values = final_palette, 
    name = "Compound Class", 
    drop = FALSE
  ) +
  scale_shape_manual(
    values = origin_shapes, 
    name = "Metabolite Origin"
  ) +
  theme_pubr() +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.text = element_text(size = 16, lineheight = 0.8), # lineheight handles wrapped text
    legend.title = element_text(size = 16, face = "bold"),
    legend.key.height = unit(0.8, "cm"), # Increase height to fit wrapped text
    legend.spacing.x = unit(0.2, "cm")
  ) +
  guides(
    fill = guide_legend(
      ncol = 6,           # Set columns to 6; with 18 items, it MUST show 3 rows
      byrow = TRUE,
      order = 1,
      override.aes = list(shape = 21, size = 5, stroke = 0.5)
    ),
    shape = guide_legend(
      nrow = 1,
      order = 2
    )
  )
unified_legend <- get_legend(legend_dummy)


####### figure S4
row_ab_ubqabundance <- plot_grid(pa + theme(legend.position = "none"),
                                 pb + theme(legend.position = "none"), 
                                 labels = c("A", "B"),
                                 label_size = 18)
row_ab_ubqabundance <- plot_grid(unified_legend,row_ab_ubqabundance,
                                 ncol=1,
                                 rel_heights=c(0.3,1))
ggsave(
  here("misc", "figs/pqn", "figs4_ubqabundance.jpg"),
  row_ab_ubqabundance,
  width = 16,
  height = 10,
  dpi = 300
)


######### figure 5
row_ab <- plot_grid(
  p1 + theme(legend.position = "none") + labs(subtitle = "XGBoost"), 
  p2 + theme(legend.position = "none") + labs(subtitle = "Random Forest"), 
  labels = c("A", "B"),
  label_size = 18, ncol = 2
)

row_c <- plot_grid(
  p3 + theme(legend.position = "none"), 
  labels = c("C"),
  label_size = 18, ncol = 1
)

row_de <- plot_grid(
  p_xgb_ubiq_scatter + theme(legend.position = "none") + labs(subtitle = "XGBoost"), 
  p_rf_ubiq_scatter 
  + theme(legend.position = "none")
  + labs(subtitle = "Random Forest"), 
  labels = c("D", "E"),
  label_size = 18, ncol = 2
)

true_final_figure <- plot_grid(
  unified_legend,
  row_ab,
  row_c,
  row_de,
  ncol = 1,
  rel_heights = c(0.6, 1, 1, 1.2, 1) 
)

ggsave(here("misc", "figs/pqn", "fig5_new.pdf"), 
       true_final_figure, width = 15, height = 15, dpi = 600, bg = "white")

