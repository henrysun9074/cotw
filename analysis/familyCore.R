library(tidyverse)
library(knitr)
library(readxl)
library(data.table)
library(vegan)
library(scales)
library(cowplot)
library(ggdendro)
library(ggridges)
library(dendextend)
library(RColorBrewer)
library(corrplot)  
library(reshape2)
library(pheatmap)
library(ggpubr)
library(forcats)
library(ggvenn)
library(ggrepel)
library(ggforce)
library(ComplexUpset)
library(grid)
library(patchwork)
library(UpSetR)
library(here)

# read in data
df <- read.csv(here("Cleaned data CSVs", "qc_data_PQN.csv"))
df$X <- NULL
met_df <- read.csv(here("Cleaned data CSVs", "merged_met_plot_df.csv"))

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
cols_bleaching <- c(
  "Bleached" = "#FF847CFF", 
  "Non-Bleached" = "#019875FF", 
  "Not Applicable" = "#D3D3D3")
cols_location <-c("#002594FF", "#E0B2CDFF", "#54C4E3FF", "#F3AA4FFF")
cols_symbiont  <- c("#D84D16FF", "#FFF800FF", "#8FDA04FF")
cols_phylum <- c("#24492EFF", "#015B58FF", "#2C6184FF", "#59629BFF", "#89689DFF", "#BA7999FF", "#E69B99FF")
cols_sclero    <- c("1" = "#DE7862FF", "0" = "#D8AF39FF")

cols_origin <- c("Host" = "#97B9CBFF", "Symbiont" = "#9057C6FF",
                 "Both" = "#FFE1BDFF", "Unknown" = "#8DC657FF")


################################################################################

### for custom compound_class
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
class_order <- c(target_classes, "Other")

origin_shapes <- c("Host" = 16, "Symbiont" = 3, "Both" = 17, "Unknown" = 8)


################################################################################

# richness by host_family - fig 3a

scler_df <- df %>% filter (df$scleractinia == 1)
scler_df$richness <- rowSums(scler_df %>% select(starts_with("x")) > 0, na.rm = TRUE)
scler_df$host_family <- reorder(scler_df$host_family, scler_df$richness, FUN = mean)

family_palette <- c("#E29191FF", "#99DD92FF", "#93D8B9FF", "#94C4D3FF", "#949ACEFF", 
                    "#B394CCFF", "#CC96B1FF", "#CCA499FF", "#DFE592FF", "#FFA560FF", 
                    "#6BFF63FF", "#65FFCCFF", "#65C4FFFF", "#656BFFFF", "#AD65FFFF", 
                    "#FF65F4FF", "#FF6584FF", "#FF6565FF")

family_counts <- scler_df %>%
  group_by(host_family) %>%
  summarise(n = n(), .groups = "drop")

scler_df <- scler_df %>%
  left_join(family_counts, by = "host_family") %>%
  mutate(family_label = paste0(host_family, " (n=", n, ")"))

scler_df$family_label <- reorder(scler_df$family_label, scler_df$richness, FUN = mean)

p_family_richness <- ggbarplot(
  scler_df, 
  x = "family_label", 
  y = "richness",
  fill = "family_label",      # Match fill to the new label
  color = "black",
  add = "mean_sd",
  error.plot = "pointrange",
  orientation = "horizontal",
  palette = family_palette,   # The 18 colors will map to the new labels
  label = FALSE
) +
  labs(
    y = "Metabolite Richness",
    x = "Host Family"
  ) +
  theme_pubr(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_family_richness)

################################################################################

## flower plot modified to show number of unique metabolites in each family
## fig 3b

draw_flower <- function(data, group_var) {
  group_summary <- data %>%
    pivot_longer(cols = starts_with("x"), names_to = "metabolite", values_to = "val") %>%
    group_by(!!sym(group_var), metabolite) %>%
    summarise(present = any(val > 0, na.rm = TRUE), .groups = "drop")
  
  total_groups <- length(unique(data[[group_var]]))
  metabolite_counts <- group_summary %>%
    group_by(metabolite) %>%
    summarise(n_groups = sum(present), .groups = "drop")
  
  core_count <- sum(metabolite_counts$n_groups == total_groups)
  unique_mets <- metabolite_counts %>% filter(n_groups == 1) %>% pull(metabolite)
  
  petal_data <- group_summary %>%
    filter(present == TRUE) %>%
    group_by(!!sym(group_var)) %>%
    summarise(
      total_n = n(),
      unique_n = sum(metabolite %in% unique_mets),
      .groups = "drop"
    ) %>%
    mutate(
      prop_unique = unique_n / total_n,
      label = paste0(!!sym(group_var), 
                     "\nTotal: ", total_n, 
                     "\nUnique: ", sprintf("%.2f%%", 100 * prop_unique))
      )
  
  n_petals <- nrow(petal_data)
  angle <- seq(0, 2 * pi, length.out = n_petals + 1)[1:n_petals]
  
  petal_data$x <- sin(angle) * 1.5
  petal_data$y <- cos(angle) * 1.5
  
  ggplot(petal_data) +
    geom_ellipse(aes(x0 = x, y0 = y, a = 0.8, b = 1.5, angle = -angle), 
                 fill = "lightblue", alpha = 0.3, color = "steelblue") +
    
    annotate("point", x = 0, y = 0, size = 45, color = "gold", 
             fill = "white", shape = 21, stroke = 2) +
    annotate("text", x = 0, y = 0, label = paste0(core_count), 
             fontface = "bold", size = 7) +
    
    geom_text(aes(x = sin(angle) * 5, y = cos(angle) * 5, label = label), 
              size = 4, lineheight = 0.9) +
    
    theme_void() +
    coord_cartesian(xlim = c(-7, 7), ylim = c(-5, 5)) 
}

df_scler <- df %>% filter(host_order == "Scleractinia", !is.na(host_family))
p_flower_family <- draw_flower(df_scler, "host_family")

################################################################################

# define helper function
get_venn_list <- function(data) {
  list(
    Host = data %>% filter(refined_origin %in% c("Host", "Both")) %>% pull(metabolite),
    Symbiont = data %>% filter(refined_origin %in% c("Symbiont", "Both")) %>% pull(metabolite)
  )
}

# test for class enrichment in family-based core metabolome with fisher's exact test
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

###fisher test for core df
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
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    arrange(p_value)
}

core_fisher <- fisher(core_df, met_df) %>% 
  filter(p_value < 0.05, odds_ratio > 1)


#######################
# venn by host_origin - fig 3c

core_filtered <- core_df %>% 
  filter(refined_origin != "Unknown")

list_core <- get_venn_list(core_filtered)

venn_fill <- c(cols_origin["Host"], cols_origin["Symbiont"])

p_core_venn <- ggvenn(
  list_core, 
  fill_color = venn_fill, 
  stroke_size = 0.5, 
  set_name_size = 5,
  text_size = 4,
  show_percentage = TRUE 
) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

print(p_core_venn)

met_df_filtered <- met_df %>%
  filter(refined_origin != "Unknown")

list_all <- get_venn_list(met_df_filtered)

### plot
cols_origin <- c("Host" = "#97B9CBFF", "Symbiont" = "#9057C6FF",
                 "Both" = "#FFE1BDFF", "Unknown" = "#8DC657FF")
venn_fill <- c(cols_origin["Host"], cols_origin["Symbiont"])

p_a <- ggvenn(list_all, fill_color = venn_fill, stroke_size = 0.5, set_name_size = 4) +
  labs() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))


################################################################################

# volcano - fig 3d
# choose only metabolites in core_df
core_abundance_df <- df %>%
  select(sample, all_of(core_metabolite_ids), scleractinia)

comm_long <- core_abundance_df %>%
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

### volcano plotting data
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

sig_threshold <- -log10(0.05)

## summary statistics
dems <- plot_data_volcano %>%
  filter(p_adj < 0.05 & log2FC > 2)
class_counts <- dems %>%
  count(compound_class, sort = TRUE)
class_counts

total_class_counts <- plot_data_volcano %>%
  count(compound_class, name = "total_in_dataset")

class_representation <- total_class_counts %>%
  left_join(class_counts, by = "compound_class") %>%
  rename(count_in_dems = n) %>%
  mutate(count_in_dems = replace_na(count_in_dems, 0)) %>%
  mutate(percent_is_dem = (count_in_dems / total_in_dataset) * 100) %>%
  arrange(desc(percent_is_dem))
class_representation

p_volcano2 <- ggplot(plot_data_volcano, aes(x = log2FC, y = neg_log_p_adj)) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey70") +
  geom_hline(yintercept = sig_threshold, linetype = "dashed", color = "grey70", linewidth = 0.8) +
  geom_point(aes(color = display_class), alpha = 0.75, size = 3.5) +
  scale_color_manual(
    name = "Compound Class",
    values = class_colors,
    breaks = classes,
    na.value = "gray60"
  ) +
  ylim(0,75) + 
  xlim(-20,20) +
  scale_shape_manual(
    name = "Metabolite Origin",
    values = origin_shapes,
    na.value = 16) +
  facet_wrap(~display_class, ncol = 4) +
  
  guides(
    color = guide_legend(ncol = 2, byrow = TRUE),
    shape = guide_legend(ncol = 1)
  ) +
  labs(
    x = "log2 Fold Change",
    y = "-log10(adj. p-value)",
  ) +
  theme_pubr() +
  theme(
    legend.position = "none",
    strip.text.x = element_text(size = 14),
    axis.title = element_text(size = 20),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

############ combine to make multipanel figure
row1 <- p_family_richness
row2 <- plot_grid(
  p_flower_family, 
  p_core_venn, 
  ncol = 2, 
  rel_widths = c(1.5, 0.7),
  labels = c("B", "C"), 
  label_size = 30,
  label_fontface = "bold"
)
row3 <- p_volcano2
final_multipanel <- plot_grid(
  row1, 
  row2, 
  row3, 
  ncol = 1,
  rel_heights = c(0.8, 0.8, 0.8),
  labels = c("A", "", "D"),  
  label_size = 30,
  label_fontface = "bold"
)

save_plot(
  filename = here("misc", "figs/pqn", "fig3.pdf"),
  plot = final_multipanel,
  base_width = 16,
  base_height = 20,
  device = "pdf",
  bg = "white"
)

################################################################################

### filters for subsequent analysis

met_df_total_all <- met_df %>%
  filter(total_ubiquity == 100)

met_df_scler_total <- met_df %>%
  filter(scler_ubiquity != 0)

met_df_nonscler_total <- met_df %>%
  filter(non_scler_ubiquity != 0)

met_df_scler_all <- met_df %>%
  filter(scler_ubiquity >= 95)

met_df_non_scler_all <- met_df %>%
  filter(non_scler_ubiquity >= 95)

# general function for volcano plots
make_volcano <- function(selected_df, outpath = NULL) {
  
  selected_metabolites <- intersect(selected_df$metabolite, colnames(df))
  
  abundance_df <- df %>%
    select(sample, all_of(selected_metabolites), scleractinia)
  
  comm_long <- abundance_df %>%
    pivot_longer(
      cols = starts_with("x"),
      names_to = "metabolite",
      values_to = "abundance"
    ) %>%
    mutate(abundance = as.numeric(as.character(abundance)))
  
  stats_data <- comm_long %>%
    select(-scleractinia) %>%
    left_join(df %>% select(sample, scleractinia), by = "sample") %>%
    filter(!is.na(scleractinia)) %>%
    mutate(group = if_else(as.character(scleractinia) == "1", "Scleractinia", "Other"))
  
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
      neg_log_p_adj = -log10(p_adj)
    )
  
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
  sig_threshold <- -log10(0.05)
  
  p <- ggplot(plot_data_volcano, aes(x = log2FC, y = neg_log_p_adj)) +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey70") +
    geom_hline(yintercept = sig_threshold, linetype = "dashed", color = "grey70", linewidth = 0.8) +
    geom_point(aes(color = display_class, shape = refined_origin), alpha = 0.75, size = 3.5) +
    scale_color_manual(
      name = "Compound Class",
      values = class_colors,
      breaks = classes,
      na.value = "gray60"
    ) +
    scale_shape_manual(
      name = "Metabolite Origin",
      values = origin_shapes,
      na.value = 16
    ) +
    ylim(0, 75) +
    xlim(-10, 10) +
    labs(
      x = "log2 Fold Change",
      y = "-log10(adj. p-value)"
    ) +
    theme_pubr() +
    theme(
      legend.position = "right",
      axis.title = element_text(size = 20),
      axis.text = element_text(size = 14)
    )
  
  if (!is.null(outpath)) {
    ggsave(outpath, p, width = 14, height = 10, dpi = 300)
  }
  
  return(p)
}

###### Figure s13A
p_volcano_scler_all <- make_volcano(
  met_df_scler_all,
  here("misc", "figs/pqn", "volcano_scler_all.jpg")
)

######### Figure s13B
p_volcano_non_scler_all <- make_volcano(
  met_df_non_scler_all,
  here("misc", "figs/pqn", "volcano_non_scler_all.jpg")
)

# combine figure s13
volcano_ubiquitous <- plot_grid(
  p_volcano_scler_all, 
  p_volcano_non_scler_all, 
  nrow = 2, 
  labels = c("A", "B"), 
  label_size = 18,
  rel_heights = c(1, 1)
)
ggsave(here("misc", "figs/pqn", "volcano_ubq.pdf"), volcano_ubiquitous, width=12,height=14,dpi=300)

################################################################################

## venn diagrams for each combination of metadata variables 

draw_venn_comparison <- function(data, group_var, custom_palette) {
  venn_list <- data %>%
    mutate(!!sym(group_var) := as.character(!!sym(group_var))) %>%
    pivot_longer(cols = starts_with("x"), names_to = "metabolite", values_to = "val") %>%
    filter(val > 0) %>%
    group_by(!!sym(group_var)) %>%
    summarise(mets = list(unique(metabolite)), .groups = "drop") %>%
    filter(lengths(mets) > 0) %>% 
    deframe()
  plot_colors <- if(!is.null(names(custom_palette))) {
    custom_palette[names(venn_list)]
  } else {
    custom_palette[1:length(venn_list)]
  }
  
  ggvenn(
    venn_list, 
    fill_color = plot_colors,
    stroke_size = 0.5, 
    set_name_size = 5,
    text_size = 4
  ) 
}

df_bleach <- df %>% filter(bleaching != "Not Applicable", !is.na(bleaching)) %>%
  filter(scleractinia == 1)

# figure S12
# scleractinian metabolites by bleaching status
venn_bleach <- draw_venn_comparison(data = df_bleach, 
                                    group_var = "bleaching", 
                                    custom_palette = cols_bleaching)

ggsave(
  here("misc", "figs/pqn", "vennbleach.jpg"),
  venn_bleach,
  width = 12,
  height = 8,
  dpi = 300
)

# figure s11
# scleractinian metabolites by location
df_loc <- df %>% filter(!is.na(location)) %>% filter(scleractinia == 1)
venn_loc <- draw_venn_comparison(df_loc, "location", cols_location)

ggsave(
  here("misc", "figs/pqn", "vennloc.jpg"),
  venn_loc,
  width = 12,
  height = 8,
  dpi = 300
)

