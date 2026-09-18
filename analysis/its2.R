library(tidyverse)
library(ggplot2)
library(forcats)
library(tidyr)
library(dplyr)
library(stringr)
library(knitr)
library(readxl)
library(data.table)
library(vegan)
library(scales)
library(ComplexUpset)
library(UpSetR)
library(rstatix)
library(ggpubr)
library(patchwork)
library(cowplot)
library(here)

# read data
df <- read.csv(here("Cleaned data CSVs", "ITS2full_PQN_Sep17.csv"))
df$X <- NULL
met_df <- read.csv(here("Cleaned data CSVs", "merged_met_plot_df.csv"))

# define color palettes
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

origin_shapes <- c("Host" = 16, "Symbiont" = 3, "Both" = 17, "Unknown" = 8)

######### define factor levels for its2

# recode ITS2.letter
# A = Symbiodinaceae
# B = Breviolum
# C = Cladocopium
# D = Durusdinium
# Keep mix and no seq the same

df <- df %>%
  mutate(
    ITS2.Letter = case_when(
      is.na(ITS2.Letter)           ~ NA_character_,
      str_detect(ITS2.Letter, ":") ~ "Mix",
      ITS2.Letter == "A" ~ "Symbiodinium",
      ITS2.Letter == "B" ~ "Breviolum",
      ITS2.Letter == "C" ~ "Cladocopium",
      ITS2.Letter == "D" ~ "Durusdinium",
      TRUE ~ ITS2.Letter   # keeps "Mix", "NoSeq", etc.
    )
  )

its2_palette <- c(
  "Breviolum"        = "#FF0000FF",
  "Cladocopium"      = "#00A08AFF",
  "Durusdinium"      = "#F2AD00FF",
  "Symbiodinium"  = "#5BBCD6FF",
  "Mix"              = "#F98400FF"
)

its2_levels <- c(
  "Symbiodinium",
  "Breviolum",
  "Cladocopium",
  "Durusdinium",
  "Mix"
)


### check no seq samples
na_breakdown <- df %>%
  filter(ITS2.Letter == "No Seq") %>%
  count(location, bleaching, name = "count") %>%
  mutate(percentage = round(count / sum(count) * 100, 1))
print(na_breakdown)

################################################################################

# Figure S7

its2_barplot <- df %>%
  filter(!is.na(ITS2.Letter), ITS2.Letter != "") %>%
  count(ITS2.Letter) %>%
  mutate(
    ITS2.Letter = factor(ITS2.Letter, 
                         levels = c("No Seq","Mix","Durusdinium","Cladocopium","Breviolum","Symbiodinium"))
  ) 

bar1 <- ggplot(its2_barplot, aes(x = ITS2.Letter, y = n, fill = ITS2.Letter)) +
  geom_col() +
  geom_text(aes(label = paste0("n=", n)), hjust = -0.2) +
  scale_fill_manual(values = its2_palette, breaks = its2_levels) +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = "Symbiont Genus",
  ) +
  theme_pubr()+
  ylim(0,150)+
  theme(legend.position = "none",
        axis.title.x = element_blank())
bar1
ggsave(
  here("misc", "figs/pqn", "bar1_figs7.pdf"), 
  plot = bar1, 
  width = 7, 
  height = 8, 
  dpi = 300
)

################################################################################

# Figure 4A

bar2_df <- df %>%
  filter(
    !is.na(host_family), host_family != "",
    !is.na(ITS2.Letter), ITS2.Letter != "", 
    ITS2.Letter != "No Seq"
  ) %>%
  count(host_family = host_family, ITS2.Letter, name = "n") %>%
  group_by(host_family) %>%
  mutate(
    total_n = sum(n),
    prop = n / total_n
  ) %>%
  ungroup() %>%
  mutate(
    host_family_label = paste0(host_family, " (n=", total_n, ")"),
    host_family_label = fct_reorder(host_family_label, total_n, .desc = TRUE)
  )

bar2 <- ggplot(bar2_df, aes(x = host_family_label, y = prop, fill = ITS2.Letter)) +
  geom_col() +
  scale_fill_manual(values = its2_palette, breaks = its2_levels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "Scleractinian Family",
    y = "Proportion of Samples",
    fill = "Symbiont Genus"
  ) + coord_flip() +
  theme_pubr() +
  theme(axis.text.y = element_text(size = 10))
bar2

################################################################################

# Figure S9A

bar3_df <- df %>%
  filter(
    !is.na(ITS2.Letter), ITS2.Letter != "",
    ITS2.Letter != "No Seq",
    !is.na(bleaching), bleaching != "Not Applicable"
  ) %>%
  count(bleaching, ITS2.Letter, name = "n") %>%
  group_by(bleaching) %>%
  mutate(
    total_n = sum(n),
    prop = n / total_n
  ) %>%
  ungroup() %>%
  mutate(
    bleaching_label = paste0(bleaching, " (n=", total_n, ")")
  )

bar3 <- ggplot(bar3_df, aes(x = bleaching_label, y = prop, fill = ITS2.Letter, alpha = bleaching)) +
  geom_col() +
  scale_fill_manual(values = its2_palette, breaks = its2_levels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_alpha_manual(
    values = c("Bleached" = 1, "Non-Bleached" = 1)  
  ) +
  labs(
    y = "Proportion of Samples",
    fill = "Symbiont Genus",
  ) + guides(alpha = "none") +
  theme_pubr() +
  theme(axis.title.x = element_blank())
bar3

################################################################################

# Figure 4B
# by location

bar4_df <- df %>%
  filter(
    !is.na(ITS2.Letter), ITS2.Letter != "",
    ITS2.Letter != "No Seq",
    !is.na(location)
  ) %>%
  count(location, ITS2.Letter, name = "n") %>%
  group_by(location) %>%
  mutate(
    total_n = sum(n),
    prop = n / total_n
  ) %>%
  ungroup() %>%
  mutate(
    loc_label = paste0(location, " (n=", total_n, ")")
  )

bar4 <- ggplot(bar4_df, aes(x = loc_label, y = prop, fill = ITS2.Letter)) +
  geom_col() +
  scale_fill_manual(values = its2_palette, breaks = its2_levels) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    y = "Proportion of Samples",
    fill = "Symbiont Genus",
  ) + guides(alpha = "none") +
  theme_pubr() +
  theme(axis.title.x = element_blank())
bar4

################################################################################
# Figure 4D
# richness 

metabolite_cols <- grep("^x", names(df), value = TRUE)

richness_df <- df %>%
  mutate(ITS2.Letter = factor(ITS2.Letter, levels = its2_levels)) |> 
  filter(
    !is.na(ITS2.Letter),
    ITS2.Letter != "No Seq",
    ITS2.Letter != "Mix"
  ) %>%
  select(sample_id, ITS2.Letter, all_of(metabolite_cols)) %>%
  rowwise() %>%
  mutate(MetabolomicRichness = sum(c_across(all_of(metabolite_cols)) > 0, na.rm = TRUE)) %>%
  ungroup()

kruskal_test_res <- richness_df %>%
  kruskal_test(MetabolomicRichness ~ ITS2.Letter)

stat.test <- richness_df %>%
  dunn_test(MetabolomicRichness ~ ITS2.Letter, p.adjust.method = "BH") %>%
  filter(p.adj < 0.05) 

max_y <- max(richness_df$MetabolomicRichness, na.rm = TRUE)
n_comparisons <- nrow(stat.test)

stat.test <- stat.test %>%
  mutate(y.position = seq(from = max_y * 1.1, 
                          by = max_y * 0.05, 
                          length.out = n_comparisons))

# plot 
richness <- ggplot(richness_df, aes(x = ITS2.Letter, y = MetabolomicRichness, fill = ITS2.Letter)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, aes(color = ITS2.Letter), median.linewidth = 2) + 
  geom_jitter(aes(color = ITS2.Letter), width = 0.2, size = 1.8, alpha = 0.3) +
  scale_fill_manual(values = its2_palette) +
  scale_color_manual(values = its2_palette) +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.01) +
  theme_pubr() + 
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

################################################################################
# Figure 4E
# Shannon entropy

shannon_index <- function(counts) {
  counts <- counts[counts > 0]   # drop zeros
  total <- sum(counts)
  if (total == 0) return(0)      # if no metabolites detected
  p <- counts / total
  -sum(p * log(p))
}

metabolite_cols <- grep("^x", names(df), value = TRUE)

# calculate metabolic entropy per sample
entropy_df <- df %>%
  mutate(ITS2.Letter = factor(ITS2.Letter, levels = its2_levels)) |>
  filter(
    !is.na(ITS2.Letter),
    ITS2.Letter != "No Seq",
    ITS2.Letter != "Mix"
  ) |> 
  select(sample_id, ITS2.Letter, all_of(metabolite_cols)) %>%
  rowwise() %>%
  mutate(MetabolicEntropy = shannon_index(c_across(all_of(metabolite_cols)))) %>%
  ungroup()

entropy_df$ITS2.Letter <- droplevels(entropy_df$ITS2.Letter)
kruskal_entropy <- entropy_df %>%
  kruskal_test(MetabolicEntropy ~ ITS2.Letter)

stat.test_entropy <- entropy_df %>%
  dunn_test(MetabolicEntropy ~ ITS2.Letter, p.adjust.method = "BH") %>%
  arrange(p.adj)

stat.test_entropy <- stat.test_entropy %>%
  filter(p.adj < 0.05)

max_y_ent <- max(entropy_df$MetabolicEntropy, na.rm = TRUE)
n_comp_ent <- nrow(stat.test_entropy)

stat.test_entropy <- stat.test_entropy %>%
  mutate(y.position = seq(from = max_y_ent * 1.05, 
                          by = max_y_ent * 0.07, 
                          length.out = n_comp_ent))

# plot
entropy <- ggplot(entropy_df, aes(x = ITS2.Letter, y = MetabolicEntropy, fill = ITS2.Letter)) +
  geom_boxplot(
    aes(color = ITS2.Letter),
    alpha = 0.8, 
    outlier.shape = NA, 
    median.linewidth = 2
  ) + 
  geom_jitter(
    aes(color = ITS2.Letter), 
    width = 0.2, 
    size = 1.8, 
    alpha = 0.3
  ) +
  scale_fill_manual(values = its2_palette) +
  scale_color_manual(values = its2_palette) +
  stat_pvalue_manual(
    stat.test_entropy, 
    label = "p.adj.signif", 
    tip.length = 0.01
  ) +
  labs(
    y = "Shannon Entropy"
  ) +
  ylim(5, 8) +
  theme_pubr() + theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank()
  )

################################################################################

## pcoas Figure 4C, S8, S9B
## filter to remove "mix" samples
keep_genera <- c("Symbiodinium", "Breviolum", "Cladocopium", "Durusdinium")

# Filter and ensure numeric data is clean
corals_4g <- df %>%
  filter(ITS2.Letter %in% keep_genera, !is.na(sample_id))

permanova_numeric_data <- corals_4g %>% select(starts_with("x"))
keep_rows <- rowSums(is.na(permanova_numeric_data)) < ncol(permanova_numeric_data)

# Final clean subsets
permanova_numeric_data2 <- permanova_numeric_data[keep_rows, , drop = FALSE]
meta2 <- corals_4g[keep_rows, , drop = FALSE] %>%
  mutate(
    ITS2.Letter = factor(ITS2.Letter, levels = its2_levels),
    location = factor(location, levels = c("Sri Lanka", "Curaçao", "Hawaii"))
  )
permanova_numeric_data2[is.na(permanova_numeric_data2)] <- 0

numeric_only_data <- meta2 %>% 
  select(starts_with("x")) %>% 
  select(where(is.numeric)) 

permanova_matrix <- as.matrix(numeric_only_data)
bray_curtis_4g <- vegdist(permanova_matrix, method = "bray")

################################################################################

# Figure S8
# nested permanova
permanova_result <- adonis2(
  bray_curtis_4g ~ ITS2.Letter / location,
  data = meta2,
  permutations = 999,
  by = "margin"
)
print(permanova_result)

pcoa_result <- cmdscale(bray_curtis_4g, eig = TRUE, k = 2)
pcoa_points <- as.data.frame(pcoa_result$points)
colnames(pcoa_points) <- c("PCoA1", "PCoA2")
pcoa_points <- bind_cols(pcoa_points, meta2)

pos_eig <- pcoa_result$eig[pcoa_result$eig > 0]
var_explained <- round(100 * pcoa_result$eig[1:2] / sum(pos_eig), 1)

perm_p_val <- min(permanova_result$`Pr(>F)`, na.rm = TRUE)
perm_r2    <- round(sum(permanova_result$R2[1]), 3)
p_label    <- if(perm_p_val <= 0.001) "p < 0.001" else paste("p =", perm_p_val)
stats_annotation_combined <- paste0("PERMANOVA: R² = ", perm_r2, ", ", p_label)

p <- ggplot(pcoa_points, aes(x = PCoA1, y = PCoA2, color = ITS2.Letter, fill = ITS2.Letter)) +
  geom_point(size = 2, alpha = 0.8, aes(color = ITS2.Letter, shape = location)) +
  stat_ellipse(
    geom = "polygon",
    alpha = 0.15,
    level = 0.95,
    type = "t",
    colour = NA
  ) +
  annotate(
    "text",
    x = -Inf,  
    y = -Inf, 
    label = stats_annotation_combined,
    hjust = -0.05,
    vjust = -0.6,
    size = 4) +
  scale_color_manual(values = its2_palette) +
  scale_fill_manual(values = its2_palette) +
  scale_shape_manual(values = c("Sri Lanka" = 15, "Curaçao" = 11, "Hawaii" = 16)) + 
  labs(
    x = paste0("PCoA1: (", var_explained[1], "%)"),
    y = paste0("PCoA2: (", var_explained[2], "%)"),
    color = "Symbiont Genus",
    fill = "Symbiont Genus",
    shape = "Location"
  ) +
  theme_pubr() +
  theme(legend.position = "right")
p
ggsave(here("misc", "figs/pqn", "its2_permanova_location.pdf"), 
       p, width = 7, height = 6, dpi = 300)
################################################################################

# Figure 4C

its2_permanova_result <- adonis2(
  bray_curtis_4g ~ ITS2.Letter,
  data = meta2,
  permutations = 999,
  by = "margin"
)
print(its2_permanova_result)

pcoa_result <- cmdscale(bray_curtis_4g, eig = TRUE, k = 2)
pcoa_points <- as.data.frame(pcoa_result$points)
colnames(pcoa_points) <- c("PCoA1", "PCoA2")
pcoa_points <- bind_cols(pcoa_points, meta2)

pos_eig <- pcoa_result$eig[pcoa_result$eig > 0]
var_explained <- round(100 * pcoa_result$eig[1:2] / sum(pos_eig), 1)

perm_p_val <- min(its2_permanova_result$`Pr(>F)`, na.rm = TRUE)
perm_r2    <- round(sum(its2_permanova_result$R2[1]), 3)
p_label    <- if(perm_p_val <= 0.001) "p < 0.001" else paste("p =", perm_p_val)
stats_annotation_combined <- paste0("PERMANOVA: R² = ", perm_r2, ", ", p_label)

p2 <- ggplot(pcoa_points, aes(x = PCoA1, y = PCoA2, color = ITS2.Letter, fill = ITS2.Letter)) +
  geom_point(size = 2, alpha = 0.8, aes(color = ITS2.Letter)) +
  stat_ellipse(
    geom = "polygon",
    alpha = 0.15,
    level = 0.95,
    type = "t",
    colour = NA
  ) +
  annotate(
    "text",
    x = -Inf,  
    y = -Inf, 
    label = stats_annotation_combined,
    hjust = -0.05,
    vjust = -0.6,
    size = 4) +
  scale_color_manual(values = its2_palette) +
  scale_fill_manual(values = its2_palette) +
  labs(
    x = paste0("PCoA1: (", var_explained[1], "%)"),
    y = paste0("PCoA2: (", var_explained[2], "%)"),
    color = "Symbiont Genus",
    fill = "Symbiont Genus"
  ) +
  theme_pubr() +
  theme(legend.position = "right")
p2

################################################################################

#Figure S9B

bleaching_permanova_result <- adonis2(
  bray_curtis_4g ~ ITS2.Letter / bleaching,
  data = meta2,
  permutations = 999,
  by = "margin"
)
print(bleaching_permanova_result)

pcoa_result <- cmdscale(bray_curtis_4g, eig = TRUE, k = 2)
pcoa_points <- as.data.frame(pcoa_result$points)
colnames(pcoa_points) <- c("PCoA1", "PCoA2")
pcoa_points <- bind_cols(pcoa_points, meta2)

pos_eig <- pcoa_result$eig[pcoa_result$eig > 0]
var_explained <- round(100 * pcoa_result$eig[1:2] / sum(pos_eig), 1)

perm_p_val <- min(bleaching_permanova_result$`Pr(>F)`, na.rm = TRUE)
perm_r2    <- round(sum(bleaching_permanova_result$R2[1]), 3)
p_label    <- if(perm_p_val <= 0.001) "p < 0.001" else paste("p =", perm_p_val)
stats_annotation_combined <- paste0("PERMANOVA: R² = ", perm_r2, ", ", p_label)

p3 <- ggplot(pcoa_points, aes(x = PCoA1, y = PCoA2, color = ITS2.Letter, fill = ITS2.Letter)) +
  geom_point(size = 2, alpha = 0.8, aes(color = ITS2.Letter, shape = bleaching)) +
  stat_ellipse(
    geom = "polygon",
    alpha = 0.15,
    level = 0.95,
    type = "t",
    colour = NA
  ) +
  annotate(
    "text",
    x = -Inf,  
    y = -Inf, 
    label = stats_annotation_combined,
    hjust = -0.05,
    vjust = -0.6,
    size = 4) +
  scale_color_manual(values = its2_palette, breaks = its2_levels) +
  scale_fill_manual(values = its2_palette, breaks = its2_levels) +
  scale_shape_manual(values = c("Bleached" = 12, "Non-Bleached" = 16)) + 
  labs(
    x = paste0("PCoA1: (", var_explained[1], "%)"),
    y = paste0("PCoA2: (", var_explained[2], "%)"),
    color = "Symbiont Genus",
    fill = "Symbiont Genus",
    shape = "Bleaching Status"
  ) +
  theme_pubr() +
  theme(legend.position = "right")
p3


################################### permanova outputs

tidy_permanova <- function(model_output, model_name) {
  model_output %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Term") %>%
    filter(!Term %in% c("Residual", "Total")) %>%
    mutate(Model = model_name) %>%
    select(Model, Term, Df, SumOfSqs, R2, F, `Pr(>F)`)
}

table_its2 <- tidy_permanova(its2_permanova_result, "ITS2 Letter")
table_its2loc   <- tidy_permanova(permanova_result, "ITS2 Letter x Location")
table_its2b <-  tidy_permanova(bleaching_permanova_result, "ITS2 Letter x Bleaching")

# Combine into one master table
full_stats_table <- bind_rows(table_its2, table_its2loc, table_its2b)

full_stats_table <- full_stats_table %>%
  mutate(
    across(c(SumOfSqs, R2, F), ~ round(., 3)),
    `Pr(>F)` = ifelse(`Pr(>F)` <= 0.001, "< 0.001", as.character(`Pr(>F)`))
  )

print(full_stats_table)
write.csv(full_stats_table, here("misc", "PERMANOVA_pqn_results_ITS2.csv"))

################################################################################

# keep Mix
upset_input <- df %>%
  filter(!is.na(ITS2.Letter), ITS2.Letter != "No Seq", ITS2.Letter != "Mix") %>%
  select(ITS2.Letter, all_of(metabolite_cols)) %>%
  pivot_longer(cols = starts_with("x", ignore.case = FALSE), names_to = "metabolite", values_to = "abundance") %>%
  group_by(ITS2.Letter, metabolite) %>%
  summarise(present = as.numeric(any(abundance > 0, na.rm = TRUE)), .groups = "drop") %>%
  pivot_wider(names_from = ITS2.Letter, values_from = present, values_fill = 0) %>%
  as.data.frame()

target_order <- c("Durusdinium", "Cladocopium", "Breviolum", "Symbiodinium")
upset_data <- upset_input[, -1]
UpSetR::upset(
  upset_data,
  sets = target_order,
  keep.order = TRUE,
  main.bar.color = "gray20",
  sets.bar.color = its2_palette[target_order],
  order.by = "freq", 
  decreasing = TRUE,
  point.size = 3.5, 
  line.size = 1.5,
  text.scale = c(1.3, 1.3, 1, 1, 1.5, 1)
)

## ComplexUpset

target_order <- c("Durusdinium", "Cladocopium", "Breviolum", "Symbiodinium")
names(its2_palette) <- its2_levels

# remove global metabolites not present in A/b/c/d/mix
upset_input_with_meta <- df %>%
  filter(ITS2.Letter %in% target_order) %>%
  select(ITS2.Letter, all_of(metabolite_cols)) %>%
  pivot_longer(
    cols = starts_with("x", ignore.case = FALSE), 
    names_to = "metabolite", 
    values_to = "abundance"
  ) %>%
  group_by(metabolite) %>%
  filter(sum(abundance, na.rm = TRUE) > 0) %>% 
  group_by(ITS2.Letter, metabolite) %>%
  summarise(
    present = as.numeric(any(abundance > 0, na.rm = TRUE)), 
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = ITS2.Letter, values_from = present, values_fill = 0) %>%
  left_join(met_df %>% select(metabolite, compound_class), by = "metabolite") %>%
  mutate(
    compound_class = trimws(as.character(compound_class)),
    display_class = if_else(
      compound_class %in% names(final_palette), 
      compound_class, 
      "Other"
    ),
    display_class = factor(display_class, levels = c(target_classes, "Other"))
  ) %>%
  as.data.frame()

upset_plot <- ComplexUpset::upset(
  upset_input_with_meta,
  target_order, 
  name = "",
  width_ratio = 0.05, 
  
  set_sizes = (
    upset_set_size() + 
      theme(
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()
      )
  ),
  
  stripes = upset_stripes(
    mapping = aes(color = group),
    colors = its2_palette
  ),
  
  base_annotations = list(
    'Intersection' = intersection_size(
      counts = TRUE,
      mapping = aes(fill = "gray20"),
      text = list(size = 2.2, vjust = -0.5) 
    ) + 
      scale_fill_identity() +
      theme_pubr() +
      theme(
        panel.grid = element_blank(), 
        axis.text.y = element_text(size = 8),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()
      )
  ),
  
  sort_intersections = 'descending',
  sort_sets = FALSE, 
  themes = upset_default_themes(text = element_text(size = 14))
)
print(upset_plot)

################################################################################

## combine plots for Figure 4

shared_legend <- get_legend(
  bar2 + 
    theme(
      legend.position = "top",
      legend.direction = "horizontal",
      legend.text = element_text(size = 18),        
      legend.title = element_text(size = 18))
)

p_bar2 <- bar2 + theme(legend.position = "none")
p_bar4 <- bar4 + theme(legend.position = "none")
p_pcoa <- p2   + theme(legend.position = "none")
p_rich <- richness + theme(legend.position = "none")
p_ent  <- entropy  + theme(legend.position = "none")

top_row <- plot_grid(
  p_bar2, p_bar4, p_pcoa,
  ncol = 3,
  rel_widths = c(0.6, 0.8, 0.8),
  labels = c("A", "B", "C"),
  label_size = 18,
  label_fontface = "bold"
)

bottom_row <- plot_grid(
  p_rich, p_ent, upset_plot,
  ncol = 3,
  rel_widths = c(0.4, 0.4, 1.5),
  labels = c("D", "E", "F"),
  label_size = 18,
  label_fontface = "bold"
)

main_plots <- plot_grid(
  top_row, 
  bottom_row, 
  ncol = 1, 
  rel_heights = c(0.8, 1)
)

final_figure <- plot_grid(
  shared_legend, 
  main_plots, 
  ncol = 1, 
  rel_heights = c(0.1, 1.9) 
)

print(final_figure)
ggsave(
  here("misc", "figs/pqn", "ITS2combined.pdf"), 
  plot = final_figure, 
  width = 18, 
  height = 12, 
  dpi = 300
)

################################################################################

# Combine Figure S9
bar3 <- bar3 + theme(legend.position = "none")
bleach_its2 <- plot_grid(bar3, p3, labels = c ("A","B"), label_size = 18,ncol = 1)
ggsave(
  here("misc", "figs/pqn", "ITS2bleaching.pdf"), 
  plot = bleach_its2, 
  width = 7, 
  height = 8, 
  dpi = 300
)

################################################################################
# Figure S10
## compound class analysis - volcano plots

target_genera <- c("Symbiodinium", "Breviolum", "Cladocopium", "Durusdinium")
genus_pairs <- combn(target_genera, 2, simplify = FALSE)

pairwise_data <- df %>%
  filter(ITS2.Letter %in% target_genera)

## pairwise combinations A vs B, A vs C, A vs D, B vs C, B vs D, C vs D
compute_pairwise_stats <- function(pair, data, met_metadata, palette) {
  group_a <- pair[1]
  group_b <- pair[2]
  
  comp_df <- data %>% filter(ITS2.Letter %in% pair)
  
  stats <- comp_df %>%
    select(ITS2.Letter, starts_with("x", ignore.case = FALSE)) %>% 
    pivot_longer(
      cols = -ITS2.Letter,        
      names_to = "metabolite", 
      values_to = "abundance"
    ) %>%
    mutate(abundance = as.numeric(abundance)) %>%
    group_by(metabolite) %>%
    summarise(
      mean_a = mean(abundance[ITS2.Letter == group_a], na.rm = TRUE),
      mean_b = mean(abundance[ITS2.Letter == group_b], na.rm = TRUE),
      p_raw  = wilcox.test(abundance ~ ITS2.Letter)$p.value,
      .groups = "drop"
    ) %>%
    mutate(
      p_adj = p.adjust(p_raw, method = "bonferroni"), 
      log2FC = log2((mean_a + 0.01) / (mean_b + 0.01)), 
      neg_log_p = -log10(p_adj),
      comparison = paste(group_a, "vs", group_b)
    ) %>%
    inner_join(met_metadata %>% select(metabolite, compound_class), by = "metabolite") %>%
    mutate(display_class = if_else(compound_class %in% names(palette), compound_class, "Other"))
  
  return(stats)
}

met_df_known <- met_df %>%
  filter(met_df$compound_class != "Unknown")
  
stats_list <- lapply(genus_pairs, function(p) {
  compute_pairwise_stats(p, pairwise_data, met_df_known, final_palette)
})

volcano_list <- lapply(stats_list, function(df) {
  ggplot(df, aes(x = log2FC, y = neg_log_p)) +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey80") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey80") +
    geom_point(aes(color = display_class), alpha = 0.6, size = 2) +
    scale_color_manual(values = final_palette, name = "Compound Class") + 
    labs(title = unique(df$comparison),
         x = "log2 Fold Change", y = "-log10(adj. p)") +
    xlim(-10, 10) +
    ylim(0, 20) +
    theme_pubr() +
    theme(legend.position = "bottom", # Positioned for extraction
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 10))
})

legend_b <- get_legend(
  volcano_list[[1]] + 
    theme(legend.position = "bottom",
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12)) +
    guides(color = guide_legend(nrow = 3, override.aes = list(size = 4))) # Makes legend "taller"
)

volcano_list_no_legend <- lapply(volcano_list, function(p) p + theme(legend.position = "none"))

pairwise_volcano_grid <- plot_grid(
  plotlist = volcano_list_no_legend,
  ncol = 3,
  nrow = 2,
  labels = "AUTO",
  label_size = 20
)

final_pairwise_plot <- plot_grid(
  pairwise_volcano_grid,
  legend_b,
  ncol = 1,
  rel_heights = c(1, 0.15) 
)

ggsave(here("misc", "figs/pqn", "its2_volcano.pdf"), 
       final_pairwise_plot, width = 18, height = 14, dpi = 300)