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
library(ggpubr)
library(ggrepel)
library(forcats)
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


#################################################################################
### Volcano plot

cols_origin <- c("Host" = "#97B9CBFF", "Symbiont" = "#9057C6FF", 
                 "Both" = "#FFE1BDFF", "Unknown" = "#8DC657FF")

met_df$refined_origin <- factor(met_df$refined_origin, 
                                levels = c("Host", "Symbiont", "Both", "Unknown"))

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
    neg_log_p_adj = -log10(p_adj)
  )

########################################
# just refined origin color

# Figure S3A

plot_data_volcano <- volcano_results %>%
  inner_join(met_df %>% select(metabolite, refined_origin), by = "metabolite")

m <- nrow(volcano_results)   
sig_threshold <- -log10(0.05 / m)

p_volcano <- ggplot(plot_data_volcano, aes(x = log2FC, y = neg_log_p_adj, color = refined_origin)) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey70") +
  geom_hline(yintercept = sig_threshold, linetype = "dashed", color = "grey70", linewidth = 0.8) +
  geom_point(alpha = 0.6, size = 2.5) +
  facet_wrap(~refined_origin, ncol = 2) +
  xlim(-25,25) +
  ylim(0,100) + 
  scale_color_manual(values = cols_origin) +
  labs(
    x = "log2 Fold Change",
    y = "-log10(adj. p-value)",
    color = "Metabolite Origin",
  ) +
  theme_pubr() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    strip.text = element_text(size = 16),
    # Increase Axis Title Text (labels)
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14)
  )
print(p_volcano)

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

###############################################################

sig_threshold <- -log10(0.05)

# build plot for Figure S3B
p_volcano2 <- ggplot(plot_data_volcano, aes(x = log2FC, y = neg_log_p_adj)) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey70") +
  geom_hline(yintercept = sig_threshold, linetype = "dashed", color = "grey70", linewidth = 0.8) +
  geom_point(aes(color = display_class), alpha = 0.75, size = 2.5) +
  scale_color_manual(
    name = "Compound Class",
    values = class_colors,
    breaks = classes,
    na.value = "gray60"
  ) +
  ylim(0,100) + 
  xlim(-25,25) +
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

#############################################

combined_volcano <- plot_grid(p_volcano, p_volcano2, ncol = 1, labels = c("A", "B"), label_size = 24, align = "hv")
ggsave(here("misc", "figs/pqn", "combined_volcano.jpg"), combined_volcano, width=14,height=20,dpi=300)


