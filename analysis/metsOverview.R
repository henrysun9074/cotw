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
library(gt)
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
    # host_order = fct_relevel(factor(host_order), "Scleractinia"),
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
cols_origin <- c("Host" = "#97B9CBFF", "Symbiont" = "#9057C6FF", 
                 "Both" = "#FFE1BDFF", "Unknown" = "#8DC657FF")

##################################################

# summary stats
sum(!is.na(df$host_order)) #542
sum(!is.na(df$host_family)) #542
sum(!is.na(df$host_species)) #479

sum(met_df$refined_origin == "Unknown") ##7789
sum(met_df$refined_origin != "Unknown") ##8579

sum(met_df$refined_origin == "Host") ##3623
sum(met_df$refined_origin == "Both") ##1977
sum(met_df$refined_origin == "Symbiont") ##2979

classification_table <- df %>%
  group_by(host_phylum, host_class, host_order, host_family, host_genus, host_species) %>%
  tally(name = "sample_count") %>%
  ungroup()

## produces a nested classification table
nested_gt_table <- classification_table %>%
  gt(groupname_col = "host_phylum") %>% # Groups by Phylum first
  tab_header(
    title = "Taxonomic Classification of Study Species",
    subtitle = "Hierarchical breakdown from Phylum to Species"
  ) %>%
  cols_label(
    host_class = "Class",
    host_order = "Order",
    host_family = "Family",
    host_genus = "Genus",
    host_species = "Species",
    sample_count = "N"
  ) %>%
  tab_options(
    row_group.font.weight = "bold",
    column_labels.font.weight = "bold"
  )
nested_gt_table

scler_df <- df %>% filter(df$scleractinia == 1, na.rm = TRUE)
out_df <- df %>% filter(df$scleractinia == 0, na.rm = TRUE)
  
taxa_summary <- out_df %>%
  summarise(
    Phyla = n_distinct(host_phylum, na.rm = TRUE),
    Classes = n_distinct(host_class, na.rm = TRUE),
    Orders = n_distinct(host_order, na.rm = TRUE),
    Families = n_distinct(host_family, na.rm = TRUE),
    Genera = n_distinct(host_genus, na.rm = TRUE),
    Species = n_distinct(host_species, na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = everything(), 
    names_to = "Taxonomic Level", 
    values_to = "Unique Count"
  )
print(taxa_summary)

################################################################################

# Bar plots and table of compounds by refined origin for Figure S1

class_counts <- met_plot_df %>%
  filter(compound_class != "Unknown") %>%
  count(compound_class, sort = TRUE) %>%
  mutate(compound_class = fct_inorder(compound_class))
kable(class_counts, col.names = c("Compound Class", "Count"))

# plot bars per compound class
plot_data <- met_plot_df %>%
  filter(compound_class != "Unknown") %>%
  group_by(display_class) %>%
  mutate(class_count = n()) %>%
  ungroup() %>%
  mutate(axis_label = paste0(display_class, " (n = ", class_count, ")"))
label_levels <- plot_data %>%
  select(display_class, axis_label) %>%
  distinct() %>%
  arrange(display_class) %>%
  pull(axis_label)

plot_data$axis_label <- factor(plot_data$axis_label, levels = label_levels)
plot_data <- plot_data %>%
  mutate(refined_origin = factor(refined_origin, 
                                 levels = c("Host", "Symbiont", "Both", "Unknown")))

ordered_classes <- levels(plot_data$display_class)
ordered_classes <- ordered_classes[ordered_classes %in% plot_data$display_class]
axis_colors <- final_palette[ordered_classes]

### Figure S1
compound_bar_plot <- ggplot(plot_data, aes(x = axis_label, fill = refined_origin)) +
  geom_bar(position = "stack", width = 0.7) +
  scale_fill_manual(values = cols_origin) +
  coord_flip() +
  labs(
    x = "Compound Class",
    y = "Number of Metabolites",
    fill = "Metabolite Origin"
  ) +
  theme_pubr() +
  theme(
    axis.text.y = element_text(color = axis_colors, face = "bold"),
    panel.grid.minor = element_blank(),
    axis.label.x = element_text(size = 14),
    legend.position = "bottom"
  )

print(compound_bar_plot)
ggsave(here("misc", "figs/pqn", "compound_barplot.jpg"),
       compound_bar_plot,
       width = 10, height = 12, dpi = 300)

#########################
# host only metabolites
class_counts_host <- met_plot_df %>%
  filter(compound_class != "Unknown") %>%
  filter(refined_origin == "Host") %>%
  count(compound_class, sort = TRUE) %>%
  mutate(compound_class = fct_inorder(compound_class))
kable(class_counts_host, col.names = c("Compound Class", "Count"))

p_class <- ggbarplot(
  class_counts_host,
  x = "compound_class",
  y = "n",
  fill = "compound_class",
  palette = final_palette,
  sort.val = "desc",
  sort.by.groups = FALSE
) +
  geom_text(
    aes(label = paste0("n=", n)),
    hjust = -0.1,   # pushes text slightly outside bars
    size = 4,
    fontface = "bold"
  ) + 
  labs(
    x = "Compound Class",
    y = "# of Metabolites"
  ) +
  theme_pubr() +
  theme(legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16)) +
  ylim(0,300) +
  coord_flip()

p_class
ggsave(here("misc", "figs/pqn", "Fig1hostmetabolites.jpg"), p_class, width = 10, height = 14, dpi = 300)

#########################
# symbiont only metabolites
class_counts_sym <- met_plot_df %>%
  filter(compound_class != "Unknown") %>%
  filter(refined_origin == "Symbiont") %>%
  count(compound_class, sort = TRUE) %>%
  mutate(compound_class = fct_inorder(compound_class))
kable(class_counts_sym, col.names = c("Compound Class", "Count"))

p_class_sym <- ggbarplot(
  class_counts_sym,
  x = "compound_class",
  y = "n",
  fill = "compound_class",
  palette = final_palette,
  sort.val = "desc",
  sort.by.groups = FALSE
) +
  geom_text(
    aes(label = paste0("n=", n)),
    hjust = -0.1,   # pushes text slightly outside bars
    size = 4,
    fontface = "bold"
  ) + 
  labs(
    x = "Compound Class",
    y = "# of Metabolites"
  ) +
  theme_pubr() +
  theme(legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16)) +
  ylim(0,300) +
  coord_flip()
p_class_sym

#########################
# symbiont only metabolites
class_counts_both <- met_plot_df %>%
  filter(compound_class != "Unknown") %>%
  filter(refined_origin == "Both") %>%
  count(compound_class, sort = TRUE) %>%
  mutate(compound_class = fct_inorder(compound_class))
kable(class_counts_both, col.names = c("Compound Class", "Count"))

p_class_both <- ggbarplot(
  class_counts_both,
  x = "compound_class",
  y = "n",
  fill = "compound_class",
  palette = final_palette,
  sort.val = "desc",
  sort.by.groups = FALSE
) +
  geom_text(
    aes(label = paste0("n=", n)),
    hjust = -0.1,   # pushes text slightly outside bars
    size = 4,
    fontface = "bold"
  ) + 
  labs(
    x = "Compound Class",
    y = "# of Metabolites"
  ) +
  theme_pubr() +
  theme(legend.position = "none",
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16)) +
  ylim(0,350) +
  coord_flip()
p_class_both
ggsave(here("misc", "figs/pqn", "Fig1hostmetabolites.jpg"), p_class_both, width = 10, height = 14, dpi = 300)

################################################################################

### get total annotated across available spectral libraries
sum(!is.na(met_df$gnps_compound_name) | 
      !is.na(met_df$compound_name) | 
      !is.na(met_df$coraldb_compound_name))
# 1400

