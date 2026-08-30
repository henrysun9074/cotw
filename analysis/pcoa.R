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
library(forcats)
library(here)

# read in data
df <- read.csv(here("Cleaned data CSVs", "qc_data_PQN.csv"))
df$X <- NULL
met_df <- read.csv(here("Cleaned data CSVs", "merged_met_plot_df.csv"))

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

permanova_numeric_data <- df %>% select(starts_with("x"))
keep_rows <- !is.na(df$host_family) & 
  complete.cases(permanova_numeric_data)

permanova_numeric_data2 <- permanova_numeric_data[keep_rows, , drop = FALSE]
bray_curtis <- vegdist(permanova_numeric_data2, method = "bray")

meta2 <- df[keep_rows, , drop = FALSE]

# ################################################################################
# #Scleractinia vs nonScleractinia PERMANOVA - Figure 2A

# compute BC dissimilarity
bray_curtis_scleractinia <- vegdist(permanova_numeric_data2, method = "bray")
scleractinia_permanova_result <- adonis2(
  bray_curtis_scleractinia ~ scleractinia,
  data = meta2,
  permutations = 999
)
print(scleractinia_permanova_result)

################ plot
pcoa_result <- cmdscale(bray_curtis_scleractinia, eig = TRUE, k = 2)
pcoa_points <- as.data.frame(pcoa_result$points)
colnames(pcoa_points) <- c("PCoA1", "PCoA2")
pcoa_points <- bind_cols(pcoa_points, meta2)

# compute percent variance explained (use only positive eigenvalues)
pos_eig <- pcoa_result$eig[pcoa_result$eig > 0]
var_explained <- round(100 * pcoa_result$eig[1:2] / sum(pos_eig), 1)

perm_p_val <- scleractinia_permanova_result$`Pr(>F)`[1]
perm_r2    <- round(scleractinia_permanova_result$R2[1], 3)
p_label <- if(perm_p_val == 0.001) "p < 0.001" else paste("p =", perm_p_val)
stats_annotation <- paste0("PERMANOVA: R² = ", perm_r2, ", ", p_label)

shapes_location <- c(16, 17, 15, 18) 
p <- ggplot(pcoa_points, aes(x = PCoA1, y = PCoA2, 
                             color = scleractinia, 
                             fill = scleractinia, 
                             shape = location)) + # Map shape to location
  geom_point(size = 3, alpha = 0.8) +
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
    label = stats_annotation,
    hjust = -0.05,
    vjust = -0.6,
    size = 4) +
  scale_color_manual(values = cols_sclero, labels = c("1" = "Scleractinia", "0" = "Other")) +
  scale_fill_manual(values = cols_sclero, labels = c("1" = "Scleractinia", "0" = "Other")) +
  # 2. Add the manual shape scale
  scale_shape_manual(values = shapes_location) + 
  labs(
    x = paste0("PCoA1: (", var_explained[1], "%)"),
    y = paste0("PCoA2: (", var_explained[2], "%)"),
    color = "Order",
    fill = "Order",
    shape = "Location"
  ) +
  theme_cowplot(font_size = 18) +
  theme(legend.position = "right")

################################################################################
#Location + bleaching status PERMANOVA - Figure S6A

# drop NAs for bleaching
keep_rows <- complete.cases(permanova_numeric_data2) & 
  !is.na(meta2$location) & 
  !is.na(meta2$bleaching)

meta_clean <- meta2[keep_rows, ]
numeric_clean <- permanova_numeric_data2[keep_rows, ]

bray_curtis_clean <- vegdist(numeric_clean, method = "bray")
locB_permanova_result <- adonis2(
  bray_curtis_clean ~ location / bleaching, 
  data = meta_clean, 
  permutations = 999
)
print(locB_permanova_result)          

perm_p_val <- locB_permanova_result$`Pr(>F)`[1]
perm_r2    <- round(locB_permanova_result$R2[1], 3)

p_label <- if (perm_p_val <= 0.001) {
  "p \u2264 0.001"
} else {
  paste0("p = ", format(perm_p_val, digits = 3))
}

permanova_label <- paste0(
  "PERMANOVA: ",
  p_label,
  ", R\u00b2 = ",
  sprintf("%.3f", perm_r2)
)
stats_annotation <- paste0("PERMANOVA: R² = ", perm_r2, ", ", p_label)

shapes_bleaching <- c(
  "Bleached" = 4,       # Cross
  "Non-Bleached" = 17,   # Triangle
  "Not Applicable" = 1   # Open Circle
)

pcoa_result2 <- cmdscale(bray_curtis_clean, eig = TRUE, k = 2)
pcoa_points2 <- as.data.frame(pcoa_result2$points)
colnames(pcoa_points2) <- c("PCoA1", "PCoA2")
pcoa_points2 <- bind_cols(pcoa_points2, meta_clean)

pos_eig <- pcoa_result2$eig[pcoa_result2$eig > 0]
var_explained <- round(100 * pcoa_result2$eig[1:2] / sum(pos_eig), 1)

p2 <- ggplot(pcoa_points2, aes(x = PCoA1, y = PCoA2, color = location)) +
  geom_point(aes(shape = bleaching), size = 3, alpha = 0.8) +
  scale_color_manual(values = cols_location) +
  scale_fill_manual(values = cols_location) +
  scale_shape_manual(values = shapes_bleaching) + 
  labs(
    x = paste0("PCoA1: (", var_explained[1], "%)"),
    y = paste0("PCoA2: (", var_explained[2], "%)"),
    color = "Location",
    fill = "Location",
    shape = "Bleaching Status"
  ) +
  annotate(
    "text",
    x = -Inf,
    y = -Inf,
    label = stats_annotation,
    hjust = -0.05,
    vjust = -0.6,
    size = 4) +
    theme(legend.position = "right") +
  theme_cowplot(font_size = 14)
# ggsave(here("misc", "figs/pqn", "pcoa_locb.jpg"), p2, width = 8, height = 6, dpi = 300)


################################################################################
## Location by symbiont potential PERMANOVA - Figure S6B

keep_rows <- !is.na(meta2$symbiont.potential)
meta_symb <- meta2[keep_rows, , drop = FALSE]
numeric_symb <- permanova_numeric_data2[keep_rows, , drop = FALSE]

bray_curtis_locS <- vegdist(numeric_symb, method = "bray")
locS_permanova_result <- adonis2(
  bray_curtis_locS ~ location / symbiont.potential, 
  data = meta_symb, 
  permutations = 999
)
print(locS_permanova_result)          

# global model stats
perm_p_val <- max(locS_permanova_result$`Pr(>F)`[1])
perm_r2    <- round(sum(locS_permanova_result$R2[1]), 3)
p_label    <- if(perm_p_val <= 0.001) "p < 0.001" else paste("p =", perm_p_val)
stats_annotation <- paste0("PERMANOVA: R² = ", perm_r2, ", ", p_label)

pcoa_result3 <- cmdscale(bray_curtis_locS, eig = TRUE, k = 2)
pcoa_points3 <- as.data.frame(pcoa_result3$points)
colnames(pcoa_points3) <- c("PCoA1", "PCoA2")
pcoa_points3 <- bind_cols(pcoa_points3, meta_symb)

pos_eig <- pcoa_result3$eig[pcoa_result3$eig > 0]
var_explained <- round(100 * pcoa_result3$eig[1:2] / sum(pos_eig), 1)

shapes_sym <- c(
  "Aposymbiotic" = 8,
  "Facultative" = 20,
  "Symbiotic" = 18
)

p3 <- ggplot(pcoa_points3, aes(x = PCoA1, y = PCoA2, color = location)) +
  geom_point(aes(shape = symbiont.potential), size = 3, alpha = 0.8) +
  scale_color_manual(values = cols_location) +
  scale_shape_manual(values = shapes_sym) + 
  labs(
    x = paste0("PCoA1: (", var_explained[1], "%)"),
    y = paste0("PCoA2: (", var_explained[2], "%)"),
    color = "Location",
    shape = "Symbiont Potential"
  ) +
  # move PERMANOVA annotation to the bottom-left
  annotate(
    "text",
    x = -Inf, y = -Inf,
    label = stats_annotation,
    hjust = -0.05, vjust = -0.6,
    size = 4
  ) +
  theme_cowplot(font_size = 14) +
  theme(legend.position = "right")

################################################################################
## Permanova bleaching by symbiont potential - Figure S6C

keep_rows <- !is.na(meta2$bleaching) & !is.na(meta2$symbiont.potential)
meta_sub <- meta2[keep_rows, , drop = FALSE]
numeric_sub <- permanova_numeric_data2[keep_rows, , drop = FALSE]

bray_curtis_sub <- vegdist(numeric_sub, method = "bray")
bleach_sym_permanova <- adonis2(
  bray_curtis_sub ~ bleaching / symbiont.potential, 
  data = meta_sub, 
  permutations = 999
)

perm_p_val <- max(bleach_sym_permanova$`Pr(>F)`[1])
perm_r2    <- round(sum(bleach_sym_permanova$R2[1]), 3)
p_label    <- if(perm_p_val <= 0.001) "p < 0.001" else paste("p =", perm_p_val)
stats_annotation <- paste0("PERMANOVA: R² = ", perm_r2, ", ", p_label)

pcoa_res <- cmdscale(bray_curtis_sub, eig = TRUE, k = 2)
pcoa_pts <- as.data.frame(pcoa_res$points)
colnames(pcoa_pts) <- c("PCoA1", "PCoA2")
pcoa_pts <- bind_cols(pcoa_pts, meta_sub)

pos_eig <- pcoa_res$eig[pcoa_res$eig > 0]
var_explained <- round(100 * pcoa_res$eig[1:2] / sum(pos_eig), 1)

p_bleach_sym <- ggplot(pcoa_pts, aes(x = PCoA1, y = PCoA2, color = bleaching)) +
  geom_point(aes(shape = symbiont.potential), size = 3, alpha = 0.8) +
  scale_color_manual(values = cols_bleaching) +
  scale_shape_manual(values = shapes_sym) + 
  labs(
    x = paste0("PCoA1: (", var_explained[1], "%)"),
    y = paste0("PCoA2: (", var_explained[2], "%)"),
    color = "Bleaching Status",
    shape = "Symbiont Potential"
  ) +
  annotate(
    "text",
    x = -Inf, y = -Inf,
    label = stats_annotation,
    hjust = -0.05, vjust = -0.6,
    size = 4
  ) +
  theme_cowplot(font_size = 14) +
  theme(legend.position = "right")

################################################################################
pcoa_supp <- plot_grid(
  p2 + theme(legend.position="right"),
  p3 + theme(legend.position = "right"), 
  p_bleach_sym + theme(legend.position = "right"), 
  labels = c("A", "B", "C"), 
  label_size = 20,
  ncol = 1,              
  align = "v",
  axis = "lr"
)

ggsave(here("misc", "figs/pqn", "pcoa_supp.jpg"), 
       pcoa_supp, 
       width = 10, 
       height = 16, 
       dpi = 300)

################################################################################
################################################################################

## Bleaching status alone PERMANOVA - Table S3

desired_levels <- c("Bleached", "Non-Bleached", "Not Applicable")
permanova_numeric_data <- df %>% select(starts_with("x"))
keep_rows <- df$bleaching %in% desired_levels & complete.cases(permanova_numeric_data)

# subset numeric matrix and metadata
permanova_numeric_data2 <- permanova_numeric_data[keep_rows, , drop = FALSE]
meta_bleach <- df[keep_rows, ] %>%
  mutate(bleaching = factor(as.character(bleaching), levels = desired_levels))

# compute Bray–Curtis dissimilarity
bray_curtis_bleach <- vegdist(permanova_numeric_data2, method = "bray")
# PERMANOVA: bleaching alone
bleaching_permanova_result <- adonis2(
  bray_curtis_bleach ~ bleaching,
  data = meta_bleach,
  permutations = 999
)
print(bleaching_permanova_result)

pcoa_result <- cmdscale(bray_curtis_bleach, eig = TRUE, k = 2)
pcoa_points <- as.data.frame(pcoa_result$points)
colnames(pcoa_points) <- c("PCoA1", "PCoA2")
pcoa_points <- bind_cols(pcoa_points, meta_bleach)

r2_val <- round(bleaching_permanova_result$R2[1], 3)
p_val  <- bleaching_permanova_result$`Pr(>F)`[1]
p_lab <- if(p_val <= 0.001) "p < 0.001" else paste0("p = ", round(p_val, 3))
stats_label <- paste0("PERMANOVA: R² = ", r2_val, ", ", p_lab)

# compute percent variance explained using only positive eigenvalues
pos_eig <- pcoa_result$eig[pcoa_result$eig > 0]
var_explained <- round(100 * pcoa_result$eig[1:2] / sum(pos_eig), 1)

p4 <- ggplot(pcoa_points, aes(x = PCoA1, y = PCoA2, color = bleaching, fill = bleaching)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(
    geom = "polygon",
    alpha = 0.15,
    level = 0.95,
    type = "t",
    colour = NA
  ) + annotate(
    "text", 
    x = -Inf, y = -Inf, 
    label = stats_label, 
    hjust = -0.1, 
    vjust = -1.2, 
    size = 4, 
    fontface = "italic"
  ) +
  scale_color_manual(values = cols_bleaching) +
  scale_fill_manual(values = cols_bleaching) +
  labs(
    x = paste0("PCoA1: (", var_explained[1], "%)"),
    y = paste0("PCoA2: (", var_explained[2], "%)"),
    color = "Bleaching Status",
    fill = "Bleaching Status"
  ) +
  theme_cowplot(font_size = 14) +
  theme(legend.position = "right")
ggsave(here("misc", "figs/pqn", "pcoa_bleaching.jpg"), p4, width = 8, height = 6, dpi = 300)


################################################################################
## Symbiont potential alone PERMANOVA - Table S3

desired_levels <- na.omit(unique(as.character(df$symbiont.potential)))
permanova_numeric_data <- df %>% select(starts_with("x"))
keep_rows <- !is.na(df$symbiont.potential) & complete.cases(permanova_numeric_data)

permanova_numeric_data2 <- permanova_numeric_data[keep_rows, , drop = FALSE]
desired_levels <- c("Aposymbiotic", "Facultative", "Symbiotic")
meta_sym <- df[keep_rows, ] %>%
  mutate(symbiont.potential = factor(as.character(symbiont.potential), 
                                     levels = desired_levels))

# Bray–Curtis dissimilarity
bray_curtis_sym <- vegdist(permanova_numeric_data2, method = "bray")

symbiont_permanova_result <- adonis2(
  bray_curtis_sym ~ symbiont.potential,
  data = meta_sym,
  permutations = 999
)
print(symbiont_permanova_result)

pcoa_result <- cmdscale(bray_curtis_sym, eig = TRUE, k = 2)
pcoa_points <- as.data.frame(pcoa_result$points)
colnames(pcoa_points) <- c("PCoA1", "PCoA2")
pcoa_points <- bind_cols(pcoa_points, meta_sym)

r2_val <- round(symbiont_permanova_result$R2[1], 3)
p_val  <- symbiont_permanova_result$`Pr(>F)`[1]

p_lab <- if(p_val <= 0.001) "p < 0.001" else paste0("p = ", round(p_val, 3))
stats_label <- paste0("PERMANOVA: R² = ", r2_val, ", ", p_lab)

p5 <- ggplot(pcoa_points, aes(x = PCoA1, y = PCoA2, color = symbiont.potential, fill = symbiont.potential)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(
    geom = "polygon",
    alpha = 0.15,
    level = 0.95,
    type = "t",
    colour = NA
  ) +
  scale_color_manual(values = cols_symbiont) +
  scale_fill_manual(values = cols_symbiont) +
  annotate(
    "text", 
    x = -Inf, y = -Inf, 
    label = stats_label, 
    hjust = -0.1, 
    vjust = -1.2, 
    size = 4, 
    fontface = "italic"
  ) +
  labs(
    x = paste0("PCoA1: (", var_explained[1], "%)"),
    y = paste0("PCoA2: (", var_explained[2], "%)"),
    color = "Symbiont Potential",
    fill = "Symbiont Potential"
  ) +
  theme_cowplot(font_size = 14) +
  theme(legend.position = "right")

################################################################################
## Location alone PERMANOVA - Table S3

desired_levels <- na.omit(unique(as.character(df$location)))
permanova_numeric_data <- df %>% select(starts_with("x"))
keep_rows <- !is.na(df$location) & complete.cases(permanova_numeric_data)

permanova_numeric_data2 <- permanova_numeric_data[keep_rows, , drop = FALSE]
meta_loc <- df[keep_rows, ] %>%
  mutate(location = factor(as.character(location), levels = desired_levels))

bray_curtis_loc <- vegdist(permanova_numeric_data2, method = "bray")

location_permanova_result <- adonis2(
  bray_curtis_loc ~ location,
  data = meta_loc,
  permutations = 999
)
print(location_permanova_result)

pcoa_result <- cmdscale(bray_curtis_loc, eig = TRUE, k = 2)
pcoa_points <- as.data.frame(pcoa_result$points)
colnames(pcoa_points) <- c("PCoA1", "PCoA2")
pcoa_points <- bind_cols(pcoa_points, meta_loc)

r2_val <- round(location_permanova_result$R2[1], 3)
p_val  <- location_permanova_result$`Pr(>F)`[1]

p_lab <- if(p_val <= 0.001) "p < 0.001" else paste0("p = ", round(p_val, 3))
stats_label <- paste0("PERMANOVA: R² = ", r2_val, ", ", p_lab)

pos_eig <- pcoa_result$eig[pcoa_result$eig > 0]
var_explained <- round(100 * pcoa_result$eig[1:2] / sum(pos_eig), 1)

p6 <- ggplot(pcoa_points, aes(x = PCoA1, y = PCoA2, color = location, fill = location)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(
    geom = "polygon",
    alpha = 0.3,
    level = 0.95,
    type = "t",
    colour = NA
  ) +
  annotate(
    "text", 
    x = -Inf, y = -Inf, 
    label = stats_label, 
    hjust = -0.1, 
    vjust = -1.2, 
    size = 4, 
    fontface = "italic"
  ) +
  scale_color_manual(values = cols_location) +
  scale_fill_manual(values = cols_location) +
  labs(
    x = paste0("PCoA1: (", var_explained[1], "%)"),
    y = paste0("PCoA2: (", var_explained[2], "%)"),
    color = "Location",
    fill = "Location"
  ) +
  theme_cowplot(font_size = 14) +
  theme(legend.position = "right")

################################################################################
################################################################################

# Box plots for Figure 2C-E

metabolite_data <- df %>%
  select(starts_with("x"))

########### abundance - unused
avg_abundance <- metabolite_data %>%
  mutate(avg_abundance = rowMeans(., na.rm = TRUE)) %>%
  select(avg_abundance)

plot_data <- df %>%
  select(host_order) %>%
  bind_cols(avg_abundance) %>%
  mutate(
    group = ifelse(host_order == "Scleractinia",
                   "1",
                   "0")
  )
plot_data$group[is.na(plot_data$group)] <- 0
unique(plot_data$group)

plot_data_clean <- plot_data %>%
  mutate(group = factor(group)) %>% # Ensure it's a factor first
  mutate(group = fct_recode(group, 
                            "Scleractinia" = "1", 
                            "Other" = "0")) %>%
  mutate(group = fct_relevel(group, "Scleractinia", "Other"))
cols_sclero_named <- c("Scleractinia" = "#DE7862FF", "Other" = "#D8AF39FF")

########### ubiquity - Figure 2E ###################

plot_data_clean <- plot_data %>%
  mutate(group = factor(group)) %>% 
  mutate(group = fct_recode(group, 
                            "Scleractinia" = "1", 
                            "Other" = "0")) %>%
  mutate(group = fct_relevel(group, "Scleractinia", "Other"))
cols_sclero_named <- c("Scleractinia" = "#DE7862FF", "Other" = "#D8AF39FF")

ubiquity_values <- metabolite_data %>%
  mutate(avg_ubiquity = rowMeans(. > 0, na.rm = TRUE)) %>%
  pull(avg_ubiquity) 
plot_data_clean$avg_ubiquity <- ubiquity_values

p_ubiquity <- ggplot(plot_data_clean, aes(x = group, y = avg_ubiquity, fill = group, color = group)) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5, show.legend = FALSE) +
  geom_boxplot(alpha = 0.7, width = 0.6, outlier.shape = NA, color = "black") + 
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.format",
    label.x = 1.4,
    symnum.args = list(
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
      symbols = c("****", "***", "**", "*", "ns")
    )
  ) +
  scale_fill_manual(values = cols_sclero_named) +
  scale_color_manual(values = cols_sclero_named) +
  labs(
    y = "Metabolite Ubiquity"
  ) +
  theme_pubr(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.title.x = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
print(p_ubiquity)

########### Richness - figure 2C ################### 
richness_values <- rowSums(metabolite_data > 0, na.rm = TRUE)
plot_data_clean$richness <- richness_values

p_richness <- ggplot(plot_data_clean, aes(x = group, y = richness, fill = group, color = group)) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5, show.legend = FALSE) +
  geom_boxplot(alpha = 0.7, width = 0.6, outlier.shape = NA, color = "black") + 
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.format",
    label.x = 1.4,
    symnum.args = list(
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
      symbols = c("****", "***", "**", "*", "ns")
    )
  ) +
  scale_fill_manual(values = cols_sclero_named) +
  scale_color_manual(values = cols_sclero_named) +
  labs(
    y = "Metabolite Richness"
  ) +
  theme_pubr(base_size = 13) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
print(p_richness)

############ shannon diversity - Figure 2D
shannon_values <- diversity(metabolite_data, index = "shannon")
plot_data_clean$shannon  <- shannon_values

p_entropy <- ggplot(plot_data_clean, aes(x = group, y = shannon, fill = group, color = group)) +
  geom_jitter(width = 0.2, alpha = 0.5, size = 1.5, show.legend = FALSE) +
  geom_boxplot(alpha = 0.7, width = 0.6, outlier.shape = NA, color = "black") + 
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.format",
    label.x = 1.4,
    symnum.args = list(
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
      symbols = c("****", "***", "**", "*", "ns")
    )
  ) +
  scale_fill_manual(values = cols_sclero_named) +
  scale_y_continuous(breaks = c(5, 6, 7), limits = c(4.5, 7.5)) +
  scale_color_manual(values = cols_sclero_named) +
  labs(
    y = "Shannon Entropy"
  ) +
  theme_pubr(base_size = 13) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )
print(p_entropy)

################################################################################
################################################################################
# Dendrograms - Figure 2B

avg_metabolite_values_family <- df |> 
  group_by(host_family) |> 
  filter(!is.na(host_family)) |> 
  summarise(across(starts_with("x"), mean, na.rm = TRUE)) |>
  ungroup()

mat_family <- avg_metabolite_values_family %>%
  as.data.frame()
rownames(mat_family) <- mat_family$host_family
mat_family <- mat_family[ , setdiff(names(mat_family), "host_family")]
mat_family[] <- lapply(mat_family, as.numeric)
mat_family <- mat_family %>% mutate_all(~ifelse(is.na(.), 0, .))

bray_curtis_family <- vegdist(sqrt(mat_family), method = 'bray')
cluster.average <- hclust(bray_curtis_family, method = 'ward.D2') 
#method = average uses upgma

dend_data <- dendro_data(cluster.average, type = "rectangle")

branch_palette <- c(
  "Ochrophyta"         = "#32CD32", 
  "Chlorophyta"        = "#00CED1",
  "Cnidaria"           = "#304530",
  "Chordata"           = "#9370DB", 
  "Porifera"           = "#ec93ed", 
  "Scleractinia_Color" = "#DE7862FF", # Scleractinia
  "grey40"             = "grey40"   # Internal nodes
)

#################################### with n=

family_counts <- df %>%
  filter(!is.na(host_family)) %>%
  group_by(host_family) %>%
  summarise(n = n()) %>%
  ungroup()

leaf_metadata <- df %>%
  filter(!is.na(host_family)) %>%
  group_by(host_family) %>%
  summarize(
    phylum = first(host_phylum),
    is_scler = first(host_order) == "Scleractinia"
  ) %>%
  left_join(family_counts, by = "host_family") %>%
  mutate(label_with_n = paste0(host_family, " (n=", n, ")")) %>%
  ungroup()

# function for visualization -- rotating internal nodes
optimize_hclust_by_phylum <- function(hc, group_lookup, max_iter = 100) {
  
  n <- length(hc$labels)
  
  get_order_groups <- function(hc) {
    labs <- hc$labels[hc$order]
    group_lookup[labs]
  }
  
  score_order <- function(hc) {
    labs <- hc$labels[hc$order]
    groups <- group_lookup[labs]
    
    phylum_order <- c(
      "Scleractinia_Color",
      "Cnidaria",
      "Porifera",
      "Chordata",
      "Chlorophyta",
      "Ochrophyta"
    )
    
    ranks <- match(groups, phylum_order)
    r <- rle(groups)
    block_counts <- table(r$values)
    split_penalty <- sum((block_counts - 1)^2, na.rm = TRUE)
    inversion_penalty <- sum(diff(ranks) < 0, na.rm = TRUE)
    
    singleton_penalty <- 0
    for (g in unique(groups)) {
      if (sum(groups == g, na.rm = TRUE) == 1) {
        pos <- which(groups == g)
        singleton_penalty <- singleton_penalty + min(pos - 1, length(groups) - pos)
      }
    }
    
    split_penalty * 10000 +
      inversion_penalty * 100 +
      singleton_penalty
  }
  
  recompute_order <- function(hc) {
    n <- length(hc$labels)
    
    get_node_order <- function(node) {
      if (node < 0) {
        return(-node)
      }
      
      left <- hc$merge[node, 1]
      right <- hc$merge[node, 2]
      
      c(get_node_order(left), get_node_order(right))
    }
    
    hc$order <- get_node_order(n - 1)
    hc
  }
  
  hc <- recompute_order(hc)
  
  for (iter in seq_len(max_iter)) {
    improved <- FALSE
    
    for (node in seq_len(nrow(hc$merge))) {
      current_score <- score_order(hc)
      
      trial <- hc
      trial$merge[node, ] <- rev(trial$merge[node, ])
      trial <- recompute_order(trial)
      
      trial_score <- score_order(trial)
      
      if (trial_score <= current_score) {
        hc <- trial
        improved <- TRUE
      }
    }
    
    if (!improved) break
  }
  
  hc
}

leaf_metadata <- leaf_metadata %>%
  mutate(
    phylum_group = if_else(is_scler, "Scleractinia_Color", phylum)
  )

group_lookup <- setNames(
  leaf_metadata$phylum_group,
  leaf_metadata$host_family
)

cluster.average <- hclust(bray_curtis_family, method = "average")

cluster.average <- optimize_hclust_by_phylum(
  cluster.average,
  group_lookup,
  max_iter = 10
)

dend <- as.dendrogram(cluster.average)
dend_data <- dendro_data(dend, type = "rectangle")

#######

dend_segments <- dend_data$segments %>%
  left_join(dend_data$labels %>% select(x, label), by = "x") %>% 
  left_join(leaf_metadata, by = c("label" = "host_family")) %>%
  mutate(
    branch_color = case_when(
      is_scler == TRUE ~ "Scleractinia_Color",
      !is.na(phylum)   ~ phylum,
      TRUE             ~ "grey40" 
    )
  )

dend_labels_phylum <- dend_data$labels %>%
  left_join(leaf_metadata, by = c("label" = "host_family")) %>%
  mutate(
    text_color_group = if_else(is_scler == TRUE, "Scleractinia_Color", phylum)
  )

# Save fig 2b dendrogram
p_dendro <- ggplot() +
  geom_segment(data = dend_segments, 
               aes(x = x, y = y, xend = xend, yend = yend, color = branch_color),
               linewidth = 0.8) +
  geom_text(data = dend_labels_phylum, 
            aes(x = x, y = y, label = label_with_n, color = text_color_group),
            hjust = -0.1, 
            size = 5, 
            show.legend = FALSE) +
  coord_flip() +
  scale_y_reverse(
    expand = c(0.5, 0), # Slightly increased expansion to fit longer text
    breaks = seq(0, 5, 1)
  ) +
  scale_color_manual(
    values = branch_palette,
    breaks = c("Scleractinia_Color","Chlorophyta", "Chordata", "Cnidaria", "Ochrophyta","Porifera"),
    labels = c("Scleractinia","Chlorophyta", "Chordata", "Cnidaria", "Ochrophyta","Porifera")
  ) +
  theme_pubr() +
  labs(y = "Bray-Curtis Distance", x = "", color = "Classification") +
  theme(
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = c(0.10, 0.5),   
    legend.background = element_blank(), 
    legend.box.background = element_blank(),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 16),
  ) +
  guides(color = guide_legend(
    override.aes = list(alpha = 1, size = 4, shape = 15)
  ))

################################################################################
# combine plots with cowplot to make final Figure 2

legend_sclero <- get_legend(
  p + 
    guides(color = guide_legend(title = "Order"), fill = guide_legend(title = "Order")) +
    theme(
      legend.position = "left", 
      legend.box.margin = margin(6, 0, 0, 6),
      legend.title = element_text(size = 24, face = "bold"), 
      legend.text = element_text(size = 20),
      legend.title.align = 0.5,
      legend.text.align = 0.5,
      legend.key.size = unit(1.2, "lines")))
clean_theme <- theme(
  legend.position = "none",
  plot.title = element_blank() 
)

p <- p + theme(legend.position = "none")
p_ubiquity_clean  <- p_ubiquity + clean_theme
p_entropy_clean  <- p_entropy + clean_theme
p_richness_clean  <- p_richness + clean_theme

row1 <- plot_grid(
  p, p_dendro,
  ncol = 2,
  rel_widths = c(1, 1.4), 
  labels = c("A", "B"),
  label_size = 24
)

row2 <- plot_grid(
  p_richness_clean,p_entropy_clean,p_ubiquity_clean,
  nrow = 1,
  rel_widths = c(0.6, 0.6, 0.6, 0.6), 
  align = 'h', axis = 'tb',
  labels = c("C", "D", "E"),
  label_size = 24
)

final_plot <- plot_grid(
  row1, 
  row2, 
  ncol = 1, 
  rel_heights = c(1.4, 1) 
)
final_plot_legend <- plot_grid(
  legend_sclero,
  final_plot,
  ncol = 2, 
  rel_widths = c(0.15, 1) 
)

print(final_plot)
ggsave(here("misc", "figs/pqn", "fig2.pdf"), final_plot_legend, width = 18, height = 12, dpi = 300)


################################################################################
################################################################################
# Curacao only - Figure S5

df_cur <- df %>% 
  filter(location == "Curaçao", !is.na(scleractinia))

num_cur  <- df_cur %>% select(starts_with("x"))
meta_cur <- df_cur

bc_cur <- vegdist(num_cur, method = "bray")
perm_cur <- adonis2(bc_cur ~ scleractinia, data = meta_cur, permutations = 999)

p_val_cur <- perm_cur[["Pr(>F)"]][1]
r2_cur    <- round(perm_cur[["R2"]][1], 3)
p_lab_cur <- if(is.na(p_val_cur)) "p = NA" else if(p_val_cur <= 0.001) "p < 0.001" else paste("p =", p_val_cur)
stats_cur <- paste0("PERMANOVA: R² = ", r2_cur, ", ", p_lab_cur)

pcoa_cur <- cmdscale(bc_cur, eig = TRUE, k = 2)
df_cur   <- cbind(as.data.frame(pcoa_cur$points), meta_cur)
var_cur  <- round(100 * pcoa_cur$eig[1:2] / sum(pcoa_cur$eig[pcoa_cur$eig > 0]), 1)

p_cur <- ggplot(df_cur, aes(x = V1, y = V2)) +
  geom_point(
    aes(
      color = as.factor(scleractinia), 
      shape = "Curaçao"  
    ), 
    size = 3, 
    alpha = 0.7
  ) +
  scale_color_manual(
    values = cols_sclero, 
    labels = c("1" = "Scleractinia", "0" = "Other"),
    name = "Order"
  ) +
  scale_shape_manual(
    values = c("Curaçao" = 16), 
    name = "Location" 
  ) +
  labs(
    x = paste0("PCoA1: (", var_cur[1], "%)"), 
    y = paste0("PCoA2: (", var_cur[2], "%)")
  ) +
  annotate(
    "text", x = Inf, y = -Inf, 
    label = stats_cur, 
    hjust =2.55, vjust = -1.2, 
    size = 3.5, fontface = "italic"
  ) +
  theme_cowplot()
p_cur
ggsave(here("misc", "figs/pqn", "pcoa_cur_scler.jpg"), p_cur, width = 8, height = 6, dpi = 300)


################################### permanova outputs across various runs

tidy_permanova <- function(model_output, model_name) {
  model_output %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Term") %>%
    filter(!Term %in% c("Residual", "Total")) %>%
    mutate(Model = model_name) %>%
    select(Model, Term, Df, SumOfSqs, R2, F, `Pr(>F)`)
}

table_sclero <- tidy_permanova(scleractinia_permanova_result, "Scleractinia")
table_cur <- tidy_permanova(perm_cur, "Scleractinia - Curacao")
table_locB   <- tidy_permanova(locB_permanova_result, "Location x Bleaching")
table_locS   <- tidy_permanova(locS_permanova_result, "Location x Symbiont")
table_bSym <-  tidy_permanova(bleach_sym_permanova, "Bleaching x Symbiont")
table_loc <- tidy_permanova(location_permanova_result, "Location")
table_b <- tidy_permanova(bleaching_permanova_result, "Bleaching")
table_sym <-  tidy_permanova(symbiont_permanova_result, "Symbiont")

full_stats_table <- bind_rows(table_sclero, table_cur, table_locB, table_locS, 
                              table_bSym, table_sym, table_loc, table_b)
full_stats_table <- full_stats_table %>%
  mutate(
    across(c(SumOfSqs, R2, F), ~ round(., 3)),
    `Pr(>F)` = ifelse(`Pr(>F)` <= 0.001, "< 0.001", as.character(`Pr(>F)`))
  )
print(full_stats_table)
write.csv(full_stats_table, here("misc", "PERMANOVA_pqn_results.csv"))