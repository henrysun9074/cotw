library(ggplot2)
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires) # need to install using devtools
library(patchwork)
library(ggspatial)
library(cowplot)
library(ggpubr)
library(maptiles)   # for satellite/terrain tiles
library(tidyterra)  # for geom_spatraster_rgb (maptiles output)

# sites
sites_df <- data.frame(
  location = c(rep("Sri Lanka (n=138)", 5), "Hawaii (n=77)",
               rep("Curaçao (n=272)", 3),
               "North Carolina (n=81)"),
  site_name = c("KYK 1", "KYK2B", "KYK 2", "KCH", "PIN 3",
                "Kaneohe Bay",
                "Piscaderabaai", "Blaubaai", "Water Factory",
                "Radio Island"),
  lat = c(8.009183, 8.011181, 8.010667, 8.861336, 8.723558,
          21.473743, 12.122715, 12.134190, 12.108721, 34.705736),
  lon = c(81.520250, 81.520194, 81.520250, 81.084431, 81.204808,
          -157.812254, -68.970977, -68.985736, -68.953637, -76.678666)
)

sites_sf <- st_as_sf(sites_df, coords = c("lon", "lat"), crs = 4326)

cols_location <- c("Curaçao (n=272)"       = "#002594FF",
                   "Hawaii (n=77)"          = "#E0B2CDFF",
                   "North Carolina (n=81)"  = "#54C4E3FF",
                   "Sri Lanka (n=138)"      = "#F3AA4FFF")

sites_global <- sites_sf %>%
  group_by(location) %>%
  summarize(geometry = st_centroid(st_union(geometry))) %>%
  ungroup()
world      <- ne_countries(scale = "large", returnclass = "sf") 
coastlines <- ne_coastline(scale = "large", returnclass = "sf")

# Draw a Robinson-projected ellipse to clip the globe to an oval shape
robin_crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

world_robin      <- st_transform(world,        crs = robin_crs)
sites_global_rob <- st_transform(sites_global, crs = robin_crs)

# Graticule ellipse 
center <- st_sfc(st_point(c(0, 0)), crs = robin_crs)

oval <- center %>%
  st_buffer(dist = 1) %>%
  st_sf() %>%
  mutate(geometry = geometry * matrix(c(18500000, 0, 0, 10000000), 2, 2)) %>%
  st_set_crs(robin_crs)

p_main <- ggplot() +
  geom_sf(data = oval, fill = "#e7f6ff", color = NA) +
  geom_sf(data = world_robin, fill = "gray60", color = "gray80",
          linewidth = 0.2, alpha = 0.7) +
  geom_sf(data = oval, fill = NA, color = "gray20", linewidth = 0.3, alpha = 0.3) +
  geom_sf(data = sites_global_rob, size = 6, shape = 21,
          fill = NA, color = "black", stroke = 1.2) +
  geom_sf(data = sites_global_rob, aes(color = location), size = 4.5,
          shape = 16, alpha = 0.95) +
  scale_color_manual(values = cols_location) +
  coord_sf(crs = robin_crs, xlim = c(-17005833, 17005833),
           ylim = c(-8625155, 8625155), expand = FALSE, datum = NA) +
  theme_void() +
  theme(plot.background  = element_rect(fill = "white", color = NA),
        legend.position  = "none",
        panel.background = element_rect(fill = "white", color = NA))
p_main

plot_inset <- function(loc_name, x_lim, y_lim, breaks_n = 3) {
  
  dat_sub <- sites_sf %>% filter(location == loc_name)
  
  bbox <- dat_sub %>%
    st_buffer(dist = 0.5) %>%  
    st_bbox()
  bbox["xmin"] <- x_lim[1]
  bbox["xmax"] <- x_lim[2]
  bbox["ymin"] <- y_lim[1]
  bbox["ymax"] <- y_lim[2]
  
  tiles <- get_tiles(bbox, provider = "Esri.WorldImagery",
                     zoom = 10, cachedir = tempdir(), verbose = FALSE)
  
  sf_use_s2(FALSE)
  world_valid <- st_make_valid(world)
  world_crop  <- st_crop(world_valid, bbox)
  
  ggplot() +
    geom_spatraster_rgb(data = tiles) +
    geom_sf(data = world_crop,
            fill = NA, color = "white", linewidth = 0.3, alpha = 0.4) +
    geom_sf(data = dat_sub, aes(color = location),
            size = 5, shape = 21,
            fill   = cols_location[loc_name],
            color  = "white",
            stroke = 1.5) +
    annotation_scale(location = "bl", width_hint = 0.35,
                     text_size = 9, text_col = "white",
                     line_col  = "white", bar_cols = c("white", "gray40")) +
    coord_sf(xlim = x_lim, ylim = y_lim, expand = FALSE, crs = 4326,
             clip = "on") +
    scale_x_continuous(breaks = scales::extended_breaks(n = breaks_n)(x_lim) %>%
                         .[. > x_lim[1] & . < x_lim[2]]) +
    scale_y_continuous(breaks = scales::extended_breaks(n = breaks_n)(y_lim) %>%
                         .[. > y_lim[1] & . < y_lim[2]]) +
    theme_bw(base_size = 10) +
    theme(
      axis.title        = element_blank(),
      axis.text.x       = element_text(size = 8, color = "black",
                                       margin = margin(t = 4)),
      axis.text.y       = element_text(size = 8, color = "black",
                                       margin = margin(r = 4)),
      plot.title        = element_blank(),
      panel.border      = element_rect(color = "black", linewidth = 1.5,
                                       fill = NA),
      panel.background  = element_rect(fill = NA, color = NA),
      panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      legend.position   = "none",
      plot.background   = element_rect(fill = NA, color = NA),
      plot.margin       = margin(10, 10, 10, 10),
      axis.ticks        = element_line(color = "black"),
      axis.ticks.length = unit(4, "pt")
    )
}

p_hi_inset  <- plot_inset("Hawaii (n=77)",         c(-158.3, -157.5), c(21.2, 21.8))
p_sl_inset  <- plot_inset("Sri Lanka (n=138)",      c(80.5,  82.0),   c(7.5,  9.5))
p_nc_inset  <- plot_inset("North Carolina (n=81)",  c(-77.5, -75.5),  c(34.5, 35.5))
p_cur_inset <- plot_inset("Curaçao (n=272)",        c(-69.15,-68.80), c(12.05,12.35))

##################### make final map for Figure 1 and save
final_map <- p_main +
  inset_element(p_hi_inset,  left = 0.05, bottom = 0.10, right = 0.35, top = 0.45) +
  inset_element(p_nc_inset, left = 0.45, bottom = 0.60, right = 0.68, top = 0.92) +
  inset_element(p_cur_inset,  left = 0.40, bottom = 0.15, right = 0.65, top = 0.50) +
  inset_element(p_sl_inset,  left = 0.70, bottom = 0.25, right = 1, top = 0.65)
final_map

ggsave(here("misc", "figs", "map.pdf"), plot = final_map,
       width = 14, height = 10, dpi = 300, bg = "white")

