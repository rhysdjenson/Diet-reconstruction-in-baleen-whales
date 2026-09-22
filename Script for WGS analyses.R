library(ggplot2)
library(simmr)
library(dplyr)
library(tidyr)
library(ggnewscale)
library(readxl)
library(compositions)
library(vegan)
library(scales)
library(reshape2)
library(colorspace)

# Stacked bar chart of prey composition of fecal samples based on whole-genome sequencing

df <- read_excel("combined_diet_percent.xlsx")

# Define krill species
krill_species <- c("Euphausia pacifica", "Thysanoessa spinifera", "Thysanoessa raschii")
df_long <- df %>%
  pivot_longer(
    cols = matches("^S\\d+$"),
    names_to = "SampleID",
    values_to = "Proportion"
  ) %>%
  mutate(
    Proportion = Proportion / 100,
    Group = ifelse(Species %in% krill_species, "Krill", "Fish")
  )

df_long <- df_long %>%
  mutate(
    SampleID = factor(SampleID, levels = paste0("S", 1:20))
  )


krill_colors <- c(
  "Euphausia pacifica"    = "#e08a8a",
  "Thysanoessa spinifera" = "#e8a3a3", 
  "Thysanoessa raschii"   = "#eeb8ae"   
)

species_order <- unique(df_long$Species)
other_species <- setdiff(species_order, krill_species)

herring <- "Clupea pallasii"
remaining_species <- setdiff(other_species, herring)
n_remaining <- length(remaining_species)

base_hues <- seq(40, 340, length.out = n_remaining) 

fish_colors <- character(n_remaining)
for (i in seq_len(n_remaining)) {
  lum <- if (i %% 2 == 0) 60 else 78  
  fish_colors[i] <- hex(polarLUV(L = lum, C = 55, H = base_hues[i]))
}
fish_colors <- setNames(fish_colors, remaining_species)

fish_colors[herring] <- "#6ca0a3"

fish_totals <- df_long %>%
  filter(!Species %in% krill_species) %>%
  group_by(Species) %>%
  summarise(Total = sum(Proportion, na.rm = TRUE)) %>%
  arrange(desc(Total))

last_species <- c("Sebastes alutus", "Sebastes caurinus", "Sebastes melanops")

fish_order <- fish_totals$Species
fish_order <- c(setdiff(fish_order, last_species), last_species)

new_levels <- c(krill_species, fish_order)
df_long$Species <- factor(df_long$Species, levels = new_levels)

df_long <- df_long %>%
  arrange(SampleID, Species) %>%
  group_by(SampleID) %>%
  mutate(
    ymax = cumsum(Proportion),
    ymin = lag(ymax, default = 0),
    xnum = as.numeric(SampleID)
  ) %>%
  ungroup()

krill_df <- df_long %>% filter(Group == "Krill")
fish_df  <- df_long %>% filter(Group == "Fish")

sample_levels <- levels(df_long$SampleID)

p_main <- ggplot() +
  geom_rect(
    data = krill_df,
    aes(xmin = xnum - 0.42, xmax = xnum + 0.42, ymin = ymin, ymax = ymax, fill = Species),
    colour = "black", linewidth = 0.4
  ) +
  scale_fill_manual(
    values = krill_colors,
    guide = guide_legend(title = "Krill", order = 1)
  ) +
  new_scale_fill() +
  geom_rect(
    data = fish_df,
    aes(xmin = xnum - 0.42, xmax = xnum + 0.42, ymin = ymin, ymax = ymax, fill = Species),
    colour = "black", linewidth = 0.4
  ) +
  scale_fill_manual(
    values = fish_colors,
    guide = guide_legend(title = "Fish", order = 2)
  ) +
  scale_x_continuous(
    breaks = seq_along(sample_levels),
    labels = sample_levels,
    expand = expansion(mult = c(0.0085, 0.0085))   # small gap on left, none on right
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Individual fecal sample", y = "Relative abundance") +
  theme_classic(base_family = "Arial") +
  theme(
    axis.text.x  = element_text(size = 14, colour = "black", face = "plain"),
    axis.text.y  = element_text(size = 14, colour = "black", face = "plain"),
    axis.title.x = element_text(size = 20, face = "plain"),
    axis.title.y = element_text(size = 20, face = "plain"),
    legend.title = element_text(size = 20, face = "plain"),
    legend.text  = element_text(size = 12, face = "italic", colour = "black"),
    panel.grid   = element_blank(),
    axis.line    = element_line(linewidth = 0.4),
    legend.spacing.y = unit(-5, "pt"),     # padding *around* each legend
    legend.box.spacing = unit(4, "pt")          # space between legend block and plot panel
  )

p_main

p_main <- p_main +
  theme(text = element_text(family = "sans", size = 16))
ggsave(
  filename = "stackedWGS.pdf",
  plot = p_main,
  device = "pdf",
  width = 12,
  height = 6,
  units = "in",
  encoding = "WinAnsi.enc"
)

# Generate a principle coordinate analysis (PCoA)

pcoadf <- read_excel("combined_diet_percent.xlsx")
pcoadf <- as.data.frame(pcoadf)
rownames(pcoadf) <- pcoadf$Species
pcoadf <- pcoadf[ , -1]   # drop the Species column
pcoadf[] <- lapply(pcoadf, as.numeric)
pcoadf<- t(pcoadf)
abund <- as.matrix(pcoadf) #Convert to matrix
str(abund)
meta <- read_excel("combined_diet_percent.xlsx", sheet = 2)  #Read the metadata sheet
meta <- as.data.frame(meta)
rownames(meta) <- meta$`Sample ID`
meta <- meta[ , -which(names(meta) == "Sample ID")]
identical(rownames(meta), rownames(abund))
keep     <- colSums(abund > 0) >= round(0.1 * nrow(abund))
abund_f  <- abund[, keep, drop = FALSE]
abund_pc <- abund_f + 1L                    # pseudocount to remove zeros
abund_clr <- clr(acomp(abund_pc))           # CLR (rows = samples)
dist_ait  <- dist(abund_clr, method = "euclidean")
pcoa <- cmdscale(dist_ait, k = 2, eig = TRUE) ## PCoA and variance
eig <- pcoa$eig
eig_pos <- eig[eig > 0]
var_explained <- eig / sum(eig_pos)         # proportion per axis
vx1 <- 100 * var_explained[1]
vx2 <- 100 * var_explained[2]
ord_df <- data.frame(
  Axis1 = pcoa$points[, 1],
  Axis2 = pcoa$points[, 2]
)
rownames(ord_df) <- rownames(abund_clr)
ord_df$Month  <- meta[rownames(ord_df), "Month"] # bring in the month as a proper column
ord_df$Group <- meta[rownames(ord_df), "Group"]  


ord_df$Month <- factor(ord_df$Month, levels = c("July", "August", "September"))


month_cols <- c(
  "July"       = "#66c2a5",
  "August"     = "#8da0cb",
  "September"  = "#fc8d62"
)

# 1) PERMANOVA
adon <- adonis2(dist_ait ~ Month, data = meta, permutations = 999)
r2_month <- adon$R2[1]
p_month  <- adon$`Pr(>F)`[1]
f_month <- adon$F[1]

# 2) PERMDISP
bd   <- betadisper(dist_ait, meta$Month)
bd_p <- permutest(bd, permutations = 999)$tab[1, "Pr(>F)"]

# 3) Annotation label
lab <- sprintf(
  "PERMANOVA (Month): R² = %.3f, p = %.3f\nPERMDISP: p = %.3f",
  r2_month, p_month, bd_p
)

# 4) PCoA plot with red/blue/green dashed ellipses
p1 <- ggplot(ord_df, aes(Axis1, Axis2)) +
  
  # 95% solid ellipses with light fill
  stat_ellipse(
    aes(group = Month, color = Month, fill = Month),
    geom = "polygon",
    type = "t",
    level = 0.95,
    linetype = "solid",
    linewidth = 0.5,
    alpha = 0.08,
    show.legend = FALSE
  ) +
  
  # Points
  geom_point(
    aes(fill = Month, shape = Month),
    size = 3,
    stroke = 0.5
  ) +
  
  scale_shape_manual(values = c("July" = 21, "August" = 24, "September" = 22), name = "Month") +
  scale_color_manual(values = month_cols) +
  scale_fill_manual(values = month_cols, name = "Month") +
  
  labs(
    x = sprintf("PCoA 1 (%.1f%%)", vx1),
    y = sprintf("PCoA 2 (%.1f%%)", vx2)
  ) +
  coord_cartesian(ylim = c(-7, 7)) +
  
  # Month R2/p annotation, bottom-left corner
  annotate(
    "text",
    x = -Inf, y = -Inf,
    label = sprintf(
      "Month: F = %s, R² = %s, p = %s",
      formatC(f_month, digits = 2, format = "g"),
      formatC(r2_month, digits = 2, format = "g"),
      formatC(p_month, digits = 2, format = "g")
    ),
    hjust = -0.05, vjust = -0.5,
    size = 5
  ) +
  
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = c(0, 1),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.margin = margin(4, 6, 4, 6),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 20),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16)
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 4))
  )

p1

p1 <- p1 +
  theme(text = element_text(family = "sans", size = 16))

ggsave(
  filename = "PCOA.pdf",
  plot = p1,
  device = "pdf",
  width = 8,
  height = 7,
  units = "in",
  encoding = "WinAnsi.enc"
)



