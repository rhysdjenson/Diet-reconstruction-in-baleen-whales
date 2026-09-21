library(ggplot2)
library(simmr)
library(dplyr)
library(tidyr)
library(scales)
library(reshape2)



# Isotopic baseline mean and SD for each location
bl_SoG23  <- list(d13C_mean = -19.61, d13C_sd = 1.35, n = 63,
                  d15N_mean =  9.3,  d15N_sd = 1.34, n2 = 63) # Pacific krill from Strait of Georgia, 2023
bl_JdF24  <- list(d13C_mean = -17.755, d13C_sd = 0.2474874, n = 2,
                  d15N_mean =  8.945,  d15N_sd = 0.1484924, n2 = 2) # Pacific krill from Juan de Fuca Strait, 2024

# Isotope values for humpback whale prey collected from Strait of Georgia, 2023
krill_SoG23   <- list(d13C_mean = -19.61, d13C_sd = 1.35, n = 63,
                      d15N_mean =  9.3,  d15N_sd = 1.34, n2 = 63)

herring_SoG23 <- list(d13C_mean = -16.61, d13C_sd = 1.05, n = 23,
                      d15N_mean = 13.46,  d15N_sd = 0.42, n2 = 23)

juv_herring_SoG23 <- list(d13C_mean = -18.25, d13C_sd = 1.12, n = 202,
                          d15N_mean = 12.54,  d15N_sd = 0.33, n2 = 202)

# Custom function: shifts the isotopic values of prey collected in SOG23 to reflect the isotopic basline of JdF24 
mk_adj <- function(prey, bl_SoG, bl_JdF) {
  # SEs of baseline means
  bl_SoG_d13C_se <- bl_SoG$d13C_sd / sqrt(bl_SoG$n)
  bl_JdF_d13C_se <- bl_JdF$d13C_sd / sqrt(bl_JdF$n)
  bl_SoG_d15N_se <- bl_SoG$d15N_sd / sqrt(bl_SoG$n2)
  bl_JdF_d15N_se <- bl_JdF$d15N_sd / sqrt(bl_JdF$n2)
  
  # Shifts and their SEs
  d13C_shift    <- bl_JdF$d13C_mean - bl_SoG$d13C_mean
  d15N_shift    <- bl_JdF$d15N_mean - bl_SoG$d15N_mean
  d13C_shift_se <- sqrt(bl_JdF_d13C_se^2 + bl_SoG_d13C_se^2)
  d15N_shift_se <- sqrt(bl_JdF_d15N_se^2 + bl_SoG_d15N_se^2)
  
  # Adjust prey means
  d13C_mean_adj <- prey$d13C_mean + d13C_shift
  d15N_mean_adj <- prey$d15N_mean + d15N_shift
  
  # Propagate uncertainty for source SD to use in MixSIAR:
  # Use prey SD (process variation) PLUS baseline shift SE (mean-alignment uncertainty)
  d13C_sd_adj <- sqrt(prey$d13C_sd^2 + d13C_shift_se^2)
  d15N_sd_adj <- sqrt(prey$d15N_sd^2 + d15N_shift_se^2)
  
  data.frame(
    d13C_mean_adj = d13C_mean_adj,
    d13C_sd_adj   = d13C_sd_adj,
    d15N_mean_adj = d15N_mean_adj,
    d15N_sd_adj   = d15N_sd_adj
  )
}

krill_adj   <- mk_adj(krill_SoG23,   bl_SoG23, bl_JdF24)
herring_adj <- mk_adj(herring_SoG23, bl_SoG23, bl_JdF24)
juv_herring_adj <- mk_adj(juv_herring_SoG23, bl_SoG23, bl_JdF24)

krill_adj # Krill adjusted to JdF24
herring_adj # Adult herring adjusted to JdF24
juv_herring_adj # Juvenile herring adjusted to JdF24

# Generate a Bayesian stable isotope mixing model (SIMM) using LE d13C and NLE d15N source and consumer data
mix <- matrix(c(
  -18.35, -17.78, -17.87, -17.93, -17.57, -16.43, -18.39, 
  -19.11, -16.47, -16.86, -17.74, -19.22, -17.37, 10.8, 10.5, 11.5,
  11.1, 11.6, 10.6, 10.4, 11.5, 11.5, 11.4, 14.3, 12.2, 13
), ncol = 2, nrow = 13) # isotope values of humpback whale fecal samples
colnames(mix) <- c("d13C", "d15N")

# The 'source' is our prey types 
source_names <- c("Krill", "Juvenile herring", "Adult herring")
source_means <- matrix(c(-17.755, -16.395, -14.755, 8.945, 12.185, 13.105), ncol =2, nrow = 3)
source_sds <- matrix(c(1.37188, 1.146278, 1.077986, 1.354668, 0.3852617, 0.464679), ncol = 2, nrow = 3)

# The 'correction' is the trophic enrichment factor between prey and feces.
correction_means <- matrix(c(-.9,-.9,-.9,2,2,2), ncol = 2, nrow = 3)
correction_sds <- matrix(c(0.96, 0.96, 0.96, 0.89, 0.89, 0.89), ncol = 2, nrow = 3)

# The isotopic composition of humpback whale feces is expected to reflect the isotopic signature 
# of recently ingested prey, modified by digestive fractionation (Sponheimer et al. 2003; Browning et al. 2014; Hückstädt et al. 2012). 
# Empirical studies of marine mammals show fecal δ¹³C values are typically 0.5–2 ‰ lower and fecal δ¹⁵N values 1–3 ‰ higher than prey, 
# consistent with isotopic discrimination during digestion.

simmr_in <- simmr_load(
  mixtures = mix,
  source_names = source_names,
  source_means = source_means,
  source_sds = source_sds,
  correction_means = correction_means,
  correction_sds = correction_sds
)

# Create an isospace plot with prey and fecal data

p <- plot(simmr_in,
          xlab = expression(paste(delta^13, "C (\u2030)", sep = "")),
          ylab = expression(paste(delta^15, "N (\u2030)", sep = "")),
          title = "",
          mix_name = "Fecal samples"
)

final_p <- p +
  theme(
    text = element_text(family = "Arial", size = 20),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", colour = "white", linewidth = 0.3),
    legend.title = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_colour_manual(
    values = c(
      "Fecal samples" = "black",
      "Krill"                = "#e41a1c",
      "Juvenile herring"     = "#377eb8",
      "Adult herring"        = "#4daf4a"
    )
  ) +
  scale_shape_manual(
    values = c(
      "Fecal samples" = 1,
      "Krill"                = 2,
      "Juvenile herring"     = 0,
      "Adult herring"        = 5
    )
  )

final_p
final_p <- final_p +
  theme(text = element_text(family = "sans", size = 16))  
final_p <- final_p +
  theme(
    axis.title.y = element_text(margin = margin(r = 8)),  
    plot.margin = margin(t = 15, r = 15, b = 15, l = 15, unit = "pt")
  )
for (i in seq_along(final_p$layers)) {
  geom_class <- class(final_p$layers[[i]]$geom)[1]
  if (geom_class %in% c("GeomErrorbarh", "GeomErrorbar")) {
    final_p$layers[[i]]$aes_params$width <- 0
    final_p$layers[[i]]$geom_params$width <- 0
  }
}
final_p
final_p <- final_p +
  scale_x_continuous(breaks = pretty_breaks(n = 4)) +
  scale_y_continuous(breaks = pretty_breaks(n = 4))

final_p

ggsave(
  filename = "isospace_plot1.pdf",
  plot = final_p,
  device = "pdf",
  width = 8,
  height = 6,
  units = "in",
  encoding = "WinAnsi.enc"
)

# Visualize prey contributions to fecal mixture (i.e., humpback whale diet)

simmr_out <- simmr_mcmc(
  simmr_in,
  mcmc_control = list(
    iter   = 1000000,
    burn   = 500000,
    thin   = 500,
    n.chain = 4
  )
)

simmr_out2 <- simmr_out$output[[1]]$BUGSoutput$sims.list$p
colnames(simmr_out2) <- simmr_out$input$source_names

df <- reshape2::melt(simmr_out2)
colnames(df) <- c("Num", "Source", "Proportion")

df$Source <- factor(df$Source, levels = c("Krill", "Juvenile herring", "Adult herring"))

label_df <- df %>%
  group_by(Source) %>%
  summarise(
    label_pos = boxplot.stats(Proportion)$stats[5] + 0.0000000001,
    .groups = "drop"
  )

contributions <- ggplot(df, aes(
  x = Source,
  y = Proportion,
  fill = Source,
  colour = Source
)) +
  geom_boxplot(
    notch = TRUE,
    outlier.size = 0.5,
    linewidth = 0.45,
    alpha = 0.6
  ) +
  geom_text(
    data = label_df,
    aes(
      x = Source,
      y = label_pos,
      label = Source
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = -0.8,
    size = 5,
    colour = "black"
  ) +
  scale_fill_manual(values = c(
    "Krill" = "#e41a1c",
    "Juvenile herring" = "#377eb8",
    "Adult herring" = "#4daf4a"
  )) +
  scale_colour_manual(values = c(
    "Krill" = "#e41a1c",
    "Juvenile herring" = "#377eb8",
    "Adult herring" = "#4daf4a"
  )) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  theme_classic(base_family = "sans", base_size = 16) +
  theme(
    axis.line = element_line(linewidth = 0.4, colour = "black"),
    axis.ticks.x = element_line(linewidth = 0.4, colour = "black"),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.title.x = element_text(margin = margin(t = 8)),
    plot.margin = margin(15, 45, 15, 15)
  ) +
  xlab("Source") +
  ylab("Proportion of diet")

contributions

ggsave(
  filename = "contributions.pdf",
  plot = contributions,
  device = "pdf",
  width = 8,
  height = 6,
  units = "in",
  encoding = "WinAnsi.enc"
)

# Isotopic discrimination factor sensitivity analysis

df <- tribble(
  ~Prey, ~Factor, ~Increase50, ~Decrease50,
  "Krill", "N15", 7.1, -14.3,
  "Krill", "C13", -11.8, 4.1,
  "Juvenile herring", "N15", -4.4, 9.9,
  "Juvenile herring", "C13", 4.8, -1.1,
  "Adult herring", "N15", -2.7, 4.4,
  "Adult herring", "C13", 7,-3
)

plot_df <- df %>%
  mutate(
    Prey = factor(Prey, levels = c("Adult herring", "Juvenile herring", "Krill")),
    Factor = factor(Factor, levels = c("C13", "N15")),
    y0 = as.numeric(Prey),
    y = case_when(
      Factor == "N15" ~ y0 + 0.18,
      Factor == "C13" ~ y0 - 0.18
    ),
    ymin = y - 0.13,
    ymax = y + 0.13
  ) %>%
  pivot_longer(
    cols = c(Increase50, Decrease50),
    names_to = "Change",
    values_to = "Diff"
  ) %>%
  mutate(
    xmin = pmin(0, Diff),
    xmax = pmax(0, Diff),
    Change = recode(
      Change,
      "Increase50" = "Increase discrimination factor by 50%",
      "Decrease50" = "Decrease discrimination factor by 50%"
    )
  )

factor_labels <- tibble(
  x = rep(-21.5, 6),
  y = c(1.18, 0.82, 2.18, 1.82, 3.18, 2.82),
  label = c(
    "Delta^{15}*N[\"feces-diet\"]",
    "Delta^{13}*C[\"feces-diet\"]",
    "Delta^{15}*N[\"feces-diet\"]",
    "Delta^{13}*C[\"feces-diet\"]",
    "Delta^{15}*N[\"feces-diet\"]",
    "Delta^{13}*C[\"feces-diet\"]"
  )
)

p <- ggplot(plot_df) +
  geom_rect(
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = Change
    ),
    colour = "black",
    linewidth = 0.5
  ) +
  geom_vline(xintercept = 0, linewidth = 0.5) +
  geom_text(
    data = factor_labels,
    aes(x = x, y = y, label = label),
    parse = TRUE,
    inherit.aes = FALSE,
    hjust = 0,
    size = 4.5
  ) +
  geom_hline(
    yintercept = c(1.5, 2.5),
    linetype = "dashed",
    colour = "grey60",
    linewidth = 0.4
  ) +
  scale_y_continuous(
    breaks = 1:3,
    labels = c("Adult herring", "Juvenile herring", "Krill"),
    limits = c(0.5, 3.5),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_x_continuous(
    limits = c(-23, 20),
    breaks = seq(-20, 20, 5),
    labels = function(x) ifelse(x > 0, paste0("+", x, "%"), paste0(x, "%"))
  ) +
  scale_fill_manual(
    values = c(
      "Decrease discrimination factor by 50%" = "#F15D59",
      "Increase discrimination factor by 50%" = "#17B02A"
    )
  ) +
  labs(
    x = "Difference from base model estimate of prey contribution (%)",
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 15) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 12),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14),
    axis.title.x = element_text(size = 15),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5.5, 5.5, 5.5, 45)
  ) +
  coord_cartesian(clip = "off")

p

ggsave(
  "SIMM_tornado.pdf",
  p,
  width = 8.5,
  height = 4.5
)