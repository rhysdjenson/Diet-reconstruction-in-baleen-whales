library(ggplot2)      
library(dplyr)        
library(tidyr)        
library(readxl)      
library(scales)       

# Pairwise comparison of SIMM versus WGS of humpback whale fecal samples

df <- read_excel("pairwisefecalsamples.xlsx")

qqnorm(df$SIAkrill)
qqline(df$SIAkrill, col = "red")
shapiro.test(df$SIAkrill)
mean(df$SIAkrill)
sd(df$SIAkrill)
median(df$SIAkrill)
IQR(df$SIAkrill)
quantile(df$SIAkrill, probs = c(0.25, 0.75))

qqnorm(df$WGSkrill)
qqline(df$WGSkrill, col = "red")
shapiro.test(df$WGSkrill)
mean(df$WGSkrill)
sd(df$WGSkrill)
median(df$WGSkrill)
IQR(df$WGSkrill)
quantile(df$WGSkrill, probs = c(0.25, 0.75))

diet_long <- df %>%
  select(ID, SIAkrill, SIAfish, WGSkrill, WGSfish) %>%
  pivot_longer(
    cols = c(SIAkrill, SIAfish, WGSkrill, WGSfish),
    names_to = "variable",
    values_to = "proportion"
  ) %>%
  mutate(
    Method = case_when(
      str_detect(variable, "^SIA") ~ "Stable Isotope Analysis",
      str_detect(variable, "^WGS") ~ "Whole Genome Sequencing"
    ),
    Prey = case_when(
      str_detect(variable, "krill") ~ "Krill",
      str_detect(variable, "fish") ~ "Fish"
    ),
    Method = factor(Method, levels = c("Stable Isotope Analysis", "Whole Genome Sequencing")),
    Prey = factor(Prey, levels = c("Krill", "Fish")),
    Method_abbr = case_when(
      Method == "Stable Isotope Analysis" ~ "SIA",
      Method == "Whole Genome Sequencing" ~ "WGS"
    ),
    x_group = factor(
      paste(Prey, Method_abbr, sep = "_"),
      levels = c("Krill_SIA", "Krill_WGS", "Fish_SIA", "Fish_WGS")
    )
  ) %>%
  drop_na(proportion)

p_krill <- wilcox.test(df$SIAkrill, df$WGSkrill, paired = TRUE, exact = TRUE)$p.value
p_fish  <- wilcox.test(df$SIAfish,  df$WGSfish,  paired = TRUE, exact = TRUE)$p.value

pairwise <- ggplot(diet_long, aes(x = x_group, y = proportion, fill = Prey)) +
  geom_violin(
    trim = FALSE, alpha = 0.55, colour = "black",
    linewidth = 0.7, width = 0.9
  ) +
  
  # Median, IQR, and 1.5x IQR whiskers
  geom_boxplot(
    width = 0.15,
    fill = "white",
    colour = "black",
    linewidth = 0.7,
    outlier.shape = NA,
    coef = 1.5,
    show.legend = FALSE
  ) +
  
  # p-values, text only, no bracket
  annotate(
    "text",
    x = 1.5, y = 0.98,
    label = paste0("p = ", formatC(p_krill, digits = 2, format = "g")),
    size = 7, vjust = 0
  ) +
  annotate(
    "text",
    x = 3.5, y = 0.98,
    label = paste0("p = ", formatC(p_fish, digits = 2, format = "g")),
    size = 7, vjust = 0
  ) +
  
  scale_x_discrete(labels = c("Krill_SIA" = "SIMM", "Krill_WGS" = "WGS", "Fish_SIA" = "SIMM", "Fish_WGS" = "WGS")) +
  scale_fill_manual(name = "Prey", values = c("Krill" = "#e08a8a", "Fish" = "#6ca0a3")) +
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.25),
    labels = scales::number_format(accuracy = 0.01),
    limits = c(0, 1.08),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(x = "Method", y = "Proportion of diet") +
  labs(x = "Method", y = "Proportion of diet") +
  theme_classic(base_size = 15) +
  theme(
    axis.title.x = element_text(face = "plain", size = 20, margin = margin(r = 10)),
    axis.text.x = element_text(face = "plain", size = 20, colour = "black"),
    axis.ticks.x = element_line(linewidth = 0.7, colour = "black"),
    axis.title.y = element_text(face = "plain", size = 20, margin = margin(r = 10)),
    axis.text.y = element_text(face = "plain", size = 16, colour = "black"),
    axis.line = element_line(linewidth = 0.7, colour = "black"),
    axis.ticks.y = element_line(linewidth = 0.7, colour = "black"),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.background = element_rect(fill = "white", colour = NA)
  )

pairwise
pairwise <- pairwise +
  theme(text = element_text(family = "sans", size = 16))
ggsave(
  filename = "Pairwise.pdf",
  plot = pairwise,
  device = "pdf",
  width = 12,
  height = 8,
  units = "in",
  encoding = "WinAnsi.enc"
)