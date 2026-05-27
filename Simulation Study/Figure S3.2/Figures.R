library(dplyr)
library(ggplot2)
library(forcats)

# Read and label timing results
dat = read.csv("time_summary.csv")
dat$method = c(
  "IPW_icdf", "OR_CPM_icdf", "AIPW_CPM_icdf", "TMLE", "AIPW",
  "IPW_Firpo", "TMLE_CPM", "TMLE_cqr", "AIPW_CPM"
)

keep = c(
  "OR_CPM_icdf", "AIPW_CPM_icdf", "TMLE", "AIPW",
  "IPW_Firpo", "TMLE_CPM", "TMLE_cqr", "AIPW_CPM"
)
df = dat %>% filter(method %in% keep)

B = 1000
df = df %>%
  mutate(
    se = sd / sqrt(B),
    lo = mean - 1.96 * se,
    hi = mean + 1.96 * se,
    method = fct_reorder(method, mean, .desc = TRUE)
  )

x_pad = 0.10 * max(df$hi)

p = ggplot(df, aes(x = mean, y = method)) +
  geom_col(width = 0.65, fill = "grey70", color = "black", linewidth = 0.25) +
  geom_errorbar(
    aes(xmin = lo, xmax = hi),
    orientation = "y",
    width = 0.18,
    linewidth = 0.35
  ) +
  geom_text(aes(label = sprintf("%.2f", mean)), hjust = -0.15, size = 3.2) +
  scale_x_continuous(limits = c(0, max(df$hi) + x_pad), expand = expansion(mult = c(0, 0))) +
  labs(
    x = "Mean elapsed time (seconds)",
    y = "Estimators"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.text.y = element_text(size = 10),
    plot.caption = element_text(hjust = 0, size = 9),
    plot.margin = margin(10, 20, 8, 10)
  )

ggsave("Fig_S3_2_comp_time.pdf", p, width = 7.2, height = 4.6, units = "in")
ggsave("Fig_S3_2_comp_time.png", p, width = 7.2, height = 4.6, units = "in", dpi = 600)
