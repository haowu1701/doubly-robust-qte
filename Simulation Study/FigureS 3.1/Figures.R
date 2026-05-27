source("Helper_functions.R")

set.seed(37023)
library(ggplot2)

n = 1000
alpha1 = 0.5
alpha2 = 0.35
beta1 = -2
beta2 = 3
delta = 2

data = generate_data(n, alpha1, alpha2, beta1, beta2, delta)

Y = data$Y
logY = log(data$Y)
cY = qnorm(pchisq(data$Y, df = 5))

df_plot = rbind(
  data.frame(Value = Y,    Scale = "Truly observed outcome"),
  data.frame(Value = logY, Scale = "Log transformation"),
  data.frame(Value = cY,   Scale = "Correct transformation")
)

df_plot$Scale = factor(
  df_plot$Scale,
  levels = c(
    "Truly observed outcome",
    "Log transformation",
    "Correct transformation"
  )
)

p = ggplot(df_plot, aes(x = Value)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 40,
    fill = "grey75",
    color = "black"
  ) +
  facet_wrap(~Scale, scales = "free_x", nrow = 1) +
  labs(
    x = "Outcome value",
    y = "Density"
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.text = element_text(size = 14),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12)
  )

p

ggsave(
  "FigureS3.1.pdf",
  plot = p,
  width = 9,
  height = 3,
  dpi = 300
)
