source('WindGithub.R')

library(grid)

# Trace plots
pdf('plots/TracePlotMu14.pdf')
par(mfrow = c(4, 2))
TraceParamsMu1 <- c(paste0('mu[', 1:4, ', 1]'), paste0('mu[', 1:4, ', 2]'))
traceplot(WindSamples[, TraceParamsMu1], col = c('#fde725', '#440154', '#21918c'))
dev.off()
pdf('plots/TracePlotMu56.pdf')
par(mfrow = c(4, 2))
TraceParamsMu1 <- c(paste0('mu[', 5:6, ', 1]'), paste0('mu[', 5:6, ', 2]'))
traceplot(WindSamples[, TraceParamsMu1], col = c('#fde725', '#440154', '#21918c'))
dev.off()

# Speed ppc
plot1 <- ggplot() +
  stat_density(
    data = df,
    geom = 'line',
    position = 'identity',
    aes(x = Wsp, colour = 'Observed', linetype = 'Observed'),
    linewidth = 1.2
  ) +
  stat_density(
    data = PredData,
    geom = 'line',
    position = 'identity',
    aes(x = WSp, colour = 'Predicted', linetype = 'Predicted'),
    linewidth = 1,
    ) +
  labs(
    x = 'Wind speed (km/h)',
    y = 'Density',
    colour = '',
    linetype = ''
    ) +
  scale_colour_manual(
    values = c(
      'Observed' = 'black',
      'Predicted' = 'blue'
      )
    ) +
  scale_linetype_manual(
    values = c(
      'Observed' = 'solid',
      'Predicted' = 'dashed'
    )
  ) +
  theme_minimal() +
  theme(
    legend.key.width = unit(0.9, 'cm'),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14)
    )
ggsave('plots/WindSpeedPpc.pdf', plot1, width = 10, height = 6)  

# Direction ppc
plot2 <- ggplot() +
  stat_density(
    data = df,
    geom = 'line',
    position = 'identity',
    aes(x = WdrRad, colour = 'Observed', linetype = 'Observed'),
    linewidth = 1.2
  ) +
  stat_density(
    data = PredData,
    geom = 'line',
    position = 'identity',
    aes(x = WdrRad, colour = 'Predicted', linetype = 'Predicted'),
    linewidth = 1
  ) +
  scale_x_continuous(
    breaks = seq(0, 2 * pi, pi / 4),
    labels = c('0°', '45°', '90°', '135°', '180°', '225°', '270°', '315°', '360°')
  ) +
  labs(
    x = 'Wind direction',
    y = 'Density',
    colour = '',
    linetype = ''
    ) +
  scale_colour_manual(
    values = c(
      'Observed' = 'black',
      'Predicted' = 'blue'
    )
    ) +
  scale_linetype_manual(
    values = c(
      'Observed' = 'solid',
      'Predicted' = 'dashed'
    )
  ) +
  theme_minimal() +
  theme(
    legend.key.width = unit(0.9, 'cm'),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
ggsave('plots/WindDirectionPpc.pdf', plot2, width = 10, height = 6)  

# Joint structure
plot3 <- ggplot() +
  geom_point(
    data = df,
    aes(x = u, y = v),
    alpha = 0.2,
    color = 'black',
    size = 0.8
  ) +
  geom_point(
    data = PredData,
    aes(x = u, y = v),
    alpha = 0.15,
    color = 'blue',
    size = 0.5
  ) +
  geom_segment(
    data = ClusterSummary %>% 
      filter(weight > 0.01) %>% 
      mutate(MeanDirDeg = factor(round(MeanDirDeg, 2))),
    aes(
      x = 0,
      y = 0,
      xend = MuU,
      yend = MuV,
      colour = MeanDirDeg
    ),
    arrow = arrow(length = unit(0.4, 'cm')),
    linewidth = 1.5,
    alpha = 1
  ) +
  scale_color_viridis_d(option = 'viridis') +
  coord_fixed() +
  labs(
    x = 'u component (E-W km/h)',
    y = 'v component (N-S km/h)',
    color = 'Direction (°)'
  ) +
  theme_minimal() +
  theme(
    legend.key.width = unit(0.9, 'cm'),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )
ggsave('plots/WindVectorsPpc.pdf', plot3, width = 10, height = 6)  

# Wind rose
PredDataRose <- PredData %>% 
  mutate(
    DirBin = cut((WdrDeg + 11.25) %% 360,
                 breaks = seq(0, 360, by = 22.5),
                 include.lowest = TRUE,
                 labels = c('N', 'NNE', 'NE', 'ENE', 'E', 'ESE', 'SE', 'SSE',
                            'S', 'SSW', 'SW','WSW', 'W', 'WNW', 'NW', 'NNW')
                 ),
    SpdBin = cut(WSp, breaks = c(0, 5, 10, 15, 20, 30, 100),
                 labels = c('0-5', '5-10', '10-15', '15-20', '20-30', '30+'))
  ) %>% 
  filter(!is.na(DirBin) & !is.na(SpdBin))

WndRoseSummary <- PredDataRose %>% 
  group_by(DirBin, SpdBin) %>% 
  summarise(count = n(), .groups = 'drop')

plot4 <- ggplot(
  WndRoseSummary,
  aes(
    x = DirBin,
    y = count,
    fill = SpdBin
  )
) +
  geom_col(position = 'stack') +
  coord_polar(start = -pi / 16) +
  scale_fill_viridis_d(option = 'viridis') +
  labs(
    x = '',
    y = 'Count',
    fill = 'Speed (km/h)'
  ) +
  theme_minimal() +
  theme(
    legend.key.width = unit(0.9, 'cm'),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )
ggsave('plots/WindRosePpc.pdf', plot4, width = 10, height = 6)

