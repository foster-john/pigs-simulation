library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(purrr)
library(mgcv)
library(lme4)

analysis_dir <- "analysis"
model_dir <- "betaSurvival_uniqueAreaTrapSnare"
path <- file.path(analysis_dir, model_dir)
density_dirs <- list.files(path)


density_dirs <- grep("density_", density_dirs, value = TRUE)

map_files2 <- function(dirs_vec, file_name){

  get_files <- function(density_dir, file_name, node){
    sim_results <- file.path(path, density_dir)
    ls <- read_rds(file.path(sim_results, file_name))
  }

  dirs_vec |>
    map(\(x) get_files(x, file_name)) |>
    list_rbind() |>
    mutate(start_density = as.factor(start_density))

}

message("Get abundance summaries")
f_name <- "abundance_summaries.rds"
df <- map_files2(density_dirs, f_name)

# create property ID for easier joining
property_ids <- df |>
  select(start_density, simulation, property) |>
  distinct() |>
  mutate(property_id = paste(start_density, simulation, property, sep = "-"))

density <- left_join(df, property_ids) |>
  select(property_id, PPNum, low_density, med_density, high_density, density, obs_flag)

message("Get take summaries")
f_name <- "take_summaries.rds"
df <- map_files2(density_dirs, f_name) |>
  left_join(property_ids)

n_methods_pp <- df |>
  select(property_id, PPNum, method) |>
  distinct() |>
  group_by(property_id, PPNum) |>
  count() |>
  rename(n_methods_used = n)

n_return <- df |>
  left_join(property_ids) |>
  select(property_id, PPNum, method) |>
  distinct() |>
  pivot_wider(names_from = method,
              values_from = method) |>
  unite(method,
        -c(property_id, PPNum),
        sep = ", ",
        na.rm = TRUE) |>
  group_by(property_id, method) |>
  mutate(return_interval = c(0, diff(PPNum))) |>
  ungroup() |>
  rename(methods_used = method) |>
  left_join(n_methods_pp)

take_joint_return <- left_join(df, n_return)

take <- take_joint_return |>
  select(property_id, PPNum, sum_take, property_area, methods_used) |>
  mutate(sum_take = sum_take / property_area) |>
  unique()

density_take <- left_join(density, take)

max_change_properties <- density_take |>
  group_by(property_id) |>
  filter(PPNum == min(PPNum) | PPNum == max(PPNum)) |>
  summarise(delta_density = diff(density),
            delta_time = diff(PPNum)) |>
  filter(delta_density == max(delta_density) |
           delta_density == min(delta_density) |
           delta_time == max(delta_time),
         abs(delta_density) > 8) |>
  pull(property_id)

save_gg <- function(dest, gg, path){
  if(!dir.exists(path)) dir.create(path, showWarnings = FALSE, recursive = TRUE)
  ggsave(
    filename = dest,
    plot = gg,
    device = "jpeg",
    path = path,
    width = 7,
    height = 5,
    units = "in",
    dpi = "retina",
    bg = "white"
  )
}

extant_properties <- density_take |>
  group_by(property_id) |>
  filter(PPNum == max(PPNum),
         density > 0) |>
  pull(property_id)

max_extant_properties <- density_take |>
  filter(property_id %in% extant_properties) |>
  group_by(property_id) |>
  filter(PPNum == min(PPNum) | PPNum == max(PPNum)) |>
  summarise(delta_density = diff(density),
            delta_time = diff(PPNum)) |>
  filter(delta_time == max(delta_time),
         abs(delta_density) > 2) |>
  pull(property_id)


plot_timeseries <- function(df, property){
  df |>
    filter(property_id == property) |>
    mutate(obs_type = if_else(obs_flag == 1, "Effort", "No effort"),
           time = 1:n()) |>
    ggplot() +
    aes(x = time, y = med_density, ymin = low_density, ymax = high_density) +
    geom_ribbon(aes(fill = "95% CI"), alpha = 0.6) +
    geom_line(aes(linetype = "Median")) +
    geom_point(aes(y = density, color = obs_type), size = 2) +
    geom_point(aes(y = sum_take, shape = methods_used), size = 2) +
    scale_fill_manual(values = "gray") +
    scale_shape_manual(na.translate = FALSE, values = 1:14) +
    labs(x = "Time",
         y = "Density (pigs / km sq.)",
         color = "Known density",
         shape = "Method(s) used",
         linetype = element_blank(),
         shape = element_blank(),
         fill = element_blank()) +
    theme_bw() +
    theme(axis.text = element_text(size = 13),
          axis.title = element_text(size = 16),
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 13))
}

properties_2_plot <- unique(c(max_change_properties, max_extant_properties))

for(i in seq_along(properties_2_plot)){
  gg <- density_take |>
    plot_timeseries(properties_2_plot[i])
  print(gg)

  save_gg(paste0("timeSeries_", properties_2_plot[i], ".jpeg"),
          gg,
          "Plots/simulationTimeSeries")

}

message("Get abundance summaries")
f_name <- "abundance_summaries.rds"
df <- map_files2(density_dirs, f_name)


xx <- lm(abundance ~ med_abundance, data = df)
xxs <- summary(xx)
rsq <- round(xxs$adj.r.squared, 2)
intercept <- xxs$coefficients[1,1]
slope <- xxs$coefficients[2,1]


gg <- df |>
  ggplot(aes(x = abundance, y = med_abundance)) +
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1, color = "1:1")) +
  geom_abline(aes(intercept = intercept, slope = slope, color = "Best fit")) +
  labs(x = "True abundance",
       y = "Median abundance",
       title = "Model performance",
       subtitle = paste0("R. sq.: ", rsq),
       color = element_blank()) +
  theme_bw() +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        title = element_text(size = 16),
        # subtitle = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13))

save_gg("observedPredictedAbundance.jpeg",
        gg,
        "Plots")



f_name <- "take_summaries.rds"
df <- map_files2(density_dirs, f_name) |>
  left_join(property_ids)

xx <- lm(take ~ med, data = df)
xxs <- summary(xx)
rsq <- round(xxs$adj.r.squared, 2)
intercept <- xxs$coefficients[1,1]
slope <- xxs$coefficients[2,1]


gg <- df |>
  ggplot(aes(x = take, y = med)) +
  geom_point() +
  geom_abline(aes(intercept = 0, slope = 1, color = "1:1")) +
  geom_smooth(method = "lm") +
  facet_wrap(~ method) +
  labs(x = "True take",
       y = "Median take",
       title = "Model performance",
       color = element_blank()) +
  theme_bw() +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size = 16),
        title = element_text(size = 16),
        # subtitle = element_text(size = 13),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 13))

save_gg("observedPredictedTake.jpeg",
        gg,
        "Plots")




