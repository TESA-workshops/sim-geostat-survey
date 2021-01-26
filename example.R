library(dplyr)
library(ggplot2)
library(sdmTMB)
library(SimSurvey)
theme_set(theme_light())

set.seed(1104)
sim <- SimSurvey::sim_abundance(ages = seq(1, 10), years = seq(1, 10)) %>%
  SimSurvey::sim_distribution(
    grid = SimSurvey::make_grid(res = c(10, 10), depth_range = c(10, 500)),
    ays_covar = SimSurvey::sim_ays_covar(phi_age = 0.8, phi_year = 0.1),
    depth_par = SimSurvey::sim_parabola(mu = 200, sigma = 30)
  )
xy_sim <- tibble::as_tibble(sim$grid_xy)
df_sim <- tibble::as_tibble(sim$sp_N)
df_sim <- left_join(df_sim, xy)
df_sim_sum <- group_by(df_sim, year, x, y, depth, cell) %>%
  summarise(N = sum(N), .groups = "drop")

survey <- SimSurvey::sim_survey(sim)
xy <- tibble::as_tibble(survey$grid_xy)
df <- tibble::as_tibble(survey$sp_N)
df <- left_join(df, xy)
df_sum <- group_by(df, year, x, y, depth, cell) %>%
  summarise(N = sum(N), .groups = "drop")

ggplot(df, aes(x, y, fill = log(N + 1))) +
  geom_tile() +
  facet_grid(age ~ year) +
  scale_fill_viridis_c()

# SA: weird edge effects?
df_sum <- dplyr::filter(df_sum, x > -110, x < 110)
ggplot(df_sum, aes(x, y, fill = log(N + 1))) +
  geom_tile() +
  facet_wrap(~year) +
  scale_fill_viridis_c()

ggplot(df_sum, aes(x, y, fill = depth)) +
  geom_tile() +
  facet_wrap(~year) +
  scale_fill_viridis_c(option = "C")

unique(df_sum$year)
range(df_sum$N)
range(df_sum$x)
range(df_sum$y)
nrow(df_sum)

dat <- df_sum
mesh <- sdmTMB::make_mesh(dat, xy_cols = c("x", "y"), cutoff = 20)
mesh$spde$mesh$n
plot(mesh)

fit <- sdmTMB(N ~ 0 + as.factor(year),
  data = dat,
  family = nbinom2(link = "log"), spde = mesh,
  include_spatial = FALSE, ar1_fields = TRUE, time = "year",
  silent = FALSE
)
print(fit)

dat$resid <- residuals(fit)

ggplot(dat, aes(x, y, colour = resid, size = abs(resid))) +
  geom_point() +
  facet_wrap(~year) +
  scale_colour_gradient2() +
  scale_size_area(max_size = 2)

hist(dat$resid)
qqnorm(dat$resid)
qqline(dat$resid)

grid_dat <- tidyr::expand_grid(x = sort(unique(dat$x)), y = sort(unique(dat$y)))
grid_dat <- purrr::map_dfr(sort(unique(dat$year)), ~ bind_cols(grid_dat, year = .))
pred <- predict(fit, newdata = grid_dat, return_tmb_object = TRUE)

dplyr::filter(df_sim_sum, N > 0, x > -110, x < 110) %>%
  ggplot(aes(x, y, fill = log(N))) +
  geom_tile() +
  facet_wrap(~year) +
  scale_fill_viridis_c()

ggplot(pred$data, aes(x, y, fill = est)) +
  geom_tile() +
  facet_wrap(~year) +
  scale_fill_viridis_c()

true <- dplyr::filter(df_sim_sum, N > 0, x > -110, x < 110) %>%
  select(x, y, year, N) %>%
  mutate(type = "True")
fitted <- pred$data %>%
  select(x, y, year, N = est) %>%
  mutate(type = "Predicted") %>%
  mutate(N = exp(N))
both <- bind_rows(true, fitted)
ggplot(both, aes(x, y, fill = log(N))) +
  geom_tile() +
  facet_grid(type ~ year) +
  scale_fill_viridis_c()

# This just looks weird:
ggplot(pred$data, aes(x, y, fill = est)) +
  geom_tile() +
  facet_wrap(~year) +
  scale_fill_viridis_c() +
  geom_point(aes(x, y, size = log(N)),
    inherit.aes = FALSE, pch = 21, alpha = 0.5,
    data = dplyr::filter(df_sim_sum, N > 0, x > -110, x < 110)
  ) +
  scale_size_area()

index <- get_index(pred)

true_abund <- dplyr::filter(df_sim_sum, x > -110, x < 110) %>%
  group_by(year) %>%
  summarise(N = sum(N)) %>%
  mutate(type = "True")
both <- index %>% mutate(type = "Estimated", N = est) %>%
  bind_rows(true_abund)

ggplot(both, aes(year, N, group = type)) +
  geom_line(aes(colour = type, lty = type)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = type), alpha = 0.3) +
  labs(x = "Year", y = "Abundance", colour = "Type", fill = "Type", lty = "Type") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")
