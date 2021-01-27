library(dplyr)
library(ggplot2)
library(sdmTMB)
library(SimSurvey)
theme_set(theme_light())

set.seed(2849)
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

survey <- SimSurvey::sim_survey(sim, n_sims = 1)
plot(colSums(survey$I), type = "o")
xy <- tibble::as_tibble(survey$grid_xy)
df <- tibble::as_tibble(survey$setdet) %>% select(x, y, set, year, N = n, tow_area)
df <- left_join(df, xy)
# df_sum <- dplyr::filter(df_sum, x > -110, x < 110)

ggplot(df, aes(x, y, colour = log(N + 1), size = N)) +
  geom_point() +
  facet_wrap(vars(year)) +
  scale_colour_viridis_c() +
  scale_size_area(max_size = 4)

# SA: weird edge effects?
ggplot(df_sim_sum, aes(x, y, fill = log(N + 1))) +
  geom_tile() +
  facet_wrap(~year) +
  scale_fill_viridis_c()

unique(df$year)
range(df$N)
range(df$x)
range(df$y)
nrow(df)

dat <- df
mesh <- sdmTMB::make_mesh(dat, xy_cols = c("x", "y"), cutoff = 15)
mesh$spde$mesh$n
plot(mesh)
dat$offset <- log(dat$tow_area)

fit <- sdmTMB(N ~ 0 + as.factor(year) + offset,
  data = dat,
  family = nbinom2(link = "log"), spde = mesh,
  include_spatial = FALSE, ar1_fields = TRUE,
  time = "year",
  silent = FALSE
)
print(fit)

dat$resid <- residuals(fit)

ggplot(dat, aes(x, y, colour = resid, size = abs(resid))) +
  geom_point() +
  facet_wrap(~year) +
  scale_colour_gradient2() +
  scale_size_area(max_size = 2)

# hist(dat$resid)
# qqnorm(dat$resid)
# qqline(dat$resid)

grid_dat <- tidyr::expand_grid(x = sort(unique(df_sim_sum$x)), y = sort(unique(df_sim_sum$y)))
grid_dat <- purrr::map_dfr(sort(unique(dat$year)), ~ bind_cols(grid_dat, year = .))
grid_dat$offset <- mean(dat$offset)
pred <- predict(fit, newdata = grid_dat, return_tmb_object = TRUE)

# TODO: adjust this for catchability:
true <- df_sim_sum %>%
  # dplyr::filter(df_sim_sum, N > 0, x > -110, x < 110) %>%
  select(x, y, year, N) %>%
  mutate(type = "True")
fitted <- pred$data %>%
  select(x, y, year, N = est) %>%
  mutate(type = "Predicted") %>%
  mutate(N = exp(N))
both <- bind_rows(true, fitted) %>%
  group_by(type) %>%
  mutate(N_scaled = N / exp(mean(log(N))))
ggplot(both, aes(x, y, fill = log(N_scaled))) +
  geom_tile() +
  facet_grid(type ~ year) +
  scale_fill_viridis_c(limits = quantile(log(both$N_scaled), probs = c(0.05, 1)))

index <- get_index(pred)

true_abund <- tibble(year = unique(df$year), N = as.numeric(colSums(survey$I))) %>%
  mutate(type = "True")
both <- index %>%
  mutate(type = "Estimated", N = est) %>%
  bind_rows(true_abund) %>%
  group_by(type) %>%
  mutate(geo_mean = exp(mean(log(N), na.rm = TRUE)),
    lwr_scaled = lwr / geo_mean, N_scaled = N / geo_mean, upr_scaled = upr / geo_mean)

ggplot(both, aes(year, N_scaled, group = type)) +
  geom_line(aes(colour = type, lty = type)) +
  geom_ribbon(aes(ymin = lwr_scaled, ymax = upr_scaled, fill = type), alpha = 0.3) +
  labs(x = "Year", y = "Abundance", colour = "Type", fill = "Type", lty = "Type") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2")
