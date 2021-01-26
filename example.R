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

xy <- tibble::as_tibble(sim$grid_xy)
df <- tibble::as_tibble(sim$sp_N)
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
max(df_sum$N)
min(df_sum$N)
range(df_sum$x)
range(df_sum$y)

dat <- df_sum
dat <- dplyr::filter(df_sum, N > 10) %>%
  # dat <- mutate(df_sum, N = if_else(N > 10, N, 0)) %>% # For Tweedie?
  mutate(N1000 = N / 1000)

ggplot(dat, aes(x, y, fill = log(N1000))) +
  geom_tile() +
  facet_wrap(~year) +
  scale_fill_viridis_c()

range(dat$N1000)

mesh <- sdmTMB::make_mesh(dat, xy_cols = c("x", "y"), cutoff = 20)
mesh$spde$mesh$n
plot(mesh)

fit <- sdmTMB(N1000 ~ 0 + as.factor(year),
  data = dat,
  family = Gamma(link = "log"), spde = mesh,
  include_spatial = FALSE, ar1_fields = TRUE, time = "year",
  silent = FALSE
)
print(fit)

dat$resid <- residuals(fit)
hist(dat$resid)
qqnorm(dat$resid)
qqline(dat$resid)

grid_dat <- tidyr::expand_grid(x = sort(unique(dat$x)), y = sort(unique(dat$y)))
grid_dat <- purrr::map_dfr(sort(unique(dat$year)), ~ bind_cols(grid_dat, year = .))
pred <- predict(fit, newdata = grid_dat, return_tmb_object = TRUE)

ggplot(dat, aes(x, y, fill = log(N1000))) +
  geom_tile() +
  facet_wrap(~year) +
  scale_fill_viridis_c()

ggplot(pred$data, aes(x, y, fill = est)) +
  geom_tile() +
  facet_wrap(~year) +
  scale_fill_viridis_c()

# This just looks weird:
ggplot(pred$data, aes(x, y, fill = est)) +
  geom_tile() +
  facet_wrap(~year) +
  scale_fill_viridis_c() +
  geom_point(aes(x, y, size = log(N1000)),
    inherit.aes = FALSE, pch = 21, alpha = 0.5, data = dat
  ) +
  scale_size_area()

index <- get_index(pred)

ggplot(index, aes(year, est)) +
  geom_line() +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  labs(x = "Year", y = "Index")
