library(dplyr)
library(ggplot2)
library(SimSurvey)
library(sdmTMB)
future::plan(future::multisession) # for the furrr package; or `plan(sequential)`

sim_and_fit <- function(iter) {
  set.seed(iter * 283028)
  sim <- sim_abundance(ages = seq(1, 10), years = seq(1, 10)) %>%
    sim_distribution(
      grid = make_grid(res = c(10, 10), depth_range = c(10, 500)),
      ays_covar = sim_ays_covar(phi_age = 0.8, phi_year = 0.1),
      depth_par = sim_parabola(mu = 200, sigma = 30)
    )

  survey <- sim_survey(sim, n_sims = 1)
  xy <- as_tibble(survey$grid_xy)
  dat <- as_tibble(survey$setdet) %>%
    select(x, y, set, year, N = n, tow_area)
  dat <- left_join(dat, xy, by = c("x", "y"))
  dat$offset <- log(dat$tow_area)

  grid <- as_tibble(select(sim$grid_xy, x, y))
  grid_dat <- tidyr::expand_grid(x = sort(unique(grid$x)), y = sort(unique(grid$y)))
  grid_dat <- purrr::map_dfr(sort(unique(dat$year)), ~ bind_cols(grid_dat, year = .))
  grid_dat$offset <- mean(dat$offset)

  mesh <- sdmTMB::make_mesh(dat, xy_cols = c("x", "y"), cutoff = 20)
  fit <- sdmTMB(N ~ 0 + as.factor(year) + offset,
    data = dat,
    family = nbinom2(link = "log"), spde = mesh,
    include_spatial = TRUE, time = "year"
  )

  pred <- predict(fit, newdata = grid_dat, return_tmb_object = TRUE)
  index <- get_index(pred)

  true_abund <- tibble(year = unique(dat$year), N = as.numeric(colSums(survey$I))) %>%
    mutate(type = "True")
  mutate(index, type = "Estimated", N = est) %>%
    bind_rows(true_abund) %>%
    mutate(iter = iter)
}

result <- furrr::future_map_dfr(seq_len(16), sim_and_fit,
  .options = furrr::furrr_options(seed = TRUE))

not_converged <- dplyr::filter(result, (bad_eig | max_gradient > 0.001) & type == "Estimated")
stopifnot(nrow(not_converged) == 0L)

result_scaled <- result %>%
  group_by(type, iter) %>%
  mutate(geo_mean = exp(mean(log(N), na.rm = TRUE)),
    lwr_scaled = lwr / geo_mean, N_scaled = N / geo_mean, upr_scaled = upr / geo_mean) %>%
  ungroup()

result_scaled %>%
  ggplot(aes(year, N_scaled, group = type)) +
  geom_line(aes(colour = type, lty = type)) +
  geom_ribbon(aes(ymin = lwr_scaled, ymax = upr_scaled, fill = type), alpha = 0.3) +
  labs(x = "Year", y = "Abundance", colour = "Type", fill = "Type", lty = "Type") +
  scale_color_manual(values = c("Estimated" = "grey30", "True" = "red")) +
  scale_fill_manual(values = c("Estimated" = "grey30", "True" = "red")) +
  facet_wrap(~iter, ncol = 4, scales = "free_y")
