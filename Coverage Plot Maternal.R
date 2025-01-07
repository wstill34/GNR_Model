simulate_reduced_hospitalizations_mat <- function(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, vm_se, em, em_se, min_coverage, max_coverage, step) {
  
  
  baseline_hospitalizations_mean_sd <- simulate_status_quo(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations)
  
  coverage_levels <- seq(min_coverage, max_coverage, by = step)
  reduced_hospitalizations_percent <- numeric(length(coverage_levels))
  reduced_hospitalizations_n <- numeric(length(coverage_levels))
  
  #simulate reduced_hospitalizations for each coverage level
  
  for (i in seq_along(coverage_levels)) {
    vm <- coverage_levels[i]
    reduced_hospitalizations_percent[i] <- (baseline_hospitalizations_mean_sd - simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, vm_se, em, em_se))/baseline_hospitalizations_mean_sd
    reduced_hospitalizations_n[i] <- baseline_hospitalizations_mean_sd - simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, vm_se, em, em_se)
    }
  
  return(data.frame(coverage = coverage_levels, reducedhosppercent = reduced_hospitalizations_percent, reducedhospn = reduced_hospitalizations_n))
  #return(data.frame(coverage = coverage_levels, reducedhospitalizations2 = reduced_hospitalizations_percentage))
  #return(c(reduced_hospitalizations_percentage))
  
}

#generate data
coverage_data_maternal <- simulate_reduced_hospitalizations_mat(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, v1, v1_se, e1, e1_se, min_coverage, max_coverage, step)

#plot coverage (x-axis) and reduced hospitalizations (y-axis)

library(ggplot2)

ggplot(coverage_data_maternal, aes(x= coverage, y = reducedhospn)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red") +
  labs(title = "Reduced Hospitalizations vs Coverage, Maternal Vaccine",
       x = "Coverage",
       y = "Reduced Hospitalizations") +
  theme_minimal()

ggplot(coverage_data_maternal, aes(x= coverage, y = reducedhosppercent)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red") +
  labs(title = "Reduced Hospitalizations vs Coverage, Maternal Vaccine",
       x = "Coverage",
       y = "Reduced Hospitalizations (percentage)") +
  theme_minimal()