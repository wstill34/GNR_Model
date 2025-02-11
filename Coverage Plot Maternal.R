# Creates function to simulate the reduction in hospitalizations for maternal vaccine based on varying coverage levels
simulate_reduced_hospitalizations_mat <- function(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, efficacy_mat, min_coverage, max_coverage, step, remaining_years) {
  
  
  baseline_hospitalizations <- status_quo_data[,1]
  baseline_deaths <- status_quo_data[,16]
  
  coverage_levels <- seq(min_coverage, max_coverage, by = step)
  reduced_hospitalizations_percent <- numeric(length(coverage_levels))
  reduced_hospitalizations_n <- numeric(length(coverage_levels))
  reduced_deaths_percent <- numeric(length(coverage_levels))
  reduced_deaths_n <- numeric(length(coverage_levels))
  
  #simulate reduced_hospitalizations for each coverage level
  
  for (i in seq_along(coverage_levels)) {
    vm <- coverage_levels[i]
    reduced_hospitalizations_percent[i] <- (baseline_hospitalizations - simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, efficacy_mat, remaining_years)[,1])/baseline_hospitalizations
    reduced_hospitalizations_n[i] <- baseline_hospitalizations - simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, efficacy_mat, remaining_years)[,1]
    reduced_deaths_percent[i] <- (baseline_deaths - simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, efficacy_mat, remaining_years)[,16])/baseline_deaths
    reduced_deaths_n[i] <- baseline_deaths - simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, efficacy_mat, remaining_years)[,16]
    }
  
  return(data.frame(coverage = coverage_levels, reduced_hosp_percent = reduced_hospitalizations_percent, reduced_hosp_n = reduced_hospitalizations_n,
                    reduced_deaths_percent = reduced_deaths_percent, reduced_deaths_n = reduced_deaths_n))
  
}

#generate data
coverage_data_maternal <- simulate_reduced_hospitalizations_mat(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, efficacy_mat, min_coverage, max_coverage, step, remaining_years)

library("writexl")
write_xlsx(coverage_data_maternal,"C:/Users/wstil/OneDrive/Desktop/Aim 3/Coverage Figures/maternal_coverage.xlsx", col_name=TRUE, format_headers=FALSE)


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