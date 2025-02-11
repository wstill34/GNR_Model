# Creates function to simulate the reduction in hospitalizations for Childhood vaccine based on varying coverage levels
simulate_reduced_hospitalizations_epi <- function(birth_cohort, periods, h, h_vec_sim, mu_gnr_sim, mu_ac, mu_gnr, n_simulations, v1, efficacy_epi, min_coverage, max_coverage, step) {
  
  
  baseline_hospitalizations <- status_quo_data[,1]
  baseline_deaths <- status_quo_data[,16]

  coverage_levels <- seq(min_coverage, max_coverage, by = step)
  reduced_hospitalizations_percent <- numeric(length(coverage_levels))
  reduced_hospitalizations_n <- numeric(length(coverage_levels))
  reduced_deaths_percent <- numeric(length(coverage_levels))
  reduced_deaths_n <- numeric(length(coverage_levels))
  
  #simulate reduced_hospitalizations for each coverage level
  
  for (i in seq_along(coverage_levels)) {
     v1 <- coverage_levels[i]
     reduced_hospitalizations_percent[i] <- (baseline_hospitalizations - simulate_epi(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, v1, efficacy_epi)[,1])/baseline_hospitalizations
     reduced_hospitalizations_n[i] <- baseline_hospitalizations - simulate_epi(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, v1, efficacy_epi)[,1]
     reduced_deaths_percent[i] <- (baseline_deaths - simulate_epi(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, v1, efficacy_epi)[,16])/baseline_deaths
     reduced_deaths_n[i] <- baseline_deaths - simulate_epi(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, v1, efficacy_epi)[,16]
  
  }
  
 return(data.frame(coverage = coverage_levels, reduced_hosp_percent = reduced_hospitalizations_percent, reduced_hosp_n = reduced_hospitalizations_n,
                   reduced_deaths_percent = reduced_deaths_percent, reduced_deaths_n = reduced_deaths_n))
}

#generate data
coverage_data_epi <- simulate_reduced_hospitalizations_epi(birth_cohort, periods, h, h_vec_sim, mu_gnr_sim, mu_ac, mu_gnr, n_simulations, v1, efficacy_epi, min_coverage, max_coverage, step)
#simulate_reduced_hospitalizations_epi(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, v1, v1_se, e1, e1_se, min_coverage, max_coverage, step)

library("writexl")
write_xlsx(coverage_data_epi,"C:/Users/wstil/OneDrive/Desktop/Aim 3/Coverage Figures/epi_coverage.xlsx", col_name=TRUE, format_headers=FALSE)


#plot coverage (x-axis) and reduced hospitalizations (y-axis)

library(ggplot2)

ggplot(coverage_data_epi, aes(x= coverage, y = reducedhospn)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red") +
  labs(title = "Reduced Hospitalizations vs Coverage, Childhood Vaccine",
       x = "Coverage",
       y = "Reduced Hospitalizations") +
  theme_minimal()
ggsave("EPI Coverage Graph.png", width = 5, height = 5)

ggplot(coverage_data_epi, aes(x= coverage, y = reducedhosppercent)) +
  geom_line(color = "blue", size = 1) +
  geom_point(color = "red") +
  labs(title = "Reduced Hospitalizations vs Coverage, Childhood Vaccine",
       x = "Coverage",
       y = "Reduced Hospitalizations (percentage)") +
  theme_minimal()
