# Creates function to simulate the reduction in hospitalizations for Childhood vaccine based on varying coverage levels
simulate_reduced_hospitalizations_epi <- function(birth_cohort, periods, h, h_vec_sim, mu_gnr_sim, mu_ac, mu_gnr, n_simulations, v1, efficacy_epi, min_coverage, max_coverage, step) {
  
  
  baseline_hospitalizations <- status_quo_data[1,1]

  coverage_levels <- seq(min_coverage, max_coverage, by = step)
  reduced_hospitalizations_percent <- numeric(length(coverage_levels))
  reduced_hospitalizations_n <- numeric(length(coverage_levels))
  
  #simulate reduced_hospitalizations for each coverage level
  
  for (i in seq_along(coverage_levels)) {
     v1 <- coverage_levels[i]
     reduced_hospitalizations_percent[i] <- (baseline_hospitalizations - simulate_epi(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, v1, efficacy_epi)[1,1])/baseline_hospitalizations
     reduced_hospitalizations_n[i] <- baseline_hospitalizations - simulate_epi(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, v1, efficacy_epi)[1,1]
     
  }
  
 return(data.frame(coverage = coverage_levels, reducedhosppercent = reduced_hospitalizations_percent, reducedhospn = reduced_hospitalizations_n))
  #return(data.frame(coverage = coverage_levels, reducedhospitalizations2 = reduced_hospitalizations_percentage))
  #return(c(reduced_hospitalizations_percentage))
  
}

#generate data
coverage_data_epi <- simulate_reduced_hospitalizations_epi(birth_cohort, periods, h, h_vec_sim, mu_gnr_sim, mu_ac, mu_gnr, n_simulations, v1, efficacy_epi, min_coverage, max_coverage, step)
#simulate_reduced_hospitalizations_epi(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, v1, v1_se, e1, e1_se, min_coverage, max_coverage, step)


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
