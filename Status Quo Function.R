
simulate_status_quo <- function (birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations) {
  

# Creates matrices to store the results from simulations
hospitalizations_sim <- matrix(0, nrow=n_simulations, ncol=periods)
#Added h_vec_sim matrix of 1000 rows, 4 columns
h_vec_sim <- matrix(0, nrow=n_simulations, ncol=periods)
#Also added mu_gnr_sim matrix of 1000 rows, 4 columns
mu_gnr_sim <- matrix(0, nrow=n_simulations, ncol=periods)
deaths_due_to_gnr_sim <- matrix(0, nrow=n_simulations, ncol=periods)
other_cause_deaths_sim <- matrix(0, nrow=n_simulations, ncol=periods)
total_deaths_sim <- matrix(0, nrow=n_simulations, ncol=periods)
survived_sim <- matrix(0, nrow=n_simulations, ncol=periods)


# Monte Carlo simulation
for (sim in 1:n_simulations) {
  remaining_cohort <- birth_cohort
 
  for (age in 1:periods) {
   # Calculates number of hospitalizations, deaths due to GNR, and other cause deaths
    h_vec_sim[sim, age] <- rbeta(1, h_vec[age], h_n[age]-h_vec[age]) #For each age group, estimates a hospitalization probability based on the mean and standard deviation and 1 trial
    hospitalizations_sim[sim, age] <- rbinom(1, remaining_cohort, h_vec_sim[sim,age]) #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
    mu_gnr_sim[sim, age] <- rbeta(1, mu_gnr_vec[age], mu_gnr_n[age] - mu_gnr_vec[age])
    deaths_due_to_gnr_sim[sim, age] <- rbinom(1, hospitalizations_sim[sim, age], mu_gnr_sim[sim,age])
    other_cause_deaths_sim[sim, age] <- rbinom(1, remaining_cohort, mu_ac[age]) - deaths_due_to_gnr_sim[sim, age]
     
   # Calculates the number of survivors
    survived_sim[sim, age] <- remaining_cohort - deaths_due_to_gnr_sim[sim, age] - other_cause_deaths_sim[sim, age]
    total_deaths_sim[sim, age] <- deaths_due_to_gnr_sim[sim, age] + other_cause_deaths_sim[sim, age]
   
    # Update the remaining cohort for the next period
    remaining_cohort <- survived_sim[sim, age]
  }
}


# Summarizes the results with confidence intervals
# First add the total hospitalizations and deaths for each simulation
total_hospitalizations_sim <- rowSums(hospitalizations_sim)
total_deaths_due_to_gnr_sim <- rowSums(deaths_due_to_gnr_sim)
total_other_cause_deaths_sim <- rowSums(other_cause_deaths_sim)
total_all_cause_deaths_sim <- rowSums(total_deaths_sim)
expected_hospitalizations <- mean(total_hospitalizations_sim)
std_dev_hospitalizations <- sd(total_hospitalizations_sim)
expected_gnr_deaths <- mean(total_deaths_due_to_gnr_sim)
std_dev_gnr_deaths <- sd(total_deaths_due_to_gnr_sim)

#return (c(expected_hospitalizations, std_dev_hospitalizations))
#return (expected_hospitalizations)

return(data.frame(expected_hosp = expected_hospitalizations, std_dev_hosp = std_dev_hospitalizations, expected_gnr_deaths = expected_gnr_deaths, std_dev_gnr_deaths = std_dev_gnr_deaths))
#return(data.frame(coverage = coverage_levels, reducedhospitalizations2 = reduced_hospitalizations_percentage))
#return(c(reduced_hospitalizations_percentage))

}

#generate data
status_quo_data <- simulate_status_quo(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations)