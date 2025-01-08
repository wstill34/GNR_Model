#Creates a function to simulate the mean and standard deviation for number of hospitalizations and deaths under maternal vaccine scenario
simulate_maternal <- function (birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, efficacy_mat) {

# Creates matrices to store the results from simulations
hospitalizations_total_sim_mat <- matrix(0, nrow=n_simulations, ncol=periods)
hospitalizations_unvax_sim_mat <- matrix(0, nrow=n_simulations, ncol=periods)
hospitalizations_vaxed_sim_mat <- matrix(0, nrow=n_simulations, ncol=periods)
#Added h_vec_sim matrix of 1000 rows, 4 columns
#h_vec_sim_mat <- matrix(0, nrow=n_simulations, ncol=periods)
#Also added mu_gnr_sim matrix of 1000 rows, 4 columns
#mu_gnr_sim_mat <- matrix(0, nrow=n_simulations, ncol=periods)
unvax_deaths_due_to_gnr_sim_mat <- matrix(0, nrow=n_simulations, ncol=periods)
vaxed_deaths_due_to_gnr_sim_mat <- matrix(0, nrow=n_simulations, ncol=periods)
total_deaths_due_to_gnr_sim_mat <- matrix(0, nrow=n_simulations, ncol=periods)
unvax_other_cause_deaths_sim_mat <- matrix(0, nrow=n_simulations, ncol=periods)
vaxed_other_cause_deaths_sim_mat <- matrix(0, nrow=n_simulations, ncol=periods)
total_other_cause_deaths_sim_mat <- matrix(0, nrow=n_simulations, ncol=periods)
total_deaths_sim_mat <- matrix(0, nrow=n_simulations, ncol=periods)
survived_sim_mat <- matrix(0, nrow=n_simulations, ncol=periods)
survived_vaxed_sim_mat <- matrix(0, nrow=n_simulations, ncol=periods)
survived_unvax_sim_mat <- matrix(0, nrow=n_simulations, ncol=periods)
#vm_sample_mat <- matrix(0, nrow=n_simulations, ncol=1)
#em_sample_mat <- matrix(0, nrow=n_simulations, ncol=1)

# Monte Carlo simulation
for (sim in 1:n_simulations) {

  # Sample v_m and e_m
  #vm_sample_mat[sim] <- rnorm(1, mean=vm, sd=vm_se)
  #em_sample_mat[sim] <- rnorm(1, mean=em, sd=em_se)
  
  # Ensure probabilities remain within [0, 1]
  #vm_sample_mat[sim] <- max(min(vm_sample_mat[sim], 1), 0)
  #em_sample_mat[sim] <- max(min(em_sample_mat[sim], 1), 0)
  
  remaining_cohort <- birth_cohort
  vaxed_population <- round(remaining_cohort*vm)
  unvax_population <- round(remaining_cohort*(1 - vm))
  
  for (age in 1:periods) {
    if (age == 1) {
      # Period 1: maternal vaccine confers protection to those who received maternal vaccine months'
      #h_vec_sim_mat[sim, age] <- rbeta(1, h_vec[age], h_n[age]-h_vec[age]) #For each age group, estimates a hospitalization probability based on the mean and standard deviation and 1 trial
      hospitalizations_unvax_sim_mat[sim, age] <- rbinom(1, unvax_population, h_vec_sim[sim,age]) #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
      hospitalizations_vaxed_sim_mat[sim, age] <- rbinom(1, vaxed_population, h_vec_sim[sim, age]*(1 - efficacy_mat[sim]))
      hospitalizations_total_sim_mat[sim, age] <- hospitalizations_unvax_sim_mat[sim, age] + hospitalizations_vaxed_sim_mat[sim, age]
      
      #mu_gnr_sim_mat[sim, age] <- rbeta(1, mu_gnr_vec[age], mu_gnr_n[age] - mu_gnr_vec[age])
      unvax_deaths_due_to_gnr_sim_mat[sim, age] <- rbinom(1, hospitalizations_unvax_sim_mat[sim, age], mu_gnr_sim[sim,age])
      vaxed_deaths_due_to_gnr_sim_mat[sim, age] <- rbinom(1, hospitalizations_vaxed_sim_mat[sim, age], mu_gnr_sim[sim,age])
      total_deaths_due_to_gnr_sim_mat[sim, age] <- unvax_deaths_due_to_gnr_sim_mat[sim, age] + vaxed_deaths_due_to_gnr_sim_mat[sim, age]
      unvax_other_cause_deaths_sim_mat[sim, age] <- rbinom(1, unvax_population, mu_ac[age]) - unvax_deaths_due_to_gnr_sim_mat[sim, age]
      vaxed_other_cause_deaths_sim_mat[sim, age] <- rbinom(1, vaxed_population, mu_ac[age]) - vaxed_deaths_due_to_gnr_sim_mat[sim, age]
      total_other_cause_deaths_sim_mat[sim, age] <- unvax_other_cause_deaths_sim_mat[sim, age] + vaxed_other_cause_deaths_sim_mat[sim, age]
      
      # Calculates the number of survivors
      survived_sim_mat[sim, age] <- remaining_cohort - total_deaths_due_to_gnr_sim_mat[sim, age] - total_other_cause_deaths_sim_mat[sim, age]
      survived_vaxed_sim_mat[sim, age] <- vaxed_population - vaxed_deaths_due_to_gnr_sim_mat[sim, age] - vaxed_other_cause_deaths_sim_mat[sim, age]
      survived_unvax_sim_mat[sim, age] <- unvax_population - unvax_deaths_due_to_gnr_sim_mat[sim, age] - unvax_other_cause_deaths_sim_mat[sim, age]
      total_deaths_sim_mat[sim, age] <- total_deaths_due_to_gnr_sim_mat[sim, age] + total_other_cause_deaths_sim_mat[sim, age]
      
      # Update the remaining cohort for the next period
      remaining_cohort <- survived_sim_mat[sim, age]
      
    } else {
      
      #h_vec_sim_mat[sim, age] <- rbeta(1, h_vec[age], h_n[age]-h_vec[age]) #For each age group, estimates a hospitalization probability based on the mean and standard deviation and 1 trial
      hospitalizations_total_sim_mat[sim, age] <- rbinom(1, remaining_cohort, h_vec_sim[sim,age]) #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
      #mu_gnr_sim_mat[sim, age] <- rbeta(1, mu_gnr_vec[age], mu_gnr_n[age] - mu_gnr_vec[age])
      total_deaths_due_to_gnr_sim_mat[sim, age] <- rbinom(1, hospitalizations_total_sim_mat[sim, age], mu_gnr_sim[sim,age])
      total_other_cause_deaths_sim_mat[sim, age] <- rbinom(1, remaining_cohort, mu_ac[age]) - total_deaths_due_to_gnr_sim_mat[sim, age]

      # Calculates the number of survivors
      survived_sim_mat[sim, age] <- remaining_cohort - total_deaths_due_to_gnr_sim_mat[sim, age] - total_other_cause_deaths_sim_mat[sim, age]
      total_deaths_sim_mat[sim, age] <- total_deaths_due_to_gnr_sim_mat[sim, age] + total_other_cause_deaths_sim_mat[sim, age]
      
      # Update the remaining cohort for the next period
      remaining_cohort <- survived_sim_mat[sim, age]
    }
  }
}


# Summarizes the results with confidence intervals
# First adds the total hospitalizations and deaths for each simulation
hospitalizations_sim_mat <- rowSums(hospitalizations_total_sim_mat)
deaths_due_to_gnr_sim_mat <- rowSums(total_deaths_due_to_gnr_sim_mat)
other_cause_deaths_sim_mat <- rowSums(total_other_cause_deaths_sim_mat)
all_cause_deaths_sim_mat <- rowSums(total_deaths_sim_mat)
expected_hospitalizations_mat <- mean(hospitalizations_sim_mat)
std_dev_hospitalizations_mat <- sd(hospitalizations_sim_mat)
expected_gnr_deaths_mat <- mean(deaths_due_to_gnr_sim_mat)
std_dev_gnr_deaths_mat <- sd(deaths_due_to_gnr_sim_mat)

return(data.frame(expected_hosp = expected_hospitalizations_mat, std_dev_hosp = std_dev_hospitalizations_mat, expected_gnr_deaths = expected_gnr_deaths_mat, std_dev_gnr_deaths = std_dev_gnr_deaths_mat))

#return (c(expected_hospitalizations_mat, std_dev_hospitalizations_mat))
#return (expected_hospitalizations_mat)
}

maternal_data <- simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, efficacy_mat)
