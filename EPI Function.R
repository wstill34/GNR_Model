#Creates a function to simulate the mean and standard deviation for number of hospitalizations and deaths under EPI/childhood vaccine scenario
simulate_epi <- function (birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, v1, v1_se, e1, e1_se) {


# Creates matrices to store the results from simulations
hospitalizations_total_sim_epi <- matrix(0, nrow=n_simulations, ncol=periods)
hospitalizations_unvax_sim_epi <- matrix(0, nrow=n_simulations, ncol=periods)
hospitalizations_vaxed_sim_epi <- matrix(0, nrow=n_simulations, ncol=periods)
#Added h_vec_sim matrix of 1000 rows, 4 columns
h_vec_sim_epi <- matrix(0, nrow=n_simulations, ncol=periods)
#Also added mu_gnr_sim matrix of 1000 rows, 4 columns
mu_gnr_sim_epi <- matrix(0, nrow=n_simulations, ncol=periods)
unvax_deaths_due_to_gnr_sim_epi <- matrix(0, nrow=n_simulations, ncol=periods)
vaxed_deaths_due_to_gnr_sim_epi <- matrix(0, nrow=n_simulations, ncol=periods)
total_deaths_due_to_gnr_sim_epi <- matrix(0, nrow=n_simulations, ncol=periods)
unvax_other_cause_deaths_sim_epi <- matrix(0, nrow=n_simulations, ncol=periods)
vaxed_other_cause_deaths_sim_epi <- matrix(0, nrow=n_simulations, ncol=periods)
total_other_cause_deaths_sim_epi <- matrix(0, nrow=n_simulations, ncol=periods)
total_deaths_sim_epi <- matrix(0, nrow=n_simulations, ncol=periods)
survived_sim_epi <- matrix(0, nrow=n_simulations, ncol=periods)
survived_vaxed_sim_epi <- matrix(0, nrow=n_simulations, ncol=periods)
survived_unvax_sim_epi <- matrix(0, nrow=n_simulations, ncol=periods)
v1_sample_epi <- matrix(0, nrow=n_simulations, ncol=1)
e1_sample_epi <- matrix(0, nrow=n_simulations, ncol=1)
#remaining_cohort <- matrix(0, nrow=n_simulations, ncol=periods)
#remaining_cohort_end <- matrix(0, nrow=n_simulations, ncol=periods)
#unvax_population <-matrix(0, nrow=n_simulations, ncol=periods)
#vaxed_population <- matrix(0, nrow=n_simulations, ncol=periods)



# Monte Carlo simulation
for (sim in 1:n_simulations) {
  
  # Sample v1 and e1
  v1_sample_epi[sim] <- rnorm(1, mean=v1, sd=v1_se)
  e1_sample_epi[sim] <- rnorm(1, mean=e1, sd=e1_se)
  #REMOVE VARIABILITY AROUND E1, EFFICACY FOR BOTH EPI AND MATERNAL VACCINES
  
  # Ensure probabilities remain within [0, 1]
  v1_sample_epi[sim] <- max(min(v1_sample_epi[sim], 1), 0)
  e1_sample_epi[sim] <- max(min(e1_sample_epi[sim], 1), 0)
  
  remaining_cohort <- birth_cohort
  vaxed_population <- round(remaining_cohort*v1_sample_epi[sim])
  unvax_population <- round(remaining_cohort*(1-v1_sample_epi[sim]))
  for (age in 1:periods) {
    if (age == 1) {
      # Period 1: no one is vaccinated before 4 months'
      h_vec_sim_epi[sim, age] <- rbeta(1, h_vec[age], h_n[age]-h_vec[age]) #For each age group, estimates a hospitalization probability based on the mean and standard deviation and 1 trial
      hospitalizations_unvax_sim_epi[sim, age] <- rbinom(1, remaining_cohort, h_vec_sim_epi[sim,age]) #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
      hospitalizations_vaxed_sim_epi[sim, age] <- 0
      hospitalizations_total_sim_epi[sim, age] <- hospitalizations_unvax_sim_epi[sim, age] + hospitalizations_vaxed_sim_epi[sim, age]
      
      mu_gnr_sim_epi[sim, age] <- rbeta(1, mu_gnr_vec[age], mu_gnr_n[age] - mu_gnr_vec[age])
      unvax_deaths_due_to_gnr_sim_epi[sim, age] <- rbinom(1, hospitalizations_unvax_sim_epi[sim, age], mu_gnr_sim_epi[sim,age])
      vaxed_deaths_due_to_gnr_sim_epi[sim, age] <- 0
      total_deaths_due_to_gnr_sim_epi[sim, age] <- unvax_deaths_due_to_gnr_sim_epi[sim, age] + vaxed_deaths_due_to_gnr_sim_epi[sim, age]
      unvax_other_cause_deaths_sim_epi[sim, age] <- rbinom(1, remaining_cohort, mu_ac[age]) - unvax_deaths_due_to_gnr_sim_epi[sim, age]
      vaxed_other_cause_deaths_sim_epi[sim, age] <- 0
      total_other_cause_deaths_sim_epi[sim, age] <- unvax_other_cause_deaths_sim_epi[sim, age] + vaxed_other_cause_deaths_sim_epi[sim, age]
      
      # Calculates the number of survivors
      survived_sim_epi[sim, age] <- remaining_cohort - total_deaths_due_to_gnr_sim_epi[sim, age] - total_other_cause_deaths_sim_epi[sim, age]
      survived_vaxed_sim_epi[sim, age] <- survived_sim_epi[sim, age] * v1_sample_epi[sim]
      survived_unvax_sim_epi[sim, age] <- survived_sim_epi[sim, age] * (1 - v1_sample_epi[sim])
      total_deaths_sim_epi[sim, age] <- total_deaths_due_to_gnr_sim_epi[sim, age] + total_other_cause_deaths_sim_epi[sim, age]
      
      # Update the remaining cohort for the next period
      remaining_cohort <- survived_sim_epi[sim, age]
      
      # Period 2 and beyond: some proportion vaccinated and protection conferred for 5 years
      vaxed_population <- round(survived_vaxed_sim_epi[sim, age])
      unvax_population <- round(survived_unvax_sim_epi[sim, age])
    } else {
      
      h_vec_sim_epi[sim, age] <- rbeta(1, h_vec[age], h_n[age]-h_vec[age]) #For each age group, estimates a hospitalization probability based on the mean and standard deviation and 1 trial
      hospitalizations_unvax_sim_epi[sim, age] <- rbinom(1, unvax_population, h_vec_sim_epi[sim,age]) #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
      hospitalizations_vaxed_sim_epi[sim, age] <- rbinom(1, vaxed_population, h_vec_sim_epi[sim,age]*(1 - e1_sample_epi[sim])) #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
      hospitalizations_total_sim_epi[sim, age] <- hospitalizations_unvax_sim_epi[sim, age] + hospitalizations_vaxed_sim_epi[sim, age]
      
      mu_gnr_sim_epi[sim, age] <- rbeta(1, mu_gnr_vec[age], mu_gnr_n[age] - mu_gnr_vec[age])
      unvax_deaths_due_to_gnr_sim_epi[sim, age] <- rbinom(1, hospitalizations_unvax_sim_epi[sim, age], mu_gnr_sim_epi[sim,age])
      vaxed_deaths_due_to_gnr_sim_epi[sim, age] <- rbinom(1, hospitalizations_vaxed_sim_epi[sim, age], mu_gnr_sim_epi[sim,age])
      total_deaths_due_to_gnr_sim_epi[sim, age] <- unvax_deaths_due_to_gnr_sim_epi[sim, age] + vaxed_deaths_due_to_gnr_sim_epi[sim, age]
      unvax_other_cause_deaths_sim_epi[sim, age] <- rbinom(1, unvax_population, mu_ac[age]) - unvax_deaths_due_to_gnr_sim_epi[sim, age]
      vaxed_other_cause_deaths_sim_epi[sim, age] <- rbinom(1, vaxed_population, mu_ac[age]) - vaxed_deaths_due_to_gnr_sim_epi[sim, age]
      total_other_cause_deaths_sim_epi[sim, age] <- unvax_other_cause_deaths_sim_epi[sim, age] + vaxed_other_cause_deaths_sim_epi[sim, age]
      
      # Calculates the number of survivors
      survived_unvax_sim_epi[sim, age] <- unvax_population - unvax_deaths_due_to_gnr_sim_epi[sim, age] - unvax_other_cause_deaths_sim_epi[sim, age]
      survived_vaxed_sim_epi[sim, age] <- vaxed_population - vaxed_deaths_due_to_gnr_sim_epi[sim, age] - vaxed_other_cause_deaths_sim_epi[sim, age]
      survived_sim_epi[sim, age] <- survived_unvax_sim_epi[sim, age] + survived_vaxed_sim_epi[sim, age]
      total_deaths_sim_epi[sim, age] <- total_deaths_due_to_gnr_sim_epi[sim, age] + total_other_cause_deaths_sim_epi[sim, age]
      
      # Update the remaining cohort for the next period
      unvax_population <- round(survived_unvax_sim_epi[sim, age])
      vaxed_population <- round(survived_vaxed_sim_epi[sim, age])
    }
  }
}


# Summarizes the results with confidence intervals
# First adds the total hospitalizations and deaths for each simulation
hospitalizations_sim_epi <- rowSums(hospitalizations_total_sim_epi)
deaths_due_to_gnr_sim_epi <- rowSums(total_deaths_due_to_gnr_sim_epi)
other_cause_deaths_sim_epi <- rowSums(total_other_cause_deaths_sim_epi)
all_cause_deaths_sim_epi <- rowSums(total_deaths_sim_epi)
expected_hospitalizations_epi <- mean(hospitalizations_sim_epi)
std_dev_hospitalizations_epi <- sd(hospitalizations_sim_epi)
expected_gnr_deaths_epi <- mean(deaths_due_to_gnr_sim_epi)
std_dev_gnr_deaths_epi <- sd(deaths_due_to_gnr_sim_epi)

return(data.frame(expected_hosp = expected_hospitalizations_epi, std_dev_hosp = std_dev_hospitalizations_epi, expected_gnr_deaths = expected_gnr_deaths_epi, std_dev_gnr_deaths = std_dev_gnr_deaths_epi))
#return (expected_hospitalizations_epi)
}

epi_data <- simulate_epi(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, v1, v1_se, e1, e1_se)
