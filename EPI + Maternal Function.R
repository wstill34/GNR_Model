#Creates a function to simulate the mean and standard deviation for number of hospitalizations and deaths under joint EPI/childhood + maternal vaccine scenario
simulate_epi_maternal <- function (birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, vm_se, em, em_se, v1, v1_se, e1, e1_se) {


    # Creates matrices to store the results from simulations
  hospitalizations_total_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  hospitalizations_unvax_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  hospitalizations_vaxed_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  #Added h_vec_sim matrix of 1000 rows, 4 columns
  h_vec_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  #Also added mu_gnr_sim matrix of 1000 rows, 4 columns
  mu_gnr_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  unvax_deaths_due_to_gnr_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  vaxed_deaths_due_to_gnr_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  total_deaths_due_to_gnr_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  unvax_other_cause_deaths_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  vaxed_other_cause_deaths_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  total_other_cause_deaths_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  total_deaths_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  survived_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  survived_vaxed_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  survived_unvax_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  v1_sample <- matrix(0, nrow=n_simulations, ncol=1)
  e1_sample <- matrix(0, nrow=n_simulations, ncol=1)
  vm_sample <- matrix(0, nrow=n_simulations, ncol=1)
  em_sample <- matrix(0, nrow=n_simulations, ncol=1)
  
  
  # Monte Carlo simulation
  for (sim in 1:n_simulations) {
    
    # Sample v_m and e_m
    vm_sample[sim] <- rnorm(1, mean=vm, sd=vm_se)
    em_sample[sim] <- rnorm(1, mean=em, sd=em_se)
    
    # Ensure probabilities remain within [0, 1]
    vm_sample[sim] <- max(min(vm_sample[sim], 1), 0)
    em_sample[sim] <- max(min(em_sample[sim], 1), 0)
    
    # Sample v1 and e1
    v1_sample[sim] <- rnorm(1, mean=v1, sd=v1_se)
    e1_sample[sim] <- rnorm(1, mean=e1, sd=e1_se)
    
    # Ensure probabilities remain within [0, 1]
    v1_sample[sim] <- max(min(v1_sample[sim], 1), 0)
    e1_sample[sim] <- max(min(e1_sample[sim], 1), 0)
    
    remaining_cohort <- birth_cohort
    vaxed_population <- round(remaining_cohort*vm_sample[sim])
    unvax_population <- round(remaining_cohort*(1 - vm_sample[sim]))
    
    for (age in 1:periods) {
      if (age == 1) {
        # Period 1: maternal vaccine confers protection to those who received maternal vaccine'
        h_vec_sim[sim, age] <- rbeta(1, h_vec[age], h_n[age]-h_vec[age]) #For each age group, estimates a hospitalization probability based on the mean and standard deviation and 1 trial
        hospitalizations_unvax_sim[sim, age] <- rbinom(1, unvax_population, h_vec_sim[sim,age]) #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
        hospitalizations_vaxed_sim[sim, age] <- rbinom(1, vaxed_population, h_vec_sim[sim, age]*(1 - em_sample[sim]))
        hospitalizations_total_sim[sim, age] <- hospitalizations_unvax_sim[sim, age] + hospitalizations_vaxed_sim[sim, age]
        
        mu_gnr_sim[sim, age] <- rbeta(1, mu_gnr_vec[age], mu_gnr_n[age] - mu_gnr_vec[age])
        unvax_deaths_due_to_gnr_sim[sim, age] <- rbinom(1, hospitalizations_unvax_sim[sim, age], mu_gnr_sim[sim,age])
        vaxed_deaths_due_to_gnr_sim[sim, age] <- rbinom(1, hospitalizations_vaxed_sim[sim, age], mu_gnr_sim[sim,age])
        total_deaths_due_to_gnr_sim[sim, age] <- unvax_deaths_due_to_gnr_sim[sim, age] + vaxed_deaths_due_to_gnr_sim[sim, age]
        unvax_other_cause_deaths_sim[sim, age] <- rbinom(1, unvax_population, mu_ac[age]) - unvax_deaths_due_to_gnr_sim[sim, age]
        vaxed_other_cause_deaths_sim[sim, age] <- rbinom(1, vaxed_population, mu_ac[age]) - vaxed_deaths_due_to_gnr_sim[sim, age]
        total_other_cause_deaths_sim[sim, age] <- unvax_other_cause_deaths_sim[sim, age] + vaxed_other_cause_deaths_sim[sim, age]
        
        # Calculates the number of survivors
        survived_sim[sim, age] <- remaining_cohort - total_deaths_due_to_gnr_sim[sim, age] - total_other_cause_deaths_sim[sim, age]
        survived_vaxed_sim[sim, age] <- vaxed_population - vaxed_deaths_due_to_gnr_sim[sim, age] - vaxed_other_cause_deaths_sim[sim, age]
        survived_unvax_sim[sim, age] <- unvax_population - unvax_deaths_due_to_gnr_sim[sim, age] - unvax_other_cause_deaths_sim[sim, age]
        total_deaths_sim[sim, age] <- total_deaths_due_to_gnr_sim[sim, age] + total_other_cause_deaths_sim[sim, age]
        
        # Update the remaining cohort for the next period
        remaining_cohort <- survived_sim[sim, age]
        vaxed_population <- round(remaining_cohort*v1_sample[sim])
        unvax_population <- round(remaining_cohort*(1-v1_sample[sim]))
        
      } else {
        h_vec_sim[sim, age] <- rbeta(1, h_vec[age], h_n[age]-h_vec[age]) #For each age group, estimates a hospitalization probability based on the mean and standard deviation and 1 trial
        hospitalizations_unvax_sim[sim, age] <- rbinom(1, unvax_population, h_vec_sim[sim,age]) #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
        hospitalizations_vaxed_sim[sim, age] <- rbinom(1, vaxed_population, h_vec_sim[sim,age]*(1 - e1_sample[sim])) #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
        hospitalizations_total_sim[sim, age] <- hospitalizations_unvax_sim[sim, age] + hospitalizations_vaxed_sim[sim, age]
        
        mu_gnr_sim[sim, age] <- rbeta(1, mu_gnr_vec[age], mu_gnr_n[age] - mu_gnr_vec[age])
        unvax_deaths_due_to_gnr_sim[sim, age] <- rbinom(1, hospitalizations_unvax_sim[sim, age], mu_gnr_sim[sim,age])
        vaxed_deaths_due_to_gnr_sim[sim, age] <- rbinom(1, hospitalizations_vaxed_sim[sim, age], mu_gnr_sim[sim,age])
        total_deaths_due_to_gnr_sim[sim, age] <- unvax_deaths_due_to_gnr_sim[sim, age] + vaxed_deaths_due_to_gnr_sim[sim, age]
        unvax_other_cause_deaths_sim[sim, age] <- rbinom(1, unvax_population, mu_ac[age]) - unvax_deaths_due_to_gnr_sim[sim, age]
        vaxed_other_cause_deaths_sim[sim, age] <- rbinom(1, vaxed_population, mu_ac[age]) - vaxed_deaths_due_to_gnr_sim[sim, age]
        total_other_cause_deaths_sim[sim, age] <- unvax_other_cause_deaths_sim[sim, age] + vaxed_other_cause_deaths_sim[sim, age]
        
        # Calculates the number of survivors
        survived_unvax_sim[sim, age] <- unvax_population - unvax_deaths_due_to_gnr_sim[sim, age] - unvax_other_cause_deaths_sim[sim, age]
        survived_vaxed_sim[sim, age] <- vaxed_population - vaxed_deaths_due_to_gnr_sim[sim, age] - vaxed_other_cause_deaths_sim[sim, age]
        survived_sim[sim, age] <- survived_unvax_sim[sim, age] + survived_vaxed_sim[sim, age]
        total_deaths_sim[sim, age] <- total_deaths_due_to_gnr_sim[sim, age] + total_other_cause_deaths_sim[sim, age]
        
        # Update the remaining cohort for the next period
        unvax_population <- round(survived_unvax_sim[sim, age])
        vaxed_population <- round(survived_vaxed_sim[sim, age])
      }
    }
  }
  
  
  # Summarizes the results with confidence intervals
  # First adds the total hospitalizations and deaths for each simulation
  hospitalizations_sim <- rowSums(hospitalizations_total_sim)
  deaths_due_to_gnr_sim <- rowSums(total_deaths_due_to_gnr_sim)
  other_cause_deaths_sim <- rowSums(total_other_cause_deaths_sim)
  all_cause_deaths_sim <- rowSums(total_deaths_sim)
  expected_hospitalizations <- mean(hospitalizations_sim)
  std_dev_hospitalizations <- sd(hospitalizations_sim)
  expected_gnr_deaths <- mean(deaths_due_to_gnr_sim)
  std_dev_gnr_deaths <- sd(deaths_due_to_gnr_sim)
  
  return(data.frame(expected_hosp = expected_hospitalizations, std_dev_hosp = std_dev_hospitalizations, expected_gnr_deaths = expected_gnr_deaths, std_dev_gnr_deaths = std_dev_gnr_deaths))
  
  #return (c(expected_hospitalizations, std_dev_hospitalizations))
}

joint_data <- simulate_epi_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, vm_se, em, em_se, v1, v1_se, e1, e1_se)
  