
simulate_status_quo_daly_less_variability <- function (birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, remaining_years, discount_rate) {
  
  
  # Creates matrices to store the results from simulations
  hospitalizations_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  #h_vec_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  #use same for each trial across scenarios
  #mu_gnr_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  #use same for each trial across 
  deaths_due_to_gnr_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  other_cause_deaths_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  total_deaths_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  survived_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  daly_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  daly_discounted_sim <- matrix(0, nrow = n_simulations, ncol=periods)
  
  
  # Monte Carlo simulation
  for (sim in 1:n_simulations) {
    remaining_cohort <- birth_cohort
    
    for (age in 1:periods) {
      # Calculates number of hospitalizations, deaths due to GNR, and other cause deaths
      #h_vec_sim[sim, age] <- rbeta(1, h_vec[age], h_n[age]-h_vec[age]) #For each age group, estimates a hospitalization probability based on the mean and standard deviation and 1 trial
      hospitalizations_sim[sim, age] <- remaining_cohort*h_vec_sim[sim,age]
      #hospitalizations_sim[sim, age] <- rbinom(1, remaining_cohort, h_vec_sim[sim,age]) #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
      #mu_gnr_sim[sim, age] <- rbeta(1, mu_gnr_vec[age], mu_gnr_n[age] - mu_gnr_vec[age])
      #deaths_due_to_gnr_sim[sim, age] <- rbinom(1, hospitalizations_sim[sim, age], mu_gnr_sim[sim,age])
      deaths_due_to_gnr_sim[sim, age] <- round(hospitalizations_sim[sim, age]*mu_gnr_sim[sim,age])
      
      #other_cause_deaths_sim[sim, age] <- round((remaining_cohort*mu_ac[sim, age]) - deaths_due_to_gnr_sim[sim, age])
      other_cause_deaths_sim[sim, age] <- round((remaining_cohort*mu_ac[age]) - deaths_due_to_gnr_sim[sim, age])
      #other_cause_deaths_sim[sim, age] <- rbinom(1, remaining_cohort, mu_ac[age]) - deaths_due_to_gnr_sim[sim, age]
      daly_sim[sim, age] <- deaths_due_to_gnr_sim[sim, age] * remaining_years[age]
      daly_discounted_sim[sim, age] <- deaths_due_to_gnr_sim[sim, age] * discounted_remaining_years[age]
      
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
  expected_hospitalizations_median <- median(total_hospitalizations_sim)
  expected_hospitalizations_age1 <- mean(hospitalizations_sim[ ,1])
  expected_hospitalizations_age2 <- mean(hospitalizations_sim[ ,2])
  expected_hospitalizations_age3 <- mean(hospitalizations_sim[ ,3])
  expected_hospitalizations_age4 <- mean(hospitalizations_sim[ ,4])
  std_dev_hospitalizations <- sd(total_hospitalizations_sim)
  expected_gnr_deaths <- mean(total_deaths_due_to_gnr_sim)
  std_dev_gnr_deaths <- sd(total_deaths_due_to_gnr_sim)
  expected_gnr_deaths_age1 <- mean(deaths_due_to_gnr_sim[ ,1])
  expected_gnr_deaths_age2 <- mean(deaths_due_to_gnr_sim[ ,2])
  expected_gnr_deaths_age3 <- mean(deaths_due_to_gnr_sim[ ,3])
  expected_gnr_deaths_age4 <- mean(deaths_due_to_gnr_sim[ ,4])
  total_daly_sim <- rowSums(daly_sim)
  total_daly_discounted_sim <- rowSums(daly_discounted_sim)
  
  total_hospitalizations_ci <- ci(total_hospitalizations_sim)
  total_gnr_deaths_ci <- ci(total_deaths_due_to_gnr_sim)
  #don't calculate this std dev
  
  #return (c(expected_hospitalizations, std_dev_hospitalizations))
  #return (expected_hospitalizations)
  return(total_daly_discounted_sim)
  #return(hist_hosp <- hospitalizations_sim)
  #return(total_gnr_deaths_ci)
  #return(data.frame(hosp_mean = expected_hospitalizations, hosp_median = expected_hospitalizations_median))
  #return(data.frame(expected_hosp = expected_hospitalizations, std_dev_hosp = std_dev_hospitalizations, expected_gnr_deaths = expected_gnr_deaths, std_dev_gnr_deaths = std_dev_gnr_deaths, age1_hosp = expected_hospitalizations_age1, age2_hosp = expected_hospitalizations_age2, age3_hosp = expected_hospitalizations_age3, age4_hosp = expected_hospitalizations_age4, gnr_deaths_age1 = expected_gnr_deaths_age1, gnr_deaths_age2 = expected_gnr_deaths_age2, gnr_deaths_age3 = expected_gnr_deaths_age3, gnr_deaths_age4 = expected_gnr_deaths_age4))
  #return(c(reduced_hospitalizations_percentage))
  #NEED TO RETURN MATRICES OR AS LIST (potentially easier) data frame where each row is a trial, but each column is what they are
  #return(data.frame(expected_hosp = expected_hospitalizations, hosp_low_ci = total_hospitalizations_ci[1], hosp_high_ci = total_hospitalizations_ci = total_hospitalizations_ci[2], 
  
  #expected_gnr_deaths = expected_gnr_deaths, gnr_deaths_low_ci = totalstd_dev_gnr_deaths = std_dev_gnr_deaths, age1_hosp = expected_hospitalizations_age1, age2_hosp = expected_hospitalizations_age2, age3_hosp = expected_hospitalizations_age3, age4_hosp = expected_hospitalizations_age4, gnr_deaths_age1 = expected_gnr_deaths_age1, gnr_deaths_age2 = expected_gnr_deaths_age2, gnr_deaths_age3 = expected_gnr_deaths_age3, gnr_deaths_age4 = expected_gnr_deaths_age4))
  
}

#generate data
status_quo_data_daly_less_variability <- simulate_status_quo_daly_less_variability(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, remaining_years, discount_rate)






#Creates a function to simulate the mean and standard deviation for number of hospitalizations and deaths under maternal vaccine scenario
simulate_maternal_daly_less_variability <- function (birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, efficacy_mat, remaining_years) {
  
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
  daly_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  daly_discounted_sim <- matrix(0, nrow = n_simulations, ncol=periods)
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
        hospitalizations_unvax_sim_mat[sim, age] <- (unvax_population*h_vec_sim[sim,age])
        #hospitalizations_unvax_sim_mat[sim, age] <- rbinom(1, unvax_population, h_vec_sim[sim,age]) #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
        hospitalizations_vaxed_sim_mat[sim, age] <- (vaxed_population*h_vec_sim[sim,age]*(1 - efficacy_mat[sim]))
        #hospitalizations_vaxed_sim_mat[sim, age] <- rbinom(1, vaxed_population, h_vec_sim[sim, age]*(1 - efficacy_mat[sim]))
        hospitalizations_total_sim_mat[sim, age] <- hospitalizations_unvax_sim_mat[sim, age] + hospitalizations_vaxed_sim_mat[sim, age]
        
        #mu_gnr_sim_mat[sim, age] <- rbeta(1, mu_gnr_vec[age], mu_gnr_n[age] - mu_gnr_vec[age])
        unvax_deaths_due_to_gnr_sim_mat[sim, age] <- (hospitalizations_unvax_sim_mat[sim, age]*mu_gnr_sim[sim,age])
        #unvax_deaths_due_to_gnr_sim_mat[sim, age] <- rbinom(1, hospitalizations_unvax_sim_mat[sim, age], mu_gnr_sim[sim,age])
        vaxed_deaths_due_to_gnr_sim_mat[sim, age] <- (hospitalizations_vaxed_sim_mat[sim, age]*mu_gnr_sim[sim,age])
        #vaxed_deaths_due_to_gnr_sim_mat[sim, age] <- rbinom(1, hospitalizations_vaxed_sim_mat[sim, age], mu_gnr_sim[sim,age])
        total_deaths_due_to_gnr_sim_mat[sim, age] <- unvax_deaths_due_to_gnr_sim_mat[sim, age] + vaxed_deaths_due_to_gnr_sim_mat[sim, age]
        unvax_other_cause_deaths_sim_mat[sim, age] <- ((unvax_population*mu_ac[age]) - unvax_deaths_due_to_gnr_sim_mat[sim, age])
        # removing for now unvax_other_cause_deaths_sim_mat[sim, age] <- rbinom(1, unvax_population, mu_ac[age]) - unvax_deaths_due_to_gnr_sim_mat[sim, age]
        vaxed_other_cause_deaths_sim_mat[sim, age] <- ((vaxed_population*mu_ac[age]) - vaxed_deaths_due_to_gnr_sim_mat[sim, age])
        # this code is sampling all-cause death: vaxed_other_cause_deaths_sim_mat[sim, age] <- rbinom(1, vaxed_population, mu_ac[age]) - vaxed_deaths_due_to_gnr_sim_mat[sim, age]
        total_other_cause_deaths_sim_mat[sim, age] <- unvax_other_cause_deaths_sim_mat[sim, age] + vaxed_other_cause_deaths_sim_mat[sim, age]
        
        # Calculates the number of survivors
        survived_sim_mat[sim, age] <- remaining_cohort - total_deaths_due_to_gnr_sim_mat[sim, age] - total_other_cause_deaths_sim_mat[sim, age]
        survived_vaxed_sim_mat[sim, age] <- vaxed_population - vaxed_deaths_due_to_gnr_sim_mat[sim, age] - vaxed_other_cause_deaths_sim_mat[sim, age]
        survived_unvax_sim_mat[sim, age] <- unvax_population - unvax_deaths_due_to_gnr_sim_mat[sim, age] - unvax_other_cause_deaths_sim_mat[sim, age]
        total_deaths_sim_mat[sim, age] <- total_deaths_due_to_gnr_sim_mat[sim, age] + total_other_cause_deaths_sim_mat[sim, age]
        daly_sim[sim, age] <- total_deaths_due_to_gnr_sim_mat[sim, age] * remaining_years[age]
        daly_discounted_sim[sim, age] <- total_deaths_due_to_gnr_sim_mat[sim, age] * discounted_remaining_years[age]
        
        # Update the remaining cohort for the next period
        remaining_cohort <- survived_sim_mat[sim, age]
        
      } else {
        
        #h_vec_sim_mat[sim, age] <- rbeta(1, h_vec[age], h_n[age]-h_vec[age]) #For each age group, estimates a hospitalization probability based on the mean and standard deviation and 1 trial
        hospitalizations_total_sim_mat[sim, age] <- (remaining_cohort*h_vec_sim[sim,age])
        #hospitalizations_total_sim_mat[sim, age] <- rbinom(1, remaining_cohort, h_vec_sim[sim,age]) #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
        #mu_gnr_sim_mat[sim, age] <- rbeta(1, mu_gnr_vec[age], mu_gnr_n[age] - mu_gnr_vec[age])
        total_deaths_due_to_gnr_sim_mat[sim, age] <- (hospitalizations_total_sim_mat[sim, age]*mu_gnr_sim[sim,age])
        #total_deaths_due_to_gnr_sim_mat[sim, age] <- rbinom(1, hospitalizations_total_sim_mat[sim, age], mu_gnr_sim[sim,age])
        total_other_cause_deaths_sim_mat[sim, age] <- ((remaining_cohort*mu_ac[age]) - total_deaths_due_to_gnr_sim_mat[sim, age])
        # samples all-cause death total_other_cause_deaths_sim_mat[sim, age] <- rbinom(1, remaining_cohort, mu_ac[age]) - total_deaths_due_to_gnr_sim_mat[sim, age]
        daly_sim[sim, age] <- total_deaths_due_to_gnr_sim_mat[sim, age] * remaining_years[age]
        daly_discounted_sim[sim, age] <- total_deaths_due_to_gnr_sim_mat[sim, age] * discounted_remaining_years[age]
        
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
  expected_hospitalizations_age1 <- mean(hospitalizations_total_sim_mat[ ,1])
  expected_hospitalizations_age2 <- mean(hospitalizations_total_sim_mat[ ,2])
  expected_hospitalizations_age3 <- mean(hospitalizations_total_sim_mat[ ,3])
  expected_hospitalizations_age4 <- mean(hospitalizations_total_sim_mat[ ,4])
  expected_gnr_deaths <- mean(deaths_due_to_gnr_sim_mat)
  expected_gnr_deaths_age1 <- mean(total_deaths_due_to_gnr_sim_mat[ ,1])
  expected_gnr_deaths_age2 <- mean(total_deaths_due_to_gnr_sim_mat[ ,2])
  expected_gnr_deaths_age3 <- mean(total_deaths_due_to_gnr_sim_mat[ ,3])
  expected_gnr_deaths_age4 <- mean(total_deaths_due_to_gnr_sim_mat[ ,4])
  daly_disc_age1 <- mean(daly_discounted_sim[ ,1])
  daly_disc_age2 <- mean(daly_discounted_sim[ ,2])
  daly_disc_age3 <- mean(daly_discounted_sim[ ,3])
  daly_disc_age4 <- mean(daly_discounted_sim[ ,4])
  daly_age1 <- mean(daly_sim[ ,1])
  daly_age2 <- mean(daly_sim[ ,2])
  daly_age3 <- mean(daly_sim[ ,3])
  daly_age4 <- mean(daly_sim[ ,4])
  total_daly_sim <- rowSums(daly_sim)
  total_daly_discounted_sim <- rowSums(daly_discounted_sim)
  expected_daly_disc <- mean(total_daly_discounted_sim)
  expected_daly <- mean(total_daly_sim)
  
  total_hospitalizations_ci <- ci(hospitalizations_sim_mat)
  total_gnr_deaths_ci <- ci(deaths_due_to_gnr_sim_mat)
  total_daly_ci <- ci(total_daly_sim)
  total_daly_discounted_ci <- ci(total_daly_discounted_sim)
  hosp_age1_ci <- ci(hospitalizations_total_sim_mat[ ,1])
  hosp_age2_ci <- ci(hospitalizations_total_sim_mat[ ,2])
  hosp_age3_ci <- ci(hospitalizations_total_sim_mat[ ,3])
  hosp_age4_ci <- ci(hospitalizations_total_sim_mat[ ,4])
  gnrdeaths_age1_ci <- ci(total_deaths_due_to_gnr_sim_mat[ ,1])
  gnrdeaths_age2_ci <- ci(total_deaths_due_to_gnr_sim_mat[ ,2])
  gnrdeaths_age3_ci <- ci(total_deaths_due_to_gnr_sim_mat[ ,3])
  gnrdeaths_age4_ci <- ci(total_deaths_due_to_gnr_sim_mat[ ,4])
  daly_disc_age1_ci <- ci(daly_discounted_sim[ ,1])
  daly_disc_age2_ci <- ci(daly_discounted_sim[ ,2])
  daly_disc_age3_ci <- ci(daly_discounted_sim[ ,3])
  daly_disc_age4_ci <- ci(daly_discounted_sim[ ,4])
  daly_age1_ci <- ci(daly_sim[ ,1])
  daly_age2_ci <- ci(daly_sim[ ,2])
  daly_age3_ci <- ci(daly_sim[ ,3])
  daly_age4_ci <- ci(daly_sim[ ,4])
  
  hosps_deaths <- cbind(hospitalizations_sim_mat, deaths_due_to_gnr_sim_mat)
  
  return(hosps_deaths)
  
  #return(total_daly_discounted_sim)
  #return(total_deaths_ci)
  #return(data.frame(expected_hosp = expected_hospitalizations_mat, std_dev_hosp = std_dev_hospitalizations_mat, expected_gnr_deaths = expected_gnr_deaths_mat, std_dev_gnr_deaths = std_dev_gnr_deaths_mat))
  #return(data.frame(expected_hosp = expected_hospitalizations_mat, std_dev_hosp = std_dev_hospitalizations_mat, expected_gnr_deaths = expected_gnr_deaths_mat, std_dev_gnr_deaths = std_dev_gnr_deaths_mat, age1_hosp = expected_hospitalizations_age1, age2_hosp = expected_hospitalizations_age2, age3_hosp = expected_hospitalizations_age3, age4_hosp = expected_hospitalizations_age4, gnr_deaths_age1 = expected_gnr_deaths_age1, gnr_deaths_age2 = expected_gnr_deaths_age2, gnr_deaths_age3 = expected_gnr_deaths_age3, gnr_deaths_age4 = expected_gnr_deaths_age4))
  
  #return (c(expected_hospitalizations_mat, std_dev_hospitalizations_mat))
  #return (expected_hospitalizations_mat)
}

maternal_data_daly_less_variability_matrix <- simulate_maternal_daly_less_variability(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, efficacy_mat, remaining_years)


#daly_diff_maternal_less_variability <- status_quo_data_daly_less_variability - maternal_data_daly_less_variability
#hist(daly_diff_maternal_less_variability)
#mean(daly_diff_maternal_less_variability)
#ci(daly_diff_maternal_less_variability)

