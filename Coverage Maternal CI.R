#Maternal coverage 0 -> 1.0 confidence intervals
#Creates a function to simulate the mean and standard deviation for number of hospitalizations and deaths under maternal vaccine scenario
sim_mat_coverage <- function (birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, efficacy_mat, remaining_years) {
  
  sq_hospitalizations <- c(status_quo_data_hosps_deaths[,1])
  sq_deaths <- c(status_quo_data_hosps_deaths[,2])
  
  #vm <- vm
  
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
  
  hosps_averted <- sq_hospitalizations - hospitalizations_sim_mat
  hosps_averted_p <- (sq_hospitalizations - hospitalizations_sim_mat)/(sq_hospitalizations)
  h_avert <- mean(hosps_averted) 
  h_avert_ci <- ci(hosps_averted)
  h_avert_p <- mean(hosps_averted_p)
  h_avert_p_ci <- ci(hosps_averted_p)
  
  deaths_averted <- sq_deaths - deaths_due_to_gnr_sim_mat
  deaths_averted_p <- (sq_deaths - deaths_due_to_gnr_sim_mat)/sq_deaths
  d_avert <- mean(deaths_averted)
  d_avert_ci <- ci(deaths_averted)
  d_avert_p <- mean(deaths_averted_p)
  d_avert_p_ci <- ci(deaths_averted_p)
  
  #hosps_deaths <- cbind(hospitalizations_sim_mat, deaths_due_to_gnr_sim_mat)
  
  #return(hosps_deaths)
  return(data.frame(coverage = vm, h_avert = h_avert, h_avert_low = h_avert_ci[1], h_avert_high = h_avert_ci[2],
                    h_avert_p = h_avert_p, h_avert_p_low = h_avert_p_ci[1], h_avert_p_high = h_avert_p_ci[2],
                    d_avert = d_avert, d_avert_low = d_avert_ci[1], d_avert_high = d_avert_ci[2],
                    d_avert_p = d_avert_p, d_avert_p_low = d_avert_p_ci[1], d_avert_p_high = d_avert_p_ci[2]))
  #return(total_daly_discounted_sim)
  #return(total_deaths_ci)
  #return(data.frame(expected_hosp = expected_hospitalizations_mat, std_dev_hosp = std_dev_hospitalizations_mat, expected_gnr_deaths = expected_gnr_deaths_mat, std_dev_gnr_deaths = std_dev_gnr_deaths_mat))
  #return(data.frame(expected_hosp = expected_hospitalizations_mat, std_dev_hosp = std_dev_hospitalizations_mat, expected_gnr_deaths = expected_gnr_deaths_mat, std_dev_gnr_deaths = std_dev_gnr_deaths_mat, age1_hosp = expected_hospitalizations_age1, age2_hosp = expected_hospitalizations_age2, age3_hosp = expected_hospitalizations_age3, age4_hosp = expected_hospitalizations_age4, gnr_deaths_age1 = expected_gnr_deaths_age1, gnr_deaths_age2 = expected_gnr_deaths_age2, gnr_deaths_age3 = expected_gnr_deaths_age3, gnr_deaths_age4 = expected_gnr_deaths_age4))
  
  #return (c(expected_hospitalizations_mat, std_dev_hospitalizations_mat))
  #return (expected_hospitalizations_mat)
}

mat_cov_data <- sim_mat_coverage(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, efficacy_mat, remaining_years)
mat_cov_data_cov <- sim_mat_coverage(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm=.801, efficacy_mat, remaining_years)
