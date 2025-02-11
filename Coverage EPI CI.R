
#EPI coverage 0 -> 1.0 confidence intervals
#Creates a function to simulate the mean and standard deviation for number of hospitalizations and deaths under EPI/childhood vaccine scenario
simulate_epi_coverage <- function (birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, v1, efficacy_epi) {
  
  sq_hospitalizations <- c(status_quo_data_hosps_deaths[,1])
  sq_deaths <- c(status_quo_data_hosps_deaths[,2])
  
  # Creates matrices to store the results from simulations
  hospitalizations_total_sim_epi <- matrix(0, nrow=n_simulations, ncol=periods)
  hospitalizations_unvax_sim_epi <- matrix(0, nrow=n_simulations, ncol=periods)
  hospitalizations_vaxed_sim_epi <- matrix(0, nrow=n_simulations, ncol=periods)
  #Added h_vec_sim matrix of 1000 rows, 4 columns
  #h_vec_sim_epi <- matrix(0, nrow=n_simulations, ncol=periods)
  #Also added mu_gnr_sim matrix of 1000 rows, 4 columns
  #mu_gnr_sim_epi <- matrix(0, nrow=n_simulations, ncol=periods)
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
  daly_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  daly_discounted_sim <- matrix(0, nrow = n_simulations, ncol=periods)
  #v1_sample_epi <- matrix(0, nrow=n_simulations, ncol=1)
  #e1_sample_epi <- matrix(0, nrow=n_simulations, ncol=1)
  #remaining_cohort <- matrix(0, nrow=n_simulations, ncol=periods)
  #remaining_cohort_end <- matrix(0, nrow=n_simulations, ncol=periods)
  #unvax_population <-matrix(0, nrow=n_simulations, ncol=periods)
  #vaxed_population <- matrix(0, nrow=n_simulations, ncol=periods)
  
  
  
  # Monte Carlo simulation
  for (sim in 1:n_simulations) {
    
    # Sample v1 and e1
    #v1_sample_epi[sim] <- rnorm(1, mean=v1, sd=v1_se)
    #e1_sample_epi[sim] <- rnorm(1, mean=e1, sd=e1_se)
    #MAKE THIS A VECTOR ON PARAMETERS PAGE, CALL WHOLE VECTOR INTO THIS FUNCTION, USE SAME VECTOR ACROSS SCENARIOS
    #REMOVE VARIABILITY AROUND E1, EFFICACY FOR BOTH EPI AND MATERNAL VACCINES
    
    # Ensure probabilities remain within [0, 1]
    #v1_sample_epi[sim] <- max(min(v1_sample_epi[sim], 1), 0)
    #e1_sample_epi[sim] <- max(min(e1_sample_epi[sim], 1), 0)
    
    remaining_cohort <- birth_cohort
    vaxed_population <- round(remaining_cohort*v1)
    unvax_population <- round(remaining_cohort*(1-v1))
    for (age in 1:periods) {
      if (age == 1) {
        # Period 1: no one is vaccinated before 4 months'
        #h_vec_sim_epi[sim, age] <- rbeta(1, h_vec[age], h_n[age]-h_vec[age]) #For each age group, estimates a hospitalization probability based on the mean and standard deviation and 1 trial
        hospitalizations_unvax_sim_epi[sim, age] <- remaining_cohort*h_vec_sim[sim,age] #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
        hospitalizations_vaxed_sim_epi[sim, age] <- 0
        hospitalizations_total_sim_epi[sim, age] <- hospitalizations_unvax_sim_epi[sim, age] + hospitalizations_vaxed_sim_epi[sim, age]
        
        #mu_gnr_sim_epi[sim, age] <- rbeta(1, mu_gnr_vec[age], mu_gnr_n[age] - mu_gnr_vec[age])
        unvax_deaths_due_to_gnr_sim_epi[sim, age] <- hospitalizations_unvax_sim_epi[sim,age]*mu_gnr_sim[sim,age]
        vaxed_deaths_due_to_gnr_sim_epi[sim, age] <- 0
        total_deaths_due_to_gnr_sim_epi[sim, age] <- unvax_deaths_due_to_gnr_sim_epi[sim, age] + vaxed_deaths_due_to_gnr_sim_epi[sim, age]
        unvax_other_cause_deaths_sim_epi[sim, age] <- remaining_cohort*mu_ac[age] - unvax_deaths_due_to_gnr_sim_epi[sim, age]
        vaxed_other_cause_deaths_sim_epi[sim, age] <- 0
        total_other_cause_deaths_sim_epi[sim, age] <- unvax_other_cause_deaths_sim_epi[sim, age] + vaxed_other_cause_deaths_sim_epi[sim, age]
        
        # Calculates the number of survivors
        survived_sim_epi[sim, age] <- remaining_cohort - total_deaths_due_to_gnr_sim_epi[sim, age] - total_other_cause_deaths_sim_epi[sim, age]
        survived_vaxed_sim_epi[sim, age] <- survived_sim_epi[sim, age] * v1
        survived_unvax_sim_epi[sim, age] <- survived_sim_epi[sim, age] * (1 - v1)
        total_deaths_sim_epi[sim, age] <- total_deaths_due_to_gnr_sim_epi[sim, age] + total_other_cause_deaths_sim_epi[sim, age]
        daly_sim[sim, age] <- total_deaths_due_to_gnr_sim_epi[sim, age] * remaining_years[age]
        daly_discounted_sim[sim, age] <- total_deaths_due_to_gnr_sim_epi[sim, age] * discounted_remaining_years[age]
        
        # Update the remaining cohort for the next period
        remaining_cohort <- survived_sim_epi[sim, age]
        
        # Period 2 and beyond: some proportion vaccinated and protection conferred for 5 years
        vaxed_population <- (survived_vaxed_sim_epi[sim, age])
        unvax_population <- (survived_unvax_sim_epi[sim, age])
      } else {
        
        #h_vec_sim_epi[sim, age] <- rbeta(1, h_vec[age], h_n[age]-h_vec[age]) #For each age group, estimates a hospitalization probability based on the mean and standard deviation and 1 trial
        hospitalizations_unvax_sim_epi[sim, age] <- unvax_population*h_vec_sim[sim,age] #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
        hospitalizations_vaxed_sim_epi[sim, age] <- vaxed_population*h_vec_sim[sim,age]*(1 - efficacy_epi[sim]) #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
        hospitalizations_total_sim_epi[sim, age] <- hospitalizations_unvax_sim_epi[sim, age] + hospitalizations_vaxed_sim_epi[sim, age]
        
        #mu_gnr_sim_epi[sim, age] <- rbeta(1, mu_gnr_vec[age], mu_gnr_n[age] - mu_gnr_vec[age])
        unvax_deaths_due_to_gnr_sim_epi[sim, age] <- hospitalizations_unvax_sim_epi[sim, age]*mu_gnr_sim[sim,age]
        vaxed_deaths_due_to_gnr_sim_epi[sim, age] <- hospitalizations_vaxed_sim_epi[sim, age]*mu_gnr_sim[sim,age]
        total_deaths_due_to_gnr_sim_epi[sim, age] <- unvax_deaths_due_to_gnr_sim_epi[sim, age] + vaxed_deaths_due_to_gnr_sim_epi[sim, age]
        unvax_other_cause_deaths_sim_epi[sim, age] <- unvax_population*mu_ac[age] - unvax_deaths_due_to_gnr_sim_epi[sim, age]
        vaxed_other_cause_deaths_sim_epi[sim, age] <- vaxed_population*mu_ac[age] - vaxed_deaths_due_to_gnr_sim_epi[sim, age]
        total_other_cause_deaths_sim_epi[sim, age] <- unvax_other_cause_deaths_sim_epi[sim, age] + vaxed_other_cause_deaths_sim_epi[sim, age]
        
        # Calculates the number of survivors
        survived_unvax_sim_epi[sim, age] <- unvax_population - unvax_deaths_due_to_gnr_sim_epi[sim, age] - unvax_other_cause_deaths_sim_epi[sim, age]
        survived_vaxed_sim_epi[sim, age] <- vaxed_population - vaxed_deaths_due_to_gnr_sim_epi[sim, age] - vaxed_other_cause_deaths_sim_epi[sim, age]
        survived_sim_epi[sim, age] <- survived_unvax_sim_epi[sim, age] + survived_vaxed_sim_epi[sim, age]
        total_deaths_sim_epi[sim, age] <- total_deaths_due_to_gnr_sim_epi[sim, age] + total_other_cause_deaths_sim_epi[sim, age]
        daly_sim[sim, age] <- total_deaths_due_to_gnr_sim_epi[sim, age] * remaining_years[age]
        daly_discounted_sim[sim, age] <- total_deaths_due_to_gnr_sim_epi[sim, age] * discounted_remaining_years[age]
        
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
  expected_hospitalizations_age1 <- mean(hospitalizations_total_sim_epi[ ,1])
  expected_hospitalizations_age2 <- mean(hospitalizations_total_sim_epi[ ,2])
  expected_hospitalizations_age3 <- mean(hospitalizations_total_sim_epi[ ,3])
  expected_hospitalizations_age4 <- mean(hospitalizations_total_sim_epi[ ,4])
  std_dev_hospitalizations_epi <- sd(hospitalizations_sim_epi)
  expected_gnr_deaths <- mean(deaths_due_to_gnr_sim_epi)
  expected_gnr_deaths_age1 <- mean(total_deaths_due_to_gnr_sim_epi[ ,1])
  expected_gnr_deaths_age2 <- mean(total_deaths_due_to_gnr_sim_epi[ ,2])
  expected_gnr_deaths_age3 <- mean(total_deaths_due_to_gnr_sim_epi[ ,3])
  expected_gnr_deaths_age4 <- mean(total_deaths_due_to_gnr_sim_epi[ ,4])
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
  
    total_hospitalizations_ci <- ci(hospitalizations_sim_epi)
  total_gnr_deaths_ci <- ci(deaths_due_to_gnr_sim_epi)
  total_daly_ci <- ci(total_daly_sim)
  total_daly_discounted_ci <- ci(total_daly_discounted_sim)
  hosp_age1_ci <- ci(hospitalizations_total_sim_epi[ ,1])
  hosp_age2_ci <- ci(hospitalizations_total_sim_epi[ ,2])
  hosp_age3_ci <- ci(hospitalizations_total_sim_epi[ ,3])
  hosp_age4_ci <- ci(hospitalizations_total_sim_epi[ ,4])
  gnrdeaths_age1_ci <- ci(total_deaths_due_to_gnr_sim_epi[ ,1])
  gnrdeaths_age2_ci <- ci(total_deaths_due_to_gnr_sim_epi[ ,2])
  gnrdeaths_age3_ci <- ci(total_deaths_due_to_gnr_sim_epi[ ,3])
  gnrdeaths_age4_ci <- ci(total_deaths_due_to_gnr_sim_epi[ ,4])
  daly_disc_age1_ci <- ci(daly_discounted_sim[ ,1])
  daly_disc_age2_ci <- ci(daly_discounted_sim[ ,2])
  daly_disc_age3_ci <- ci(daly_discounted_sim[ ,3])
  daly_disc_age4_ci <- ci(daly_discounted_sim[ ,4])
  daly_age1_ci <- ci(daly_sim[ ,1])
  daly_age2_ci <- ci(daly_sim[ ,2])
  daly_age3_ci <- ci(daly_sim[ ,3])
  daly_age4_ci <- ci(daly_sim[ ,4])
  
  hosps_averted <- sq_hospitalizations - hospitalizations_sim_epi
  hosps_averted_p <- (sq_hospitalizations - hospitalizations_sim_epi)/(sq_hospitalizations)
  h_avert <- mean(hosps_averted) 
  h_avert_ci <- ci(hosps_averted)
  h_avert_p <- mean(hosps_averted_p)
  h_avert_p_ci <- ci(hosps_averted_p)
  
  deaths_averted <- sq_deaths - deaths_due_to_gnr_sim_epi
  deaths_averted_p <- (sq_deaths - deaths_due_to_gnr_sim_epi)/sq_deaths
  d_avert <- mean(deaths_averted)
  d_avert_ci <- ci(deaths_averted)
  d_avert_p <- mean(deaths_averted_p)
  d_avert_p_ci <- ci(deaths_averted_p)
  
  return(data.frame(coverage = v1, h_avert = h_avert, h_avert_low = h_avert_ci[1], h_avert_high = h_avert_ci[2],
                    h_avert_p = h_avert_p, h_avert_p_low = h_avert_p_ci[1], h_avert_p_high = h_avert_p_ci[2],
                    d_avert = d_avert, d_avert_low = d_avert_ci[1], d_avert_high = d_avert_ci[2],
                    d_avert_p = d_avert_p, d_avert_p_low = d_avert_p_ci[1], d_avert_p_high = d_avert_p_ci[2]))
  #return(total_hospitalizations_ci)
  #return(data.frame(expected_hosp = expected_hospitalizations_epi, std_dev_hosp = std_dev_hospitalizations_epi, expected_gnr_deaths = expected_gnr_deaths_epi, std_dev_gnr_deaths = std_dev_gnr_deaths_epi))
  #return(data.frame(expected_hosp = expected_hospitalizations_epi, std_dev_hosp = std_dev_hospitalizations_epi, expected_gnr_deaths = expected_gnr_deaths_epi, std_dev_gnr_deaths = std_dev_gnr_deaths_epi, age1_hosp = expected_hospitalizations_age1, age2_hosp = expected_hospitalizations_age2, age3_hosp = expected_hospitalizations_age3, age4_hosp = expected_hospitalizations_age4, gnr_deaths_age1 = expected_gnr_deaths_age1, gnr_deaths_age2 = expected_gnr_deaths_age2, gnr_deaths_age3 = expected_gnr_deaths_age3, gnr_deaths_age4 = expected_gnr_deaths_age4))
  #return(data.frame(expected_hosp = expected_hospitalizations_epi, hosp_low_ci = total_hospitalizations_ci[1], hosp_high_ci = total_hospitalizations_ci[2],
                    #hosp_age1 = expected_hospitalizations_age1, hosp_age1_low_ci = hosp_age1_ci[1], hosp_age1_high_ci = hosp_age1_ci[2],
                    #hosp_age2 = expected_hospitalizations_age2, hosp_age2_low_ci = hosp_age2_ci[1], hosp_age2_high_ci = hosp_age2_ci[2],
                    #hosp_age3 = expected_hospitalizations_age3, hosp_age3_low_ci = hosp_age3_ci[1], hosp_age3_high_ci = hosp_age3_ci[2],
                    #hosp_age4 = expected_hospitalizations_age4, hosp_age4_low_ci = hosp_age4_ci[1], hosp_age4_high_ci = hosp_age4_ci[2],
                    #gnr_deaths = expected_gnr_deaths, gnr_deaths_low_ci = total_gnr_deaths_ci[1], gnr_deaths_high_ci = total_gnr_deaths_ci[2],
                    #gnr_deaths_age1 = expected_gnr_deaths_age1, deaths_age1_low_ci = gnrdeaths_age1_ci[1], deaths_age1_high_ci = gnrdeaths_age1_ci[2],
                    #gnr_deaths_age2 = expected_gnr_deaths_age2, deaths_age2_low_ci = gnrdeaths_age2_ci[1], deaths_age2_high_ci = gnrdeaths_age2_ci[2],
                    #gnr_deaths_age3 = expected_gnr_deaths_age3, deaths_age3_low_ci = gnrdeaths_age3_ci[1], deaths_age3_high_ci = gnrdeaths_age3_ci[2],
                    #gnr_deaths_age4 = expected_gnr_deaths_age4, deaths_age4_low_ci = gnrdeaths_age4_ci[1], deaths_age4_high_ci = gnrdeaths_age4_ci[2],
                    #daly_disc = expected_daly_disc, daly_disc_low_ci = total_daly_discounted_ci[1], daly_disc_high_ci = total_daly_discounted_ci[2],
                    #daly_disc_age1 = daly_disc_age1, daly_disc_age1_low_ci = daly_disc_age1_ci[1], daly_disc_age1_high_ci = daly_disc_age1_ci[2],
                    #daly_disc_age2 = daly_disc_age2, daly_disc_age2_low_ci = daly_disc_age2_ci[1], daly_disc_age2_high_ci = daly_disc_age2_ci[2],
                    #daly_disc_age3 = daly_disc_age3, daly_disc_age3_low_ci = daly_disc_age3_ci[1], daly_disc_age3_high_ci = daly_disc_age3_ci[2],
                    #daly_disc_age4 = daly_disc_age4, daly_disc_age4_low_ci = daly_disc_age4_ci[1], daly_disc_age4_high_ci = daly_disc_age4_ci[2],
                    #daly = expected_daly, daly_low_ci = total_daly_ci[1], daly_high_ci = total_daly_ci[2],
                    #daly_age1 = daly_age1, daly_age1_low_ci = daly_age1_ci[1], daly_age1_high_ci = daly_age1_ci[2],
                    #daly_age2 = daly_age2, daly_age2_low_ci = daly_age2_ci[1], daly_age2_high_ci = daly_age2_ci[2],
                    #daly_age3 = daly_age3, daly_age3_low_ci = daly_age3_ci[1], daly_age3_high_ci = daly_age3_ci[2],
                    #daly_age4 = daly_age4, daly_age4_low_ci = daly_age4_ci[1], daly_age4_high_ci = daly_age4_ci[2]))
}

epi_data_coverage <- simulate_epi_coverage(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, v1=1, efficacy_epi)


#library("writexl")
#write_xlsx(epi_data,"C:/Users/wstil/OneDrive/Desktop/Aim 3/GNR Model/epi_output.xlsx", col_name=TRUE, format_headers=FALSE)