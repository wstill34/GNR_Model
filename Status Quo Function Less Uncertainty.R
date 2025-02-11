
simulate_status_quo_hosps_deaths <- function (birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, remaining_years, discount_rate) {
  
  
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
      deaths_due_to_gnr_sim[sim, age] <- hospitalizations_sim[sim, age]*mu_gnr_sim[sim,age]
      
      #other_cause_deaths_sim[sim, age] <- round((remaining_cohort*mu_ac[sim, age]) - deaths_due_to_gnr_sim[sim, age])
      other_cause_deaths_sim[sim, age] <- (remaining_cohort*mu_ac[age]) - deaths_due_to_gnr_sim[sim, age]
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
  
  total_hospitalizations_ci <- ci(total_hospitalizations_sim)
  total_gnr_deaths_ci <- ci(total_deaths_due_to_gnr_sim)
  total_daly_ci <- ci(total_daly_sim)
  
  total_daly_discounted_ci <- ci(total_daly_discounted_sim)
  hosp_age1_ci <- ci(hospitalizations_sim[ ,1])
  hosp_age2_ci <- ci(hospitalizations_sim[ ,2])
  hosp_age3_ci <- ci(hospitalizations_sim[ ,3])
  hosp_age4_ci <- ci(hospitalizations_sim[ ,4])
  gnrdeaths_age1_ci <- ci(deaths_due_to_gnr_sim[ ,1])
  gnrdeaths_age2_ci <- ci(deaths_due_to_gnr_sim[ ,2])
  gnrdeaths_age3_ci <- ci(deaths_due_to_gnr_sim[ ,3])
  gnrdeaths_age4_ci <- ci(deaths_due_to_gnr_sim[ ,4])
  daly_disc_age1_ci <- ci(daly_discounted_sim[ ,1])
  daly_disc_age2_ci <- ci(daly_discounted_sim[ ,2])
  daly_disc_age3_ci <- ci(daly_discounted_sim[ ,3])
  daly_disc_age4_ci <- ci(daly_discounted_sim[ ,4])
  daly_age1_ci <- ci(daly_sim[ ,1])
  daly_age2_ci <- ci(daly_sim[ ,2])
  daly_age3_ci <- ci(daly_sim[ ,3])
  daly_age4_ci <- ci(daly_sim[ ,4])
  
  hosps_deaths <- cbind(total_hospitalizations_sim, total_deaths_due_to_gnr_sim, total_daly_discounted_sim)
  
  
  #return (expected_hospitalizations)
  #RETURNS NUMBER OF DALYS
  #return(total_daly_sim)
  #return(hist_hosp <- hospitalizations_sim)
  #return(total_gnr_deaths_ci)
  #return(data.frame(hosp_mean = expected_hospitalizations, hosp_median = expected_hospitalizations_median))
  #return(data.frame(expected_hosp = expected_hospitalizations, std_dev_hosp = std_dev_hospitalizations, expected_gnr_deaths = expected_gnr_deaths, std_dev_gnr_deaths = std_dev_gnr_deaths, age1_hosp = expected_hospitalizations_age1, age2_hosp = expected_hospitalizations_age2, age3_hosp = expected_hospitalizations_age3, age4_hosp = expected_hospitalizations_age4, gnr_deaths_age1 = expected_gnr_deaths_age1, gnr_deaths_age2 = expected_gnr_deaths_age2, gnr_deaths_age3 = expected_gnr_deaths_age3, gnr_deaths_age4 = expected_gnr_deaths_age4))
  #return(c(reduced_hospitalizations_percentage))
  #NEED TO RETURN MATRICES OR AS LIST (potentially easier) data frame where each row is a trial, but each column is what they are
  return(hosps_deaths)
  #return(data.frame(expected_hosp = expected_hospitalizations, hosp_low_ci = total_hospitalizations_ci[1], hosp_high_ci = total_hospitalizations_ci[2],
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

#generate data
status_quo_data_hosps_deaths <- simulate_status_quo_hosps_deaths(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, remaining_years, discount_rate)


#library("writexl")
#write_xlsx(status_quo_data,"C:/Users/wstil/OneDrive/Desktop/Aim 3/GNR Model/status_quo_output.xlsx", col_name=TRUE, format_headers=FALSE)