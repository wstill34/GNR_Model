#Creates a function to simulate the status quo outputs
simulate_sq <- function (birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations) {
  

# Creates matrices to store the results from simulations
hospitalizations_sim <- matrix(0, nrow=n_simulations, ncol=periods)
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
    deaths_due_to_gnr_sim[sim, age] <- hospitalizations_sim[sim,age]*mu_gnr_sim[sim,age]
    other_cause_deaths_sim[sim, age] <- remaining_cohort*mu_ac[age] - deaths_due_to_gnr_sim[sim, age]
     
   # Calculates the number of survivors
    survived_sim[sim, age] <- remaining_cohort - deaths_due_to_gnr_sim[sim, age] - other_cause_deaths_sim[sim, age]
    total_deaths_sim[sim, age] <- deaths_due_to_gnr_sim[sim, age] + other_cause_deaths_sim[sim, age]
    daly_sim[sim, age] <- deaths_due_to_gnr_sim[sim, age] * remaining_years[age]
    daly_discounted_sim[sim, age] <- deaths_due_to_gnr_sim[sim, age] * discounted_remaining_years[age]
   
    # Update the remaining cohort for the next period
    remaining_cohort <- survived_sim[sim, age]
  }
}


# Summarizes the results with confidence intervals
# Add the total hospitalizations, deaths, and DALYs for each simulation
total_hospitalizations_sim <- rowSums(hospitalizations_sim)
total_deaths_due_to_gnr_sim <- rowSums(deaths_due_to_gnr_sim)
total_other_cause_deaths_sim <- rowSums(other_cause_deaths_sim)
total_all_cause_deaths_sim <- rowSums(total_deaths_sim)
total_daly_sim <- rowSums(daly_sim)
total_daly_discounted_sim <- rowSums(daly_discounted_sim)

#Expected hospitalizations, overall and by age group
expected_hospitalizations <- mean(total_hospitalizations_sim)
expected_hospitalizations_age1 <- mean(hospitalizations_sim[ ,1])
expected_hospitalizations_age2 <- mean(hospitalizations_sim[ ,2])
expected_hospitalizations_age3 <- mean(hospitalizations_sim[ ,3])
expected_hospitalizations_age4 <- mean(hospitalizations_sim[ ,4])

#Expected gnr deaths, overall and by age group
expected_gnr_deaths <- mean(total_deaths_due_to_gnr_sim)
expected_gnr_deaths_age1 <- mean(deaths_due_to_gnr_sim[ ,1])
expected_gnr_deaths_age2 <- mean(deaths_due_to_gnr_sim[ ,2])
expected_gnr_deaths_age3 <- mean(deaths_due_to_gnr_sim[ ,3])
expected_gnr_deaths_age4 <- mean(deaths_due_to_gnr_sim[ ,4])

#Expected DALYs, overall and by age group
expected_daly <- mean(total_daly_sim)
daly_age1 <- mean(daly_sim[ ,1])
daly_age2 <- mean(daly_sim[ ,2])
daly_age3 <- mean(daly_sim[ ,3])
daly_age4 <- mean(daly_sim[ ,4])
expected_daly_disc <- mean(total_daly_discounted_sim)
daly_disc_age1 <- mean(daly_discounted_sim[ ,1])
daly_disc_age2 <- mean(daly_discounted_sim[ ,2])
daly_disc_age3 <- mean(daly_discounted_sim[ ,3])
daly_disc_age4 <- mean(daly_discounted_sim[ ,4])

#Confidence Intervals
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
daly_age1_ci <- ci(daly_sim[ ,1])
daly_age2_ci <- ci(daly_sim[ ,2])
daly_age3_ci <- ci(daly_sim[ ,3])
daly_age4_ci <- ci(daly_sim[ ,4])
daly_disc_age1_ci <- ci(daly_discounted_sim[ ,1])
daly_disc_age2_ci <- ci(daly_discounted_sim[ ,2])
daly_disc_age3_ci <- ci(daly_discounted_sim[ ,3])
daly_disc_age4_ci <- ci(daly_discounted_sim[ ,4])

hosps_deaths_dalys <- cbind(total_hospitalizations_sim, total_deaths_due_to_gnr_sim, total_daly_sim, total_daly_discounted_sim)

#TO OUTPUT PEDIATRIC DATA
# return(data.frame(expected_hosp = expected_hospitalizations, hosp_low_ci = total_hospitalizations_ci[1], hosp_high_ci = total_hospitalizations_ci[2],
# hosp_age1 = expected_hospitalizations_age1, hosp_age1_low_ci = hosp_age1_ci[1], hosp_age1_high_ci = hosp_age1_ci[2],
# hosp_age2 = expected_hospitalizations_age2, hosp_age2_low_ci = hosp_age2_ci[1], hosp_age2_high_ci = hosp_age2_ci[2],
# hosp_age3 = expected_hospitalizations_age3, hosp_age3_low_ci = hosp_age3_ci[1], hosp_age3_high_ci = hosp_age3_ci[2],
# hosp_age4 = expected_hospitalizations_age4, hosp_age4_low_ci = hosp_age4_ci[1], hosp_age4_high_ci = hosp_age4_ci[2],
# gnr_deaths = expected_gnr_deaths, gnr_deaths_low_ci = total_gnr_deaths_ci[1], gnr_deaths_high_ci = total_gnr_deaths_ci[2],
# gnr_deaths_age1 = expected_gnr_deaths_age1, deaths_age1_low_ci = gnrdeaths_age1_ci[1], deaths_age1_high_ci = gnrdeaths_age1_ci[2],
# gnr_deaths_age2 = expected_gnr_deaths_age2, deaths_age2_low_ci = gnrdeaths_age2_ci[1], deaths_age2_high_ci = gnrdeaths_age2_ci[2],
# gnr_deaths_age3 = expected_gnr_deaths_age3, deaths_age3_low_ci = gnrdeaths_age3_ci[1], deaths_age3_high_ci = gnrdeaths_age3_ci[2],
# gnr_deaths_age4 = expected_gnr_deaths_age4, deaths_age4_low_ci = gnrdeaths_age4_ci[1], deaths_age4_high_ci = gnrdeaths_age4_ci[2],
# daly_disc = expected_daly_disc, daly_disc_low_ci = total_daly_discounted_ci[1], daly_disc_high_ci = total_daly_discounted_ci[2],
# daly_disc_age1 = daly_disc_age1, daly_disc_age1_low_ci = daly_disc_age1_ci[1], daly_disc_age1_high_ci = daly_disc_age1_ci[2],
# daly_disc_age2 = daly_disc_age2, daly_disc_age2_low_ci = daly_disc_age2_ci[1], daly_disc_age2_high_ci = daly_disc_age2_ci[2],
# daly_disc_age3 = daly_disc_age3, daly_disc_age3_low_ci = daly_disc_age3_ci[1], daly_disc_age3_high_ci = daly_disc_age3_ci[2],
# daly_disc_age4 = daly_disc_age4, daly_disc_age4_low_ci = daly_disc_age4_ci[1], daly_disc_age4_high_ci = daly_disc_age4_ci[2],
# daly = expected_daly, daly_low_ci = total_daly_ci[1], daly_high_ci = total_daly_ci[2],
# daly_age1 = daly_age1, daly_age1_low_ci = daly_age1_ci[1], daly_age1_high_ci = daly_age1_ci[2],
# daly_age2 = daly_age2, daly_age2_low_ci = daly_age2_ci[1], daly_age2_high_ci = daly_age2_ci[2],
# daly_age3 = daly_age3, daly_age3_low_ci = daly_age3_ci[1], daly_age3_high_ci = daly_age3_ci[2],
# daly_age4 = daly_age4, daly_age4_low_ci = daly_age4_ci[1], daly_age4_high_ci = daly_age4_ci[2]))

#OUTPUT MATRIX OF HOSPITALIZATIONS, DEATHS, AND DALYS LOST
#return(hosps_deaths_dalys)

}

#generate data
#OUTPUT PRIMARY DATA
#sq_data <- simulate_sq(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations)

#OUTPUT MATRIX OF HOSPITALIZATIONS, DEATHS, AND DALYS LOST
sq_hosps_deaths_dalys <- simulate_sq(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations)

#EXPORT SQ DATA
library("writexl")
write_xlsx(sq_data,"C:/Users/wstil/OneDrive/Desktop/Aim 3/Hosp Deaths DALY Figure/sq_output.xlsx", col_name=TRUE, format_headers=FALSE)