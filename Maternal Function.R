#Creates a function to simulate the mean and standard deviation for number of hospitalizations and deaths under maternal vaccine scenario
simulate_maternal <- function (birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, efficacy_mat, remaining_years, mat_cost_per_dose) {
  
  sq_hospitalizations <- c(status_quo_hosps_deaths_dalys[,1])
  sq_hosp_cost <- sq_hospitalizations*157.50
  sq_deaths <- c(status_quo_hosps_deaths_dalys[,2])
  sq_daly_undisc <- c(status_quo_hosps_deaths_dalys[,3])
  sq_daly_disc <- c(status_quo_hosps_deaths_dalys[,4])
  
  

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
  vaxed_population <- (remaining_cohort*vm)
  unvax_population <- (remaining_cohort*(1 - vm))
  
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
      unvax_other_cause_deaths_sim_mat[sim, age] <- (unvax_population*mu_ac[age]) - unvax_deaths_due_to_gnr_sim_mat[sim, age]
      # removing for now unvax_other_cause_deaths_sim_mat[sim, age] <- rbinom(1, unvax_population, mu_ac[age]) - unvax_deaths_due_to_gnr_sim_mat[sim, age]
      vaxed_other_cause_deaths_sim_mat[sim, age] <- (vaxed_population*mu_ac[age]) - vaxed_deaths_due_to_gnr_sim_mat[sim, age]
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
      total_other_cause_deaths_sim_mat[sim, age] <- (remaining_cohort*mu_ac[age]) - total_deaths_due_to_gnr_sim_mat[sim, age]
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

hosps_deaths_dalys <- cbind(hospitalizations_sim_mat, deaths_due_to_gnr_sim_mat, total_daly_sim, total_daly_discounted_sim)

mat_hosp_cost <- hospitalizations_sim_mat*157.50
mat_doses_admin <- birth_cohort*vm
mat_dose_cost <- mat_cost_per_dose*mat_doses_admin
mat_total_cost <- mat_hosp_cost + mat_dose_cost
incr_cost <- mat_total_cost - sq_hosp_cost
#DISCOUNTED
daly_avert <- sq_daly_disc - total_daly_discounted_sim
#UNDISCOUNTED (USE THIS)
#daly_avert <- sq_daly_undisc - total_daly_sim
icer <- incr_cost / daly_avert

expected_daly_avert <- median(daly_avert)
daly_avert_ci <- ci(daly_avert)
expected_icer <- median(icer)
icer_ci <- ci(icer)


#TO OUTPUT PRIMARY MATERNAL DATA
# return(data.frame(expected_hosp = expected_hospitalizations_mat, hosp_low_ci = total_hospitalizations_ci[1], hosp_high_ci = total_hospitalizations_ci[2],
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

#TO OUTPUT ICER VARIABLES
return(data.frame(daly_avert = expected_daly_avert, daly_low = daly_avert_ci[1], daly_high = daly_avert_ci[2],
icer = expected_icer, icer_low = icer_ci[1], icer_high = icer_ci[2]))

#TO OUTPUT DALY LOST MATRIX
#return(total_daly_sim)

#TO OUTPUT INCREMENTAL COST MATRIX
#return(incr_cost)

#TO OUTPUT DALY AVERT MATRIX
#return(daly_avert)

#TO OUTPUT MATERNAL COST MATRIX
#return(mat_total_cost)

#TO OUTPUT ICER MATRIX
#return(icer)

}
#TO OUTPUT PRIMARY MATERNAL DATA
#maternal_data_final <- simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, efficacy_mat, remaining_years, mat_cost_per_dose)

#TO OUTPUT MATRIX OF HOSPITALIZATIONS, DEATHS, AND DALYS LOST
#maternal_hosps_deaths_dalys <- simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, efficacy_mat, remaining_years, mat_cost_per_dose)

#TO OUTPUT ICER VARIABLES
maternal_icer_new <-simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, efficacy_mat, remaining_years, mat_cost_per_dose)

#TO OUTPUT DALY LOST MATRIX
#maternal_daly_lost <- simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, efficacy_mat, remaining_years, mat_cost_per_dose)

#TO OUTPUT INCR COST MATRIX
#maternal_incr_cost <- simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, efficacy_mat, remaining_years, mat_cost_per_dose)

#TO OUTPUT DALY AVERT MATRIX
#maternal_daly_avert <- simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, efficacy_mat, remaining_years, mat_cost_per_dose)

#TO OUTPUT ICER MATRIX
#maternal_icer_matrix <- simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, efficacy_mat, remaining_years, mat_cost_per_dose)

#TO OUTPUT MATERNAL COST MATRIX
#maternal_total_cost <- simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, efficacy_mat, remaining_years, mat_cost_per_dose)

#EXPORT MATERNAL DATA
#library("writexl")
#write_xlsx(maternal_data_final,"C:/Users/wstil/OneDrive/Desktop/Aim 3/Hosp Deaths DALY Figure/maternal_output.xlsx", col_name=TRUE, format_headers=FALSE)