#Creates a function to simulate the mean and standard deviation for number of hospitalizations and deaths under EPI/childhood vaccine scenario
simulate_epi <- function (birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, v1, efficacy_epi) {

  sq_hospitalizations <- c(status_quo_hosps_deaths_dalys[,1])
  sq_hosp_cost <- sq_hospitalizations*157.50
  sq_deaths <- c(status_quo_hosps_deaths_dalys[,2])
  sq_daly_undisc <- c(status_quo_hosps_deaths_dalys[,3])
  #sq_daly_disc <- c(status_quo_hosps_deaths_dalys[,4])
  mat_daly_undisc <- c(maternal_hosps_deaths_dalys[,3])
  #mat_daly_disc <- c(maternal_hosps_deaths_dalys[,4])
  mat_total_cost <- maternal_total_cost
  
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
  #vaxed_population <- round(remaining_cohort*v1)
  #unvax_population <- round(remaining_cohort*(1-v1))
  for (age in 1:periods) {
    if (age == 1) {
      # Period 1: no one is vaccinated before 4 months'
      #h_vec_sim_epi[sim, age] <- rbeta(1, h_vec[age], h_n[age]-h_vec[age]) #For each age group, estimates a hospitalization probability based on the mean and standard deviation and 1 trial
      hospitalizations_unvax_sim_epi[sim, age] <- remaining_cohort*h_vec_sim[sim,age]
      #hospitalizations_unvax_sim_epi[sim, age] <- rbinom(1, remaining_cohort, h_vec_sim[sim,age]) #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
      hospitalizations_vaxed_sim_epi[sim, age] <- 0
      hospitalizations_total_sim_epi[sim, age] <- hospitalizations_unvax_sim_epi[sim, age] + hospitalizations_vaxed_sim_epi[sim, age]
      
      #mu_gnr_sim_epi[sim, age] <- rbeta(1, mu_gnr_vec[age], mu_gnr_n[age] - mu_gnr_vec[age])
      unvax_deaths_due_to_gnr_sim_epi[sim, age] <- hospitalizations_unvax_sim_epi[sim,age]*mu_gnr_sim[sim,age]
      #unvax_deaths_due_to_gnr_sim_epi[sim, age] <- rbinom(1, hospitalizations_unvax_sim_epi[sim, age], mu_gnr_sim[sim,age])
      vaxed_deaths_due_to_gnr_sim_epi[sim, age] <- 0
      total_deaths_due_to_gnr_sim_epi[sim, age] <- unvax_deaths_due_to_gnr_sim_epi[sim, age] + vaxed_deaths_due_to_gnr_sim_epi[sim, age]
      unvax_other_cause_deaths_sim_epi[sim, age] <- remaining_cohort*mu_ac[age] - unvax_deaths_due_to_gnr_sim_epi[sim, age]
      #unvax_other_cause_deaths_sim_epi[sim, age] <- rbinom(1, remaining_cohort, mu_ac[age]) - unvax_deaths_due_to_gnr_sim_epi[sim, age]
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
      vaxed_population <- survived_vaxed_sim_epi[sim, age]
      unvax_population <- survived_unvax_sim_epi[sim, age]
    } else {
      
      #h_vec_sim_epi[sim, age] <- rbeta(1, h_vec[age], h_n[age]-h_vec[age]) #For each age group, estimates a hospitalization probability based on the mean and standard deviation and 1 trial
      hospitalizations_unvax_sim_epi[sim, age] <- unvax_population*h_vec_sim[sim,age] #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
      hospitalizations_vaxed_sim_epi[sim, age] <- vaxed_population*h_vec_sim[sim,age]*(1 - efficacy_epi[sim])
      #hospitalizations_unvax_sim_epi[sim, age] <- rbinom(1, unvax_population, h_vec_sim[sim,age]) #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
      #hospitalizations_vaxed_sim_epi[sim, age] <- rbinom(1, vaxed_population, h_vec_sim[sim,age]*(1 - efficacy_epi[sim])) #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
      hospitalizations_total_sim_epi[sim, age] <- hospitalizations_unvax_sim_epi[sim, age] + hospitalizations_vaxed_sim_epi[sim, age]
      
      #mu_gnr_sim_epi[sim, age] <- rbeta(1, mu_gnr_vec[age], mu_gnr_n[age] - mu_gnr_vec[age])
      unvax_deaths_due_to_gnr_sim_epi[sim, age] <- hospitalizations_unvax_sim_epi[sim, age]*mu_gnr_sim[sim,age]
      vaxed_deaths_due_to_gnr_sim_epi[sim, age] <- hospitalizations_vaxed_sim_epi[sim, age]*mu_gnr_sim[sim,age]
      #unvax_deaths_due_to_gnr_sim_epi[sim, age] <- rbinom(1, hospitalizations_unvax_sim_epi[sim, age], mu_gnr_sim[sim,age])
      #vaxed_deaths_due_to_gnr_sim_epi[sim, age] <- rbinom(1, hospitalizations_vaxed_sim_epi[sim, age], mu_gnr_sim[sim,age])
      total_deaths_due_to_gnr_sim_epi[sim, age] <- unvax_deaths_due_to_gnr_sim_epi[sim, age] + vaxed_deaths_due_to_gnr_sim_epi[sim, age]
      unvax_other_cause_deaths_sim_epi[sim, age] <- unvax_population*mu_ac[age] - unvax_deaths_due_to_gnr_sim_epi[sim, age]
      vaxed_other_cause_deaths_sim_epi[sim, age] <- vaxed_population*mu_ac[age] - vaxed_deaths_due_to_gnr_sim_epi[sim, age]
      #unvax_other_cause_deaths_sim_epi[sim, age] <- rbinom(1, unvax_population, mu_ac[age]) - unvax_deaths_due_to_gnr_sim_epi[sim, age]
      #vaxed_other_cause_deaths_sim_epi[sim, age] <- rbinom(1, vaxed_population, mu_ac[age]) - vaxed_deaths_due_to_gnr_sim_epi[sim, age]
      total_other_cause_deaths_sim_epi[sim, age] <- unvax_other_cause_deaths_sim_epi[sim, age] + vaxed_other_cause_deaths_sim_epi[sim, age]
      
      # Calculates the number of survivors
      survived_unvax_sim_epi[sim, age] <- unvax_population - unvax_deaths_due_to_gnr_sim_epi[sim, age] - unvax_other_cause_deaths_sim_epi[sim, age]
      survived_vaxed_sim_epi[sim, age] <- vaxed_population - vaxed_deaths_due_to_gnr_sim_epi[sim, age] - vaxed_other_cause_deaths_sim_epi[sim, age]
      survived_sim_epi[sim, age] <- survived_unvax_sim_epi[sim, age] + survived_vaxed_sim_epi[sim, age]
      total_deaths_sim_epi[sim, age] <- total_deaths_due_to_gnr_sim_epi[sim, age] + total_other_cause_deaths_sim_epi[sim, age]
      daly_sim[sim, age] <- total_deaths_due_to_gnr_sim_epi[sim, age] * remaining_years[age]
      daly_discounted_sim[sim, age] <- total_deaths_due_to_gnr_sim_epi[sim, age] * discounted_remaining_years[age]
      
      # Update the remaining cohort for the next period
      unvax_population <- survived_unvax_sim_epi[sim, age]
      vaxed_population <- survived_vaxed_sim_epi[sim, age]
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

hosps_deaths_dalys <- cbind(hospitalizations_sim_epi, deaths_due_to_gnr_sim_epi, total_daly_sim, total_daly_discounted_sim)

epi_hosp_cost <- hospitalizations_sim_epi*157.50
epi_cost_per_dose <- 1
epi_doses_admin <- birth_cohort*(1-mu_ac[1])*v1*3 #coverage applied to those who survive period 1
epi_dose_cost <- epi_cost_per_dose*epi_doses_admin
epi_total_cost <- epi_hosp_cost + epi_dose_cost
incr_cost <- epi_total_cost - sq_hosp_cost
incr_cost_epi_mat <- epi_total_cost - mat_total_cost
#daly_avert_epi_sq <- sq_daly_disc - total_daly_discounted_sim
daly_avert_epi_sq <- sq_daly_undisc - total_daly_sim
daly_avert_epi_mat <- mat_daly_undisc - total_daly_sim
icer_epi_sq <- incr_cost / daly_avert_epi_sq
icer_epi_mat <- incr_cost_epi_mat / daly_avert_epi_mat

expected_daly_avert <- mean(daly_avert_epi_sq)
#expected_daly_avert <- median(daly_avert_epi_mat)
daly_avert_ci <- ci(daly_avert_epi_sq)
#daly_avert_ci <- ci(daly_avert_epi_mat)
expected_icer <- median(icer_epi_sq)
#expected_icer <- median(icer_epi_mat)
icer_ci <- ci(icer_epi_sq)
#icer_ci <- ci(icer_epi_mat)


#TO OUTPUT PEDIATRIC DATA
# return(data.frame(expected_hosp = expected_hospitalizations_epi, hosp_low_ci = total_hospitalizations_ci[1], hosp_high_ci = total_hospitalizations_ci[2],
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
#return(data.frame(daly_avert = expected_daly_avert, daly_low = daly_avert_ci[1], daly_high = daly_avert_ci[2],
                 # icer = expected_icer, icer_low = icer_ci[1], icer_high = icer_ci[2]))

#TO OUTPUT DALY LOST MATRIX
#return(total_daly_discounted_sim)

#TO OUTPUT INCREMENTAL COST MATRIX
#return(incr_cost)

#TO OUTPUT DALY AVERT MATRIX
return(daly_avert_epi_sq)

#TO OUTPUT ICER MATRIX
#return(icer_epi_sq)
}

#TO OUTPUT PEDIATRIC DATA
#epi_data_final <- simulate_epi(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, v1, efficacy_epi)

#TO OUTPUT MATRIX OF HOSPITALIZATIONS, DEATHS, AND DALYS LOST
#epi_hosps_deaths_dalys <- simulate_epi(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, v1, efficacy_epi)

#TO OUTPUT ICER VARIABLES
#epi_icer_sq <- simulate_epi(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, v1, efficacy_epi)

#TO OUTPUT DALY LOST MATRIX
#epi_daly_lost <- simulate_epi(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, v1, efficacy_epi)

#TO OUTPUT INCR COST MATRIX
#epi_incr_cost <- simulate_epi(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, v1, efficacy_epi)

#TO OUTPUT DALY AVERT MATRIX
epi_daly_avert_matrix <- simulate_epi(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, v1, efficacy_epi)
  
#EXPORT PEDIATRIC DATA
#library("writexl")
#write_xlsx(epi_data_final,"C:/Users/wstil/OneDrive/Desktop/Aim 3/Hosp Deaths DALY Figure/epi_output.xlsx", col_name=TRUE, format_headers=FALSE)