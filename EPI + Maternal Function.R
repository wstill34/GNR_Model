#Creates a function to simulate the joint epi/maternal vaccine scenario outputs
simulate_joint <- function (birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, coverage_mat, coverage_epi, efficacy_mat, efficacy_epi) {

  sq_hospitalizations <- c(sq_hosps_deaths_dalys[,1])
  sq_hosp_cost <- sq_hospitalizations*157.50
  sq_deaths <- c(sq_hosps_deaths_dalys[,2])
  sq_daly_undisc <- c(sq_hosps_deaths_dalys[,3])
  sq_daly_disc <- c(sq_hosps_deaths_dalys[,4])
  mat_daly_undisc <- c(maternal_hosps_deaths_dalys[,3])
  mat_daly_disc <- c(maternal_hosps_deaths_dalys[,4])
  mat_total_cost <- c(maternal_icer_matrix[,4])
  
  # Creates matrices to store the results from simulations
  hospitalizations_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  hospitalizations_unvax_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  hospitalizations_vaxed_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  unvax_deaths_due_to_gnr_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  vaxed_deaths_due_to_gnr_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  deaths_due_to_gnr_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  unvax_other_cause_deaths_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  vaxed_other_cause_deaths_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  other_cause_deaths_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  total_deaths_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  survived_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  survived_vaxed_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  survived_unvax_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  daly_sim <- matrix(0, nrow=n_simulations, ncol=periods)
  daly_discounted_sim <- matrix(0, nrow = n_simulations, ncol=periods)
  

  # Monte Carlo simulation
  for (sim in 1:n_simulations) {
    
    remaining_cohort <- birth_cohort
    vaxed_population <- remaining_cohort * coverage_mat
    unvax_population <- remaining_cohort * (1 - coverage_mat)
    
    for (age in 1:periods) {
      if (age == 1) {
        
        # Period 1: maternal vaccine confers protection to those who received maternal vaccine
        hospitalizations_unvax_sim[sim, age] <- (unvax_population*h_vec_sim[sim,age])
        hospitalizations_vaxed_sim[sim, age] <- (vaxed_population*h_vec_sim[sim,age]*(1 - efficacy_mat[sim]))
        hospitalizations_sim[sim, age] <- hospitalizations_unvax_sim[sim, age] + hospitalizations_vaxed_sim[sim, age]
        
        unvax_deaths_due_to_gnr_sim[sim, age] <- (hospitalizations_unvax_sim[sim, age]*mu_gnr_sim[sim,age])
        vaxed_deaths_due_to_gnr_sim[sim, age] <- (hospitalizations_vaxed_sim[sim, age]*mu_gnr_sim[sim,age])
        deaths_due_to_gnr_sim[sim, age] <- unvax_deaths_due_to_gnr_sim[sim, age] + vaxed_deaths_due_to_gnr_sim[sim, age]
        unvax_other_cause_deaths_sim[sim, age] <- (unvax_population*mu_ac[age]) - unvax_deaths_due_to_gnr_sim[sim, age]
        vaxed_other_cause_deaths_sim[sim, age] <- (vaxed_population*mu_ac[age]) - vaxed_deaths_due_to_gnr_sim[sim, age]
        other_cause_deaths_sim[sim, age] <- unvax_other_cause_deaths_sim[sim, age] + vaxed_other_cause_deaths_sim[sim, age]
        
        # Calculates the number of survivors
        survived_sim[sim, age] <- remaining_cohort - deaths_due_to_gnr_sim[sim, age] - other_cause_deaths_sim[sim, age]
        survived_vaxed_sim[sim, age] <- vaxed_population - vaxed_deaths_due_to_gnr_sim[sim, age] - vaxed_other_cause_deaths_sim[sim, age]
        survived_unvax_sim[sim, age] <- unvax_population - unvax_deaths_due_to_gnr_sim[sim, age] - unvax_other_cause_deaths_sim[sim, age]
        total_deaths_sim[sim, age] <- deaths_due_to_gnr_sim[sim, age] + other_cause_deaths_sim[sim, age]
        daly_sim[sim, age] <- deaths_due_to_gnr_sim[sim, age] * remaining_years[age]
        daly_discounted_sim[sim, age] <- deaths_due_to_gnr_sim[sim, age] * discounted_remaining_years[age]
        
        # Update the remaining cohort for the next period
        remaining_cohort <- survived_sim[sim, age]
        vaxed_population <- remaining_cohort * coverage_epi
        unvax_population <- remaining_cohort * (1-coverage_epi)
        
      } else {
        
        hospitalizations_unvax_sim[sim, age] <- unvax_population*h_vec_sim[sim,age] #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
        hospitalizations_vaxed_sim[sim, age] <- vaxed_population*h_vec_sim[sim,age]*(1 - efficacy_epi[sim])
        hospitalizations_sim[sim, age] <- hospitalizations_unvax_sim[sim, age] + hospitalizations_vaxed_sim[sim, age]
        
        unvax_deaths_due_to_gnr_sim[sim, age] <- hospitalizations_unvax_sim[sim, age]*mu_gnr_sim[sim,age]
        vaxed_deaths_due_to_gnr_sim[sim, age] <- hospitalizations_vaxed_sim[sim, age]*mu_gnr_sim[sim,age]
        deaths_due_to_gnr_sim[sim, age] <- unvax_deaths_due_to_gnr_sim[sim, age] + vaxed_deaths_due_to_gnr_sim[sim, age]
        unvax_other_cause_deaths_sim[sim, age] <- unvax_population*mu_ac[age] - unvax_deaths_due_to_gnr_sim[sim, age]
        vaxed_other_cause_deaths_sim[sim, age] <- vaxed_population*mu_ac[age] - vaxed_deaths_due_to_gnr_sim[sim, age]
        other_cause_deaths_sim[sim, age] <- unvax_other_cause_deaths_sim[sim, age] + vaxed_other_cause_deaths_sim[sim, age]
        
        # Calculates the number of survivors
        survived_unvax_sim[sim, age] <- unvax_population - unvax_deaths_due_to_gnr_sim[sim, age] - unvax_other_cause_deaths_sim[sim, age]
        survived_vaxed_sim[sim, age] <- vaxed_population - vaxed_deaths_due_to_gnr_sim[sim, age] - vaxed_other_cause_deaths_sim[sim, age]
        survived_sim[sim, age] <- survived_unvax_sim[sim, age] + survived_vaxed_sim[sim, age]
        total_deaths_sim[sim, age] <- deaths_due_to_gnr_sim[sim, age] + other_cause_deaths_sim[sim, age]
        daly_sim[sim, age] <- deaths_due_to_gnr_sim[sim, age] * remaining_years[age]
        daly_discounted_sim[sim, age] <- deaths_due_to_gnr_sim[sim, age] * discounted_remaining_years[age]
        
        # Update the remaining cohort for the next period
        unvax_population <- survived_unvax_sim[sim, age]
        vaxed_population <- survived_vaxed_sim[sim, age]
        
      }
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
  
joint_hosp_cost <- total_hospitalizations_sim*157.50
mat_cost_per_dose <- 1.47
mat_doses_admin <- birth_cohort*coverage_mat
mat_dose_cost <- mat_cost_per_dose*mat_doses_admin
epi_cost_per_dose <- 1.47
epi_doses_per_child <- 3
epi_doses_admin <- birth_cohort*(1-mu_ac[1])*coverage_epi*epi_doses_per_child #coverage applied to those who survive period 1
epi_dose_cost <- epi_cost_per_dose*epi_doses_admin
joint_dose_cost <- mat_dose_cost + epi_dose_cost
joint_total_cost <- joint_hosp_cost + joint_dose_cost

incr_cost_joint_sq <- joint_total_cost - sq_hosp_cost
incr_cost_joint_mat <- joint_total_cost - mat_total_cost
#daly_avert_joint_sq <- sq_daly_undisc - total_daly_sim
daly_avert_joint_sq <- sq_daly_disc - total_daly_discounted_sim
#daly_avert_joint_mat <- mat_daly_undisc - total_daly_sim
daly_avert_joint_mat <- mat_daly_disc - total_daly_discounted_sim
icer_joint_sq <- incr_cost_joint_sq / daly_avert_joint_sq
icer_joint_mat <- incr_cost_joint_mat / daly_avert_joint_mat

#Comparing Joint to Status Quo
expected_daly_avert_joint_sq <- mean(daly_avert_joint_sq)
daly_avert_joint_sq_ci <- ci(daly_avert_joint_sq)
expected_icer_joint_sq <- mean(incr_cost_joint_sq)/expected_daly_avert_joint_sq
icer_joint_sq_ci <- ci(icer_joint_sq)
#Comparing Joint to Maternal
expected_daly_avert_joint_mat <- mean(daly_avert_joint_mat)
daly_avert_joint_mat_ci <- ci(daly_avert_joint_mat)
expected_icer_joint_mat <- mean(incr_cost_joint_mat)/expected_daly_avert_joint_mat
icer_joint_mat_ci <- ci(icer_joint_mat)

icer_matrix_joint_sq <-cbind(incr_cost_joint_sq, daly_avert_joint_sq, icer_joint_sq)
icer_matrix_joint_mat <- cbind(incr_cost_joint_mat, daly_avert_joint_mat, icer_joint_mat)
 
  #TO OUTPUT JOINT PRIMARY DATA
  # return(data.frame(expected_hosp = expected_hospitalizations, hosp_low_ci = total_hospitalizations_ci[1], hosp_high_ci = total_hospitalizations_ci[2],
  #                   hosp_age1 = expected_hospitalizations_age1, hosp_age1_low_ci = hosp_age1_ci[1], hosp_age1_high_ci = hosp_age1_ci[2],
  #                   hosp_age2 = expected_hospitalizations_age2, hosp_age2_low_ci = hosp_age2_ci[1], hosp_age2_high_ci = hosp_age2_ci[2],
  #                   hosp_age3 = expected_hospitalizations_age3, hosp_age3_low_ci = hosp_age3_ci[1], hosp_age3_high_ci = hosp_age3_ci[2],
  #                   hosp_age4 = expected_hospitalizations_age4, hosp_age4_low_ci = hosp_age4_ci[1], hosp_age4_high_ci = hosp_age4_ci[2],
  #                   gnr_deaths = expected_gnr_deaths, gnr_deaths_low_ci = total_gnr_deaths_ci[1], gnr_deaths_high_ci = total_gnr_deaths_ci[2],
  #                   gnr_deaths_age1 = expected_gnr_deaths_age1, deaths_age1_low_ci = gnrdeaths_age1_ci[1], deaths_age1_high_ci = gnrdeaths_age1_ci[2],
  #                   gnr_deaths_age2 = expected_gnr_deaths_age2, deaths_age2_low_ci = gnrdeaths_age2_ci[1], deaths_age2_high_ci = gnrdeaths_age2_ci[2],
  #                   gnr_deaths_age3 = expected_gnr_deaths_age3, deaths_age3_low_ci = gnrdeaths_age3_ci[1], deaths_age3_high_ci = gnrdeaths_age3_ci[2],
  #                   gnr_deaths_age4 = expected_gnr_deaths_age4, deaths_age4_low_ci = gnrdeaths_age4_ci[1], deaths_age4_high_ci = gnrdeaths_age4_ci[2],
  #                   daly_disc = expected_daly_disc, daly_disc_low_ci = total_daly_discounted_ci[1], daly_disc_high_ci = total_daly_discounted_ci[2],
  #                   daly_disc_age1 = daly_disc_age1, daly_disc_age1_low_ci = daly_disc_age1_ci[1], daly_disc_age1_high_ci = daly_disc_age1_ci[2],
  #                   daly_disc_age2 = daly_disc_age2, daly_disc_age2_low_ci = daly_disc_age2_ci[1], daly_disc_age2_high_ci = daly_disc_age2_ci[2],
  #                   daly_disc_age3 = daly_disc_age3, daly_disc_age3_low_ci = daly_disc_age3_ci[1], daly_disc_age3_high_ci = daly_disc_age3_ci[2],
  #                   daly_disc_age4 = daly_disc_age4, daly_disc_age4_low_ci = daly_disc_age4_ci[1], daly_disc_age4_high_ci = daly_disc_age4_ci[2],
  #                   daly = expected_daly, daly_low_ci = total_daly_ci[1], daly_high_ci = total_daly_ci[2],
  #                   daly_age1 = daly_age1, daly_age1_low_ci = daly_age1_ci[1], daly_age1_high_ci = daly_age1_ci[2],
  #                   daly_age2 = daly_age2, daly_age2_low_ci = daly_age2_ci[1], daly_age2_high_ci = daly_age2_ci[2],
  #                   daly_age3 = daly_age3, daly_age3_low_ci = daly_age3_ci[1], daly_age3_high_ci = daly_age3_ci[2],
  #                   daly_age4 = daly_age4, daly_age4_low_ci = daly_age4_ci[1], daly_age4_high_ci = daly_age4_ci[2]))

  #OUTPUT MATRIX OF HOSPITALIZATIONS, DEATHS, AND DALYS LOST
  #return(hosps_deaths_dalys)

  # OUTPUT JOINT ICER DATA
  return(data.frame(daly_avert_joint_sq = expected_daly_avert_joint_sq, daly_joint_sq_low = daly_avert_joint_sq_ci[1], daly_joint_sq_high = daly_avert_joint_sq_ci[2],
                    icer_joint_sq = expected_icer_joint_sq, icer_joint_sq_low = icer_joint_sq_ci[1], icer_joint_sq_high = icer_joint_sq_ci[2],
                    daly_avert_joint_mat = expected_daly_avert_joint_mat, daly_joint_mat_low = daly_avert_joint_mat_ci[1], daly_joint_mat_high = daly_avert_joint_mat_ci[2],
                    icer_joint_mat = expected_icer_joint_mat, icer_joint_mat_low = icer_joint_mat_ci[1], icer_joint_mat_high = icer_joint_mat_ci[2]))

  # OUTPUT ICER MATRIX - JOINT VS SQ
  #return(icer_matrix_joint_sq)

  # OUTPUT ICER MATRIX - JOINT VS MAT
  #return(icer_matrix_joint_mat)
}

#TO OUTPUT JOINT PRIMARY DATA
#joint_data <- simulate_joint(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, coverage_mat, coverage_epi, efficacy_mat, efficacy_epi)

#TO OUTPUT MATRIX OF HOSPITALIZATIONS, DEATHS, AND DALYS LOST
#joint_hosps_deaths_dalys <- simulate_joint(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, coverage_mat, coverage_epi, efficacy_mat, efficacy_epi)

#TO OUTPUT ICER VARIABLES
joint_icer_data <- simulate_joint(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, coverage_mat, coverage_epi, efficacy_mat, efficacy_epi)

#TO OUTPUT ICER MATRIX JOINT VS SQ
#joint_sq_icer_matrix <- simulate_joint(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, coverage_mat, coverage_epi, efficacy_mat, efficacy_epi)
  
#TO OUTPUT ICER MATRIX JOINT VS MAT
#joint_mat_icer_matrix <- simulate_joint(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, coverage_mat, coverage_epi, efficacy_mat, efficacy_epi)


#EXPORT JOINT DATA
#library("writexl")
#write_xlsx(joint_data,"C:/Users/wstil/OneDrive/Desktop/Aim 3/Hosp Deaths DALY Figure/joint_output.xlsx", col_name=TRUE, format_headers=FALSE)
