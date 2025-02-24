#Creates a function to simulate the maternal vaccine outputs
simulate_maternal <- function (birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, coverage_mat, efficacy_mat, remaining_years, mat_cost_per_dose) {
  
  sq_hospitalizations <- c(sq_hosps_deaths_dalys[,1])
  sq_hosp_cost <- sq_hospitalizations*157.50
  sq_deaths <- c(sq_hosps_deaths_dalys[,2])
  sq_daly_undisc <- c(sq_hosps_deaths_dalys[,3])
  sq_daly_disc <- c(sq_hosps_deaths_dalys[,4])
  
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
  vaxed_population <- remaining_cohort*coverage_mat
  unvax_population <- remaining_cohort*(1 - coverage_mat)
  
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
      
    } else {
      
      hospitalizations_sim[sim, age] <- (remaining_cohort*h_vec_sim[sim,age])
      deaths_due_to_gnr_sim[sim, age] <- (hospitalizations_sim[sim, age]*mu_gnr_sim[sim,age])
      other_cause_deaths_sim[sim, age] <- (remaining_cohort*mu_ac[age]) - deaths_due_to_gnr_sim[sim, age]
      daly_sim[sim, age] <- deaths_due_to_gnr_sim[sim, age] * remaining_years[age]
      daly_discounted_sim[sim, age] <- deaths_due_to_gnr_sim[sim, age] * discounted_remaining_years[age]

      # Calculates the number of survivors
      survived_sim[sim, age] <- remaining_cohort - deaths_due_to_gnr_sim[sim, age] - other_cause_deaths_sim[sim, age]
      total_deaths_sim[sim, age] <- deaths_due_to_gnr_sim[sim, age] + other_cause_deaths_sim[sim, age]
      
      # Update the remaining cohort for the next period
      remaining_cohort <- survived_sim[sim, age]
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

mat_hosp_cost <- total_hospitalizations_sim*157.50
mat_doses_admin <- birth_cohort*coverage_mat
#mat_cost_per_dose <- 1.47
mat_dose_cost <- mat_cost_per_dose*mat_doses_admin
mat_total_cost <- mat_hosp_cost + mat_dose_cost
incr_cost <- mat_total_cost - sq_hosp_cost
#DISCOUNTED
daly_avert <- sq_daly_disc - total_daly_discounted_sim
#UNDISCOUNTED (USE THIS)
#daly_avert <- sq_daly_undisc - total_daly_sim
icer <- incr_cost / daly_avert

expected_incr_cost <- mean(incr_cost)
expected_daly_avert <- mean(daly_avert)
daly_avert_ci <- ci(daly_avert)
expected_icer <- median(icer)
expected_icer2 <- expected_incr_cost/expected_daly_avert
icer_ci <- ci(icer)
expected_icer2 <- expected_incr_cost/expected_daly_avert

icer_matrix_mat_sq <-cbind(incr_cost, daly_avert, icer, mat_total_cost)

reduced_hosp_n <- sq_hospitalizations - total_hospitalizations_sim
reduced_hosp_percent <- (sq_hospitalizations - total_hospitalizations_sim)/sq_hospitalizations
expected_reduced_hosp <- mean(reduced_hosp_n)
reduced_hosp_ci <- ci(reduced_hosp_n)
expected_reduced_hosp_percent <- mean(reduced_hosp_percent)
reduced_hosp_percent_ci <- ci(reduced_hosp_percent)

reduced_deaths_n <- sq_deaths - total_deaths_due_to_gnr_sim
reduced_deaths_percent <- (sq_deaths - total_deaths_due_to_gnr_sim)/sq_deaths
expected_reduced_deaths <- mean(reduced_deaths_n)
reduced_deaths_ci <- ci(reduced_deaths_n)
expected_reduced_deaths_percent <- mean(reduced_deaths_percent)
reduced_deaths_percent_ci <- ci(reduced_deaths_percent)



#TO OUTPUT PRIMARY MATERNAL DATA
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

#OUTPUT ICER DATA
return(data.frame(daly_avert = expected_daly_avert, daly_low = daly_avert_ci[1], daly_high = daly_avert_ci[2],
icer = expected_icer2, icer_low = icer_ci[1], icer_high = icer_ci[2], icer_median = expected_icer))

# OUTPUT ICER MATRIX
#return(icer_matrix_mat_sq)

#FOR ICER OVER COST ANALYSIS, FIGURE 3
#return(data.frame(incr_cost = expected_incr_cost, daly_avert = expected_daly_avert, icer_median = expected_icer, icer2 = expected_icer2,
#                  icer_low_ci = icer_ci[1], icer_high_ci = icer_ci[2]))

#FOR COVERAGE ANALYSIS, FIGURE 2
#return(data.frame(coverage = coverage_mat, hosp_n = expected_reduced_hosp, hosp_low = reduced_hosp_ci[1], hosp_high = reduced_hosp_ci[2],
#                  hosp_percent = expected_reduced_hosp_percent, hosp_percent_low = reduced_hosp_percent_ci[1], hosp_percent_high = reduced_hosp_percent_ci[2],
#                  death_n = expected_reduced_deaths, deaths_low = reduced_deaths_ci[1], deaths_high = reduced_deaths_ci[2],
#                  death_percent = expected_reduced_deaths_percent, death_percent_low = reduced_deaths_percent_ci[1], death_percent_high = reduced_deaths_percent_ci[2]))

}
#TO OUTPUT PRIMARY MATERNAL DATA
#maternal_data <- simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, coverage_mat, efficacy_mat, remaining_years, mat_cost_per_dose)

#TO OUTPUT MATRIX OF HOSPITALIZATIONS, DEATHS, AND DALYS LOST
#maternal_hosps_deaths_dalys <- simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, coverage_mat, efficacy_mat, remaining_years, mat_cost_per_dose)

#TO OUTPUT ICER VARIABLES
maternal_icer <-simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, coverage_mat, efficacy_mat, remaining_years, mat_cost_per_dose)

#TO OUTPUT ICER MATRIX
#maternal_icer_matrix <- simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, coverage_mat, efficacy_mat, remaining_years, mat_cost_per_dose)

#maternal_icer_check <- simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, coverage_mat, efficacy_mat, remaining_years, mat_cost_per_dose)

#FIGURE 2 MATERNAL COVERAGE
#maternal_coverage_change <- simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, coverage_mat=1, efficacy_mat, remaining_years, mat_cost_per_dose)

#EXPORT MATERNAL DATA
#library("writexl")
#write_xlsx(maternal_data,"C:/Users/wstil/OneDrive/Desktop/Aim 3/Hosp Deaths DALY Figure/maternal_output.xlsx", col_name=TRUE, format_headers=FALSE)

