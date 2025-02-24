#ICER function for range of dose cost, Figure 3 of paper
#Maternal ICER function
simulate_ICER_mat <- function(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, coverage_mat, efficacy_mat, remaining_years, mat_cost_per_dose, min_cost, max_cost, step_cost) {
  
  min_cost <- 0.25
  max_cost <- 5
  step_cost <- 0.25
  cost_levels <- seq(min_cost, max_cost, by = step_cost)
  
  icer <- numeric(length(cost_levels))
  icer_low_ci <- numeric(length(cost_levels))
  icer_high_ci <- numeric(length(cost_levels))
  icer_median <- numeric(length(cost_levels))
  
    for (i in seq_along(cost_levels)) {
      mat_cost_per_dose <- cost_levels[i]
      icer[i] <- simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, coverage_mat, efficacy_mat, remaining_years, mat_cost_per_dose)[,4]
      icer_low_ci[i] <- simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, coverage_mat, efficacy_mat, remaining_years, mat_cost_per_dose)[,5]
      icer_high_ci[i] <- simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, coverage_mat, efficacy_mat, remaining_years, mat_cost_per_dose)[,6]
      icer_median[i] <- simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, coverage_mat, efficacy_mat, remaining_years, mat_cost_per_dose)[,7]
    }
  
  return(data.frame(cost = cost_levels, icer = icer, icer_low = icer_low_ci, icer_high = icer_high_ci, icer_median = icer_median))
  
}

icer_mat_cost <- simulate_ICER_mat(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, coverage_mat, efficacy_mat, remaining_years, mat_cost_per_dose, min_cost, max_cost, step_cost)


library("writexl")
write_xlsx(icer_mat_cost,"C:/Users/wstil/OneDrive/Desktop/Aim 3/ICER Figure/icer_mat_undisc.xlsx", col_name=TRUE, format_headers=FALSE)

#EPI ICER function

simulate_ICER_epi <- function(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, v1, efficacy_epi, remaining_years, min_cost, max_cost, step_cost) {
  
  baseline_DALY <- status_quo_data[,31]
  epi_DALY <- epi_data[,31]
  
  min_cost <- 0.25
  max_cost <- 5
  step_cost <- 0.25
  cost_levels <- seq(min_cost, max_cost, by = step_cost)
  
  sq_total_cost <- status_quo_data[,1]*157.50
  
  epi_dose_cost <- cost_levels*birth_cohort*(1-mu_ac[1])*v1*3
  epi_hosp_cost <- epi_data[,1]*157.50
  epi_total_cost <- epi_dose_cost + epi_hosp_cost
  epi_incr_cost <- epi_total_cost - sq_total_cost
  
  icer <- numeric(length(cost_levels))
  
  
  for (i in seq_along(cost_levels)) {
    icer[i] <- (epi_incr_cost[i]/(baseline_DALY - epi_DALY))
  }
  
  return(data.frame(cost = cost_levels, icer = icer))
  
}

icer_epi <- simulate_ICER_epi(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, v1, efficacy_epi, remaining_years, min_cost, max_cost, step_cost)


library("writexl")
write_xlsx(icer_epi,"C:/Users/wstil/OneDrive/Desktop/Aim 3/GNR Model/icer_epi.xlsx", col_name=TRUE, format_headers=FALSE)

#Joint ICER function

simulate_ICER_joint <- function(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, efficacy_mat, v1, efficacy_epi, remaining_years, min_cost, max_cost, step_cost) {
  
  baseline_DALY <- status_quo_data[,31]
  joint_DALY <- joint_data[,31]
  
  min_cost <- 0.25
  max_cost <- 5
  step_cost <- 0.25
  cost_levels <- seq(min_cost, max_cost, by = step_cost)
  
  sq_total_cost <- status_quo_data[,1]*157.50
  
  joint_dose_cost <- cost_levels*birth_cohort*vm + cost_levels*birth_cohort*(1-mu_ac[1])*v1*3
  joint_hosp_cost <- joint_data[,1]*157.50
  joint_total_cost <- joint_dose_cost + joint_hosp_cost
  joint_incr_cost <- joint_total_cost - sq_total_cost
  
  icer <- numeric(length(cost_levels))
  
  
  for (i in seq_along(cost_levels)) {
    icer[i] <- (joint_incr_cost[i]/(baseline_DALY - joint_DALY))
  }
  
  return(data.frame(cost = cost_levels, icer = icer))
  
}

icer_joint <- simulate_ICER_joint(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, vm, efficacy_mat, v1, efficacy_epi, remaining_years, min_cost, max_cost, step_cost)

library("writexl")
write_xlsx(icer_joint,"C:/Users/wstil/OneDrive/Desktop/Aim 3/GNR Model/icer_joint.xlsx", col_name=TRUE, format_headers=FALSE)
