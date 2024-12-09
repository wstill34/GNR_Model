# Sets the initial parameters birth cohort (number of individuals at time 0) and 
# periods of time (4 periods: 0-<4 months, 4-<12 months, 12-<24 months, 24-<60 months)
birth_cohort <- 112680
periods <- 4

h_vec <- c(71, 20, 12, 21) # GNR hospitalized patients by age group
h_n <- c(126449, 252898, 502662, 1038045) # n for age group

mu_gnr_vec <- c(41, 7, 7, 8)
mu_gnr_n <- c(71, 19, 12, 19)

# Sets age-specific parameters for each period
h <- h_vec / h_n  # Probability of hospitalization due to GNR infections
mu_ac <- c(0.0348, 0.0348, 0.004725, 0.02835) # All-cause mortality probability
mu_gnr <- mu_gnr_vec / mu_gnr_n    # Probability of death given hospitalization due to GNR infection
length <- c("0 - <4 months ", "4 - <12 months", "12 - <24 months", "24 - 59 months")

# Sets maternal vaccination parameters
vm <- 0.43  # Proportion of mothers vaccinated
vm_se <- 0.043 # Standard error for v_m
em <- 0.75  # Efficacy of maternal vaccine
em_se <- 0.075 # Standard error for e_m

# Defines number of simulations for Monte Carlo model
n_simulations <- 1000

# Creates matrices to store the results from simulations
hospitalizations_total_sim <- matrix(0, nrow=n_simulations, ncol=periods)
hospitalizations_unvax_sim <- matrix(0, nrow=n_simulations, ncol=periods)
hospitalizations_vaxed_sim <- matrix(0, nrow=n_simulations, ncol=periods)
#Added h_vec_sim matrix of 1000 rows, 4 columns
h_vec_sim <- matrix(0, nrow=n_simulations, ncol=periods)
#Also added mu_gnr_sim matrix of 1000 rows, 4 columns
mu_gnr_sim <- matrix(0, nrow=n_simulations, ncol=periods)
unvax_deaths_due_to_gnr_sim <- matrix(0, nrow=n_simulations, ncol=periods)
vaxed_deaths_due_to_gnr_sim <- matrix(0, nrow=n_simulations, ncol=periods)
total_deaths_due_to_gnr_sim <- matrix(0, nrow=n_simulations, ncol=periods)
unvax_other_cause_deaths_sim <- matrix(0, nrow=n_simulations, ncol=periods)
vaxed_other_cause_deaths_sim <- matrix(0, nrow=n_simulations, ncol=periods)
total_other_cause_deaths_sim <- matrix(0, nrow=n_simulations, ncol=periods)
total_deaths_sim <- matrix(0, nrow=n_simulations, ncol=periods)
survived_sim <- matrix(0, nrow=n_simulations, ncol=periods)
survived_vaxed_sim <- matrix(0, nrow=n_simulations, ncol=periods)
survived_unvax_sim <- matrix(0, nrow=n_simulations, ncol=periods)
vm_sample <- matrix(0, nrow=n_simulations, ncol=1)
em_sample <- matrix(0, nrow=n_simulations, ncol=1)

# Monte Carlo simulation
for (sim in 1:n_simulations) {

  # Sample v_m and e_m
  vm_sample[sim] <- rnorm(1, mean=vm, sd=vm_se)
  em_sample[sim] <- rnorm(1, mean=em, sd=em_se)
  
  # Ensure probabilities remain within [0, 1]
  vm_sample[sim] <- max(min(vm_sample[sim], 1), 0)
  em_sample[sim] <- max(min(em_sample[sim], 1), 0)
  
  remaining_cohort <- birth_cohort
  vaxed_population <- round(remaining_cohort*vm_sample[sim])
  unvax_population <- round(remaining_cohort*(1 - vm_sample[sim]))
  
  for (age in 1:periods) {
    if (age == 1) {
      # Period 1: no one is vaccinated before 4 months'
      h_vec_sim[sim, age] <- rbeta(1, h_vec[age], h_n[age]-h_vec[age]) #For each age group, estimates a hospitalization probability based on the mean and standard deviation and 1 trial
      hospitalizations_unvax_sim[sim, age] <- rbinom(1, unvax_population, h_vec_sim[sim,age]) #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
      hospitalizations_vaxed_sim[sim, age] <- rbinom(1, vaxed_population, h_vec_sim[sim, age]*(1 - em_sample[sim]))
      hospitalizations_total_sim[sim, age] <- hospitalizations_unvax_sim[sim, age] + hospitalizations_vaxed_sim[sim, age]
      
      mu_gnr_sim[sim, age] <- rbeta(1, mu_gnr_vec[age], mu_gnr_n[age] - mu_gnr_vec[age])
      unvax_deaths_due_to_gnr_sim[sim, age] <- rbinom(1, hospitalizations_unvax_sim[sim, age], mu_gnr_sim[sim,age])
      vaxed_deaths_due_to_gnr_sim[sim, age] <- rbinom(1, hospitalizations_vaxed_sim[sim, age], mu_gnr_sim[sim,age])
      total_deaths_due_to_gnr_sim[sim, age] <- unvax_deaths_due_to_gnr_sim[sim, age] + vaxed_deaths_due_to_gnr_sim[sim, age]
      unvax_other_cause_deaths_sim[sim, age] <- rbinom(1, unvax_population, mu_ac[age]) - unvax_deaths_due_to_gnr_sim[sim, age]
      vaxed_other_cause_deaths_sim[sim, age] <- rbinom(1, vaxed_population, mu_ac[age]) - vaxed_deaths_due_to_gnr_sim[sim, age]
      total_other_cause_deaths_sim[sim, age] <- unvax_other_cause_deaths_sim[sim, age] + vaxed_other_cause_deaths_sim[sim, age]
      
      # Calculates the number of survivors
      survived_sim[sim, age] <- remaining_cohort - total_deaths_due_to_gnr_sim[sim, age] - total_other_cause_deaths_sim[sim, age]
      survived_vaxed_sim[sim, age] <- vaxed_population - vaxed_deaths_due_to_gnr_sim[sim, age] - vaxed_other_cause_deaths_sim[sim, age]
      survived_unvax_sim[sim, age] <- unvax_population - unvax_deaths_due_to_gnr_sim[sim, age] - unvax_other_cause_deaths_sim[sim, age]
      total_deaths_sim[sim, age] <- total_deaths_due_to_gnr_sim[sim, age] + total_other_cause_deaths_sim[sim, age]
      
      # Update the remaining cohort for the next period
      remaining_cohort <- survived_sim[sim, age]
      
    } else {
      
      h_vec_sim[sim, age] <- rbeta(1, h_vec[age], h_n[age]-h_vec[age]) #For each age group, estimates a hospitalization probability based on the mean and standard deviation and 1 trial
      hospitalizations_total_sim[sim, age] <- rbinom(1, remaining_cohort, h_vec_sim[sim,age]) #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
      mu_gnr_sim[sim, age] <- rbeta(1, mu_gnr_vec[age], mu_gnr_n[age] - mu_gnr_vec[age])
      total_deaths_due_to_gnr_sim[sim, age] <- rbinom(1, hospitalizations_total_sim[sim, age], mu_gnr_sim[sim,age])
      total_other_cause_deaths_sim[sim, age] <- rbinom(1, remaining_cohort, mu_ac[age]) - total_deaths_due_to_gnr_sim[sim, age]

      # Calculates the number of survivors
      survived_sim[sim, age] <- remaining_cohort - total_deaths_due_to_gnr_sim[sim, age] - total_other_cause_deaths_sim[sim, age]
      total_deaths_sim[sim, age] <- total_deaths_due_to_gnr_sim[sim, age] + total_other_cause_deaths_sim[sim, age]
      
      # Update the remaining cohort for the next period
      remaining_cohort <- survived_sim[sim, age]
    }
  }
}

  
# Summarizes the results with confidence intervals
# First adds the total hospitalizations and deaths for each simulation
hospitalizations_sim <- rowSums(hospitalizations_total_sim)
deaths_due_to_gnr_sim <- rowSums(total_deaths_due_to_gnr_sim)
other_cause_deaths_sim <- rowSums(total_other_cause_deaths_sim)
all_cause_deaths_sim <- rowSums(total_deaths_sim)

# Establishes the confidence interval function
ci <- function(x) quantile(x, probs=c(0.025, 0.975))

total_hospitalizations_ci <- ci(hospitalizations_sim)
total_deaths_due_to_gnr_ci <- ci(deaths_due_to_gnr_sim)
total_other_cause_deaths_ci <- ci(other_cause_deaths_sim)
total_all_cause_deaths_ci <- ci(all_cause_deaths_sim)


# Print the results
cat("Total hospitalizations due to GNR infections over 5 years:", mean(hospitalizations_sim), "\n")
cat("95% CI for total hospitalizations:", total_hospitalizations_ci, "\n\n")

cat("Total deaths due to GNR infections over 5 years:", mean(deaths_due_to_gnr_sim), "\n")
cat("95% CI for total deaths due to GNR:", total_deaths_due_to_gnr_ci, "\n\n")

cat("Total deaths due to other causes over 5 years:", mean(other_cause_deaths_sim), "\n")
cat("95% CI for total deaths due to other causes:", total_other_cause_deaths_ci, "\n\n")

cat("Total all-cause deaths over 5 years:", mean(all_cause_deaths_sim), "\n")
cat("95% CI for total all-cause deaths:", total_all_cause_deaths_ci, "\n\n")

cat("Number of individuals surviving at the end of 5 years:", mean(survived_sim[, periods]), "\n")
cat("95% CI for survivors:", ci(survived_sim[, periods]), "\n\n")

# Displays results by year
cat("\nResults by period (mean values):\n")
cat("Period\tHospitalizations\tGNR Deaths\tOther-cause Deaths\tTotal Deaths\tSurvivors\n")
for (age in 1:periods) {
  cat(age, "\t\t", round(mean(hospitalizations_total_sim[, age])), "\t\t", round(mean(total_deaths_due_to_gnr_sim[, age])), "\t\t", mean(total_other_cause_deaths_sim[, age]), "\t\t", mean(total_deaths_sim[, age]), "\t\t", mean(survived_sim[, age]), "\n")
}

