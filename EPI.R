# Sets the initial parameters birth cohort (number of individuals at time 0)
# periods of time (4 periods: 0-<4 months, 4-<12 months, 12-<24 months, 24-<60 months )
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

# Sets vaccination parameters
v1 <- 0.8   # Proportion vaccinated at end of period 1
v1_se <- 0.08 # Standard error for v1, for now 10% of v1
e1 <- 0.7   # Vaccine efficacy in reducing hospitalizations
e1_se <- 0.07  # Standard error for e1, for now 10% of e1

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
v1_sample <- matrix(0, nrow=n_simulations, ncol=1)
e1_sample <- matrix(0, nrow=n_simulations, ncol=1)
#remaining_cohort <- matrix(0, nrow=n_simulations, ncol=periods)
#remaining_cohort_end <- matrix(0, nrow=n_simulations, ncol=periods)
#unvax_population <-matrix(0, nrow=n_simulations, ncol=periods)
#vaxed_population <- matrix(0, nrow=n_simulations, ncol=periods)



# Monte Carlo simulation
for (sim in 1:n_simulations) {

  # Sample v1 and e1
  v1_sample[sim] <- rnorm(1, mean=v1, sd=v1_se)
  e1_sample[sim] <- rnorm(1, mean=e1, sd=e1_se)
  
  # Ensure probabilities remain within [0, 1]
  v1_sample[sim] <- max(min(v1_sample[sim], 1), 0)
  e1_sample[sim] <- max(min(e1_sample[sim], 1), 0)
  
  remaining_cohort <- birth_cohort
  vaxed_population <- remaining_cohort*v1_sample[sim]
  unvax_population <- remaining_cohort*(1-v1_sample[sim])
  for (age in 1:periods) {
    if (age == 1) {
      # Period 1: no one is vaccinated before 4 months'
      h_vec_sim[sim, age] <- rbeta(1, h_vec[age], h_n[age]-h_vec[age]) #For each age group, estimates a hospitalization probability based on the mean and standard deviation and 1 trial
      hospitalizations_unvax_sim[sim, age] <- rbinom(1, remaining_cohort, h_vec_sim[sim,age]) #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
      hospitalizations_vaxed_sim[sim, age] <- 0
      hospitalizations_total_sim[sim, age] <- hospitalizations_unvax_sim[sim, age] + hospitalizations_vaxed_sim[sim, age]
      
      mu_gnr_sim[sim, age] <- rbeta(1, mu_gnr_vec[age], mu_gnr_n[age] - mu_gnr_vec[age])
      unvax_deaths_due_to_gnr_sim[sim, age] <- rbinom(1, hospitalizations_unvax_sim[sim, age], mu_gnr_sim[sim,age])
      vaxed_deaths_due_to_gnr_sim[sim, age] <- 0
      total_deaths_due_to_gnr_sim[sim, age] <- unvax_deaths_due_to_gnr_sim[sim, age] + vaxed_deaths_due_to_gnr_sim[sim, age]
      unvax_other_cause_deaths_sim[sim, age] <- rbinom(1, remaining_cohort, mu_ac[age]) - unvax_deaths_due_to_gnr_sim[sim, age]
      vaxed_other_cause_deaths_sim[sim, age] <- 0
      total_other_cause_deaths_sim[sim, age] <- unvax_other_cause_deaths_sim[sim, age] + vaxed_other_cause_deaths_sim[sim, age]
      
      # Calculates the number of survivors
      survived_sim[sim, age] <- remaining_cohort - total_deaths_due_to_gnr_sim[sim, age] - total_other_cause_deaths_sim[sim, age]
      survived_vaxed_sim[sim, age] <- survived_sim[sim, age] * v1_sample[sim]
      survived_unvax_sim[sim, age] <- survived_sim[sim, age] * (1 - v1_sample[sim])
      total_deaths_sim[sim, age] <- total_deaths_due_to_gnr_sim[sim, age] + total_other_cause_deaths_sim[sim, age]
      
      # Update the remaining cohort for the next period
      remaining_cohort <- survived_sim[sim, age]
      
      # Period 2 and beyond: some proportion vaccinated and protection conferred for 5 years
      vaxed_population <- survived_vaxed_sim[sim, age]
      unvax_population <- survived_unvax_sim[sim, age]
    } else {

      h_vec_sim[sim, age] <- rbeta(1, h_vec[age], h_n[age]-h_vec[age]) #For each age group, estimates a hospitalization probability based on the mean and standard deviation and 1 trial
      hospitalizations_unvax_sim[sim, age] <- rbinom(1, unvax_population, h_vec_sim[sim,age]) #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
      hospitalizations_vaxed_sim[sim, age] <- rbinom(1, vaxed_population, h_vec_sim[sim,age]*(1 - e1_sample[sim])) #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
      hospitalizations_total_sim[sim, age] <- hospitalizations_unvax_sim[sim, age] + hospitalizations_vaxed_sim[sim, age]
      
      mu_gnr_sim[sim, age] <- rbeta(1, mu_gnr_vec[age], mu_gnr_n[age] - mu_gnr_vec[age])
      unvax_deaths_due_to_gnr_sim[sim, age] <- rbinom(1, hospitalizations_unvax_sim[sim, age], mu_gnr_sim[sim,age])
      vaxed_deaths_due_to_gnr_sim[sim, age] <- rbinom(1, hospitalizations_vaxed_sim[sim, age], mu_gnr_sim[sim,age])
      total_deaths_due_to_gnr_sim[sim, age] <- unvax_deaths_due_to_gnr_sim[sim, age] + vaxed_deaths_due_to_gnr_sim[sim, age]
      unvax_other_cause_deaths_sim[sim, age] <- rbinom(1, unvax_population, mu_ac[age]) - unvax_deaths_due_to_gnr_sim[sim, age]
      vaxed_other_cause_deaths_sim[sim, age] <- rbinom(1, vaxed_population, mu_ac[age]) - vaxed_deaths_due_to_gnr_sim[sim, age]
      total_other_cause_deaths_sim[sim, age] <- unvax_other_cause_deaths_sim[sim, age] + vaxed_other_cause_deaths_sim[sim, age]
      
      # Calculates the number of survivors
      survived_unvax_sim[sim, age] <- unvax_population - unvax_deaths_due_to_gnr_sim[sim, age] - unvax_other_cause_deaths_sim[sim, age]
      survived_vaxed_sim[sim, age] <- vaxed_population - vaxed_deaths_due_to_gnr_sim[sim, age] - vaxed_other_cause_deaths_sim[sim, age]
      survived_sim[sim, age] <- survived_unvax_sim[sim, age] + survived_vaxed_sim[sim, age]
      total_deaths_sim[sim, age] <- total_deaths_due_to_gnr_sim[sim, age] + total_other_cause_deaths_sim[sim, age]
      
      # Update the remaining cohort for the next period
      unvax_population <- survived_unvax_sim[sim, age]
      vaxed_population <- survived_vaxed_sim[sim, age]
    }
  }
}

# Summarizes the results with confidence intervals
# First adds the total hospitalizations and deaths for each simulation
total_hospitalizations_sim <- rowSums(hospitalizations_sim)
total_deaths_due_to_gnr_sim <- rowSums(deaths_due_to_gnr_sim)
total_other_cause_deaths_sim <- rowSums(other_cause_deaths_sim)
total_all_cause_deaths_sim <- rowSums(total_deaths_sim)

# Establishes the confidence interval function
ci <- function(x) quantile(x, probs=c(0.025, 0.975))

total_hospitalizations_ci <- ci(total_hospitalizations_sim)
total_deaths_due_to_gnr_ci <- ci(total_deaths_due_to_gnr_sim)
total_other_cause_deaths_ci <- ci(total_other_cause_deaths_sim)
total_all_cause_deaths_ci <- ci(total_all_cause_deaths_sim)

# Prints the results
cat("Total hospitalizations due to GNR infections over 5 years:", mean(total_hospitalizations_sim), "\n")
cat("95% CI for total hospitalizations:", total_hospitalizations_ci, "\n\n")

cat("Total deaths due to GNR infections over 5 years:", mean(total_deaths_due_to_gnr_sim), "\n")
cat("95% CI for total deaths due to GNR:", total_deaths_due_to_gnr_ci, "\n\n")

cat("Total deaths due to other causes over 5 years:", mean(total_other_cause_deaths_sim), "\n")
cat("95% CI for total deaths due to other causes:", total_other_cause_deaths_ci, "\n\n")

cat("Total all-cause deaths over 5 years:", mean(total_all_cause_deaths_sim), "\n")
cat("95% CI for total all-cause deaths:", total_all_cause_deaths_ci, "\n\n")

cat("Number of individuals surviving at the end of 5 years:", mean(survived_sim[, periods]), "\n")
cat("95% CI for survivors:", ci(survived_sim[, periods]), "\n\n")

# Display results by age group
cat("\nResults by period (mean values):\n")
cat("Period\tHospitalizations\tGNR Deaths\tOther-cause Deaths\tTotal Deaths\tSurvivors\n")
for (age in 1:periods) {
  cat(age, "\t", mean(hospitalizations_sim[, age]), "\t\t", mean(deaths_due_to_gnr_sim[, age]), "\t\t", mean(other_cause_deaths_sim[, age]), "\t\t", mean(total_deaths_sim[, age]), "\t\t", mean(survived_sim[, age]), "\n")
}
