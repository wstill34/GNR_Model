# Sets the initial parameters birth cohort (number of individuals at time 0) and 
# periods of time (4 periods: 0-<4 months, 4-<12 months, 12-<24 months, 24-<60 months)
birth_cohort <- 112680
periods <- 4

# Sets age-specific parameters for each period
h <- c(71/126449, 20/252898, 12/502662, 21/1038045)  # Probability of hospitalization due to GNR infections
h_se <- c(0.0000432584, 0.0000229607, 0.0000154507, 0.0000246348)  # Standard error for h
mu_ac <- c(0.0348, 0.0348, 0.004725, 0.02835) # All-cause mortality probability
mu_gnr <- c(41/71, 7/19, 7/12, 8/19)    # Probability of death given hospitalization due to GNR infection
mu_gnr_se <- c(0.13874, 0.175219, 0.341565, 0.195646)    # Standard error for mu2, for now just 10% of GNR mortality probability, discuss with Meagan
length <- c("0 - <4 months ", "4 - <12 months", "12 - <24 months", "24 - 59 months")

# Defines number of simulations for Monte Carlo model
n_simulations <- 1000

# Creates matrices to store the results from simulations
hospitalizations_sim <- matrix(0, nrow=n_simulations, ncol=periods)
#Added h_vec_sim matrix of 1000 rows, 4 columns
h_vec_sim <- matrix(0, nrow=n_simulations, ncol=periods)
mu_gnr_sim <- matrix(0, nrow=n_simulations, ncol=periods)
deaths_due_to_gnr_sim <- matrix(0, nrow=n_simulations, ncol=periods)
other_cause_deaths_sim <- matrix(0, nrow=n_simulations, ncol=periods)
total_deaths_sim <- matrix(0, nrow=n_simulations, ncol=periods)
survived_sim <- matrix(0, nrow=n_simulations, ncol=periods)

# Monte Carlo simulation
for (sim in 1:n_simulations) {
  remaining_cohort <- birth_cohort
 
  for (age in 1:periods) {
   # Calculates number of hospitalizations, deaths due to GNR, and other cause deaths
    h_vec_sim[sim, age] <- rnorm(1, mean = h[age], sd = h_se[age]) #For each age group, estimates a hospitalization probability based on the mean and standard deviation and 1 trial
    #QUESTION FOR MEAGAN - what about negative h_vec_sim values? do we need to add a lower limit of 0? No right, because the binomial distribution starts at 0 so can't have negative hospitalizations
    hospitalizations_sim[sim, age] <- rbinom(1, remaining_cohort, h_vec_sim[age]) #for each age group, estimates number of hospitalized babies based on h_vec_sim and binomial distribution
    mu_gnr_sim[sim, age] <- rnorm(1, mean= mu_gnr[age], sd = mu_gnr_se[age])
    deaths_due_to_gnr_sim[sim, age] <- rbinom(1, hospitalizations_sim[sim, age], mu_gnr_sim[age])
    #QUESTION FOR MEAGAN - do we need to create a normal distribution for deaths as well, followed by binomial, or just binomial here sufficient?
    other_cause_deaths_sim[sim, age] <- rbinom(1, remaining_cohort, mu_ac[age]) - deaths_due_to_gnr_sim[sim, age]
     
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

# Establishes the confidence interval function
ci <- function(x) quantile(x, probs=c(0.025, 0.975))


total_hospitalizations_ci <- ci(total_hospitalizations_sim)
total_deaths_due_to_gnr_ci <- ci(total_deaths_due_to_gnr_sim)
total_other_cause_deaths_ci <- ci(total_other_cause_deaths_sim)
total_all_cause_deaths_ci <- ci(total_all_cause_deaths_sim)

# Print the results
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

# Displays results by age group
cat("\nResults by period (mean values):\n")
cat("Period\tDuration\tHospitalizations\tGNR Deaths\tOther-cause Deaths\tSurvivors\n")
for (age in 1:periods) {
  cat(age, "\t", length[age], "\t", round(mean(hospitalizations_sim[, age])), "\t", round(mean(deaths_due_to_gnr_sim[, age])), "\t", round(mean(other_cause_deaths_sim[, age])), "\t", round(mean(survived_sim[, age])), "\n")
}
#QUESTION FOR MEAGAN, results are more hospitalizations than I would expect, discuss reasons why/if this is a probable error?
