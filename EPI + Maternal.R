# Sets the initial parameters birth cohort (number of individuals at time 0) and 
# periods of time (4 periods: 0-<4 months, 4-<12 months, 12-<24 months, 24-<60 months)
birth_cohort <- 112680
periods <- 4

# Sets age-specific parameters for each period
h <- c(0.0002003, 0.0000448, 0.0000279, 0.0000491)  # Probability of hospitalization due to GNR infections
h_se <- c(0.00004476, 0.00002117, 0.00001669, 0.00002216)  # Standard error for h
mu1 <- c(0.0348, 0.0348, 0.004725, 0.02835) # All-cause mortality probability
mu1_se <- c(0.00348, 0.00348, 0.0004725, 0.002835) # Standard error for mu1, for now just 10% of all-cause mortality probability, discuss with Meagan
mu2 <- c(0.5658, 0.3529, 0.5714, 0.3529)    # Probability of death given hospitalization due to GNR infection
mu2_se <- c(0.05658, 0.03529, 0.05714, 0.03529)    # Standard error for mu2, for now just 10% of GNR mortality probability, discuss with Meagan


# Sets vaccination parameters (childhood vaccine)
v1 <- 0.8   # Proportion vaccinated at end of year 1
v1_se <- 0.08 # Standard error for v1
e1 <- 0.7   # Vaccine efficacy in reducing hospitalizations
e1_se <- 0.07  # Standard error for e1

# Sets maternal vaccination parameters
v_m <- 0.43  # Proportion of mothers vaccinated
v_m_se <- 0.043 # Standard error for v_m
e_m <- 0.75  # Efficacy of maternal vaccine
e_m_se <- 0.075 # Standard error for e_m

# Defines number of simulations for Monte Carlo model
n_simulations <- 1000

# Creates matrices to store the results from simulations
hospitalizations_sim <- matrix(0, nrow=n_simulations, ncol=periods)
deaths_due_to_gnr_sim <- matrix(0, nrow=n_simulations, ncol=periods)
other_cause_deaths_sim <- matrix(0, nrow=n_simulations, ncol=periods)
total_deaths_sim <- matrix(0, nrow=n_simulations, ncol=periods)
survived_sim <- matrix(0, nrow=n_simulations, ncol=periods)

# Monte Carlo simulation
for (sim in 1:n_simulations) {
  remaining_cohort <- birth_cohort
  
  # Sample v1, e1, v_m, and e_m
  v1_sample <- rnorm(1, mean=v1, sd=v1_se)
  e1_sample <- rnorm(1, mean=e1, sd=e1_se)
  v_m_sample <- rnorm(1, mean=v_m, sd=v_m_se)
  e_m_sample <- rnorm(1, mean=e_m, sd=e_m_se)
  
  # Ensure probabilities remain within [0, 1]
  v1_sample <- max(min(v1_sample, 1), 0)
  e1_sample <- max(min(e1_sample, 1), 0)
  v_m_sample <- max(min(v_m_sample, 1), 0)
  e_m_sample <- max(min(e_m_sample, 1), 0)
  
  for (age in 1:periods) {
    # Sample h, mu1, and mu2 with uncertainty
    h_sample <- rnorm(1, mean=h[age], sd=h_se[age])
    mu1_sample <- rnorm(1, mean=mu1[age], sd=mu1_se[age])
    mu2_sample <- rnorm(1, mean=mu2[age], sd=mu2_se[age])
    
    # Ensure probabilities remain within [0, 1]
    h_sample <- max(min(h_sample, 1), 0)
    mu1_sample <- max(min(mu1_sample, 1), 0)
    mu2_sample <- max(min(mu2_sample, 1), 0)
    
    if (age == 1) {
      # Period 1: applies maternal vaccine protection to this group only
      vaccinated_mothers <- birth_cohort * v_m_sample
      unvaccinated_mothers <- birth_cohort * (1 - v_m_sample)
      
      # Babies born to vaccinated mothers have reduced hospitalization risk
      hospitalizations_vaccinated <- vaccinated_mothers * h_sample * (1 - e_m_sample)
      hospitalizations_unvaccinated <- unvaccinated_mothers * h_sample
      
      # Total hospitalizations in Period 1
      hospitalizations_sim[sim, age] <- hospitalizations_vaccinated + hospitalizations_unvaccinated
    } else {
      # Period 2 and beyond: applies childhood vaccine protection to these groups only
      vaccinated_population <- remaining_cohort * v1_sample
      unvaccinated_population <- remaining_cohort * (1 - v1_sample)
      
      # Vaccinated group has reduced probability of hospitalization
      hospitalizations_vaccinated <- vaccinated_population * h_sample * (1 - e1_sample)
      hospitalizations_unvaccinated <- unvaccinated_population * h_sample
      
      # Total hospitalizations
      hospitalizations_sim[sim, age] <- hospitalizations_vaccinated + hospitalizations_unvaccinated
    }
    
    # Calculates deaths due to GNR and other causes
    deaths_due_to_gnr_sim[sim, age] <- hospitalizations_sim[sim, age] * mu2_sample
    other_cause_deaths_sim[sim, age] <- remaining_cohort * mu1_sample - deaths_due_to_gnr_sim[sim, age]
    total_deaths_sim[sim, age] <- deaths_due_to_gnr_sim[sim, age] + other_cause_deaths_sim[sim, age]
    
    # Calculates the number of survivors
    survived_sim[sim, age] <- remaining_cohort - total_deaths_sim[sim, age]
    
    # Updates the remaining cohort for the next period
    remaining_cohort <- survived_sim[sim, age]
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

# Displays results by year
cat("\nResults by period (mean values):\n")
cat("Period\tHospitalizations\tGNR Deaths\tOther-cause Deaths\tTotal Deaths\tSurvivors\n")
for (age in 1:periods) {
  cat(age, "\t", mean(hospitalizations_sim[, age]), "\t\t", mean(deaths_due_to_gnr_sim[, age]), "\t\t", mean(other_cause_deaths_sim[, age]), "\t\t", mean(total_deaths_sim[, age]), "\t\t", mean(survived_sim[, age]), "\n")
}
