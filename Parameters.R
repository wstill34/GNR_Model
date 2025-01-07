#Establishes parameters for all scenarios

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

# Defines number of simulations for Monte Carlo model
n_simulations <- 1000

# Sets childhood vaccination parameters
v1 <- 0.8   # Proportion vaccinated at end of period 1
v1_se <- 0.08 # Standard error for v1, for now 10% of v1
e1 <- 0.7   # Vaccine efficacy in reducing hospitalizations
e1_se <- 0.07  # Standard error for e1, for now 10% of e1

# Sets maternal vaccination parameters
vm <- 0.43  # Proportion of mothers vaccinated
vm_se <- 0.043 # Standard error for v_m
em <- 0.75  # Efficacy of maternal vaccine
em_se <- 0.075 # Standard error for e_m

# Sets range of coverage options for coverage vs reduced hospitalizations graphs
min_coverage <- 0.0
step <- 0.1
max_coverage <- 1.0


