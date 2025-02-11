#Establishes parameters for all scenarios

birth_cohort <- 112680
periods <- 4

h_vec <- c(71, 20, 12, 21) # GNR hospitalized patients by age group
h_n <- c(126449, 252898, 502662, 1038045) # n for age group

mu_gnr_vec <- c(41, 7, 7, 8)
mu_gnr_n <- c(71, 19, 12, 19)

# Sets age-specific parameters for each period
h <- h_vec / h_n  # Probability of hospitalization due to GNR infections
mu_gnr <- mu_gnr_vec / mu_gnr_n    # Probability of death given hospitalization due to GNR infection
mu_ac <- c(0.0348, 0.0348, 0.004725, 0.02835) # All-cause mortality probability
length <- c("0 - <4 months ", "4 - <12 months", "12 - <24 months", "24 - 59 months")

# Defines number of simulations for Monte Carlo model
n_simulations <- 10000

# Sets childhood vaccination parameters
#80.1% of children born 7 to 12 months after introduction received three doses of pentavalent rotavirus vaccine by 9 months of age
v1 <- 0.801   # Proportion vaccinated at end of period 1
#v1_se <- 0.08 # Standard error for v1, for now 10% of v1
#e1 <- 0.7   # Vaccine efficacy in reducing hospitalizations
#e1_se <- 0.07  # Standard error for e1, for now 10% of e1

# Sets maternal vaccination parameters
#In Mali, according to the 2018 Demographic and Health Survey (DHS), only about 43% of pregnant women had at least four antenatal care visits, which is the internationally recommended standard
vm <- 0.433  # Proportion of mothers vaccinated
#vm_se <- 0.043 # Standard error for v_m
#em <- 0.75  # Efficacy of maternal vaccine
#em_se <- 0.075 # Standard error for e_m

# Sets range of coverage options for coverage vs reduced hospitalizations graphs
min_coverage <- 0.0
step <- 0.1
max_coverage <- 1.0

#Based on PCV vaccine pooled estimate from https://pmc.ncbi.nlm.nih.gov/articles/instance/7332418/bin/41586_2020_2238_MOESM1_ESM.pdf
#81.2% efficacy (95% CI: 63.1%, 90.5%)

#For uncertainty distribution for EPI efficacy:
a <- rnorm(n_simulations, mean = log(.188), sd = ((log(.369)-log(.095))/3.98))
efficacy_epi <- 1 - exp(a)
hist(efficacy_epi)
mean(efficacy_epi)

#Based on maternal influenza vaccine estimate from https://pubmed.ncbi.nlm.nih.gov/31259452/
#56.8% efficacy (95% CI: 25.0, 75.1)

#For uncertainty distribution for Maternal efficacy:
b <- rnorm(n_simulations, mean = log(.432), sd = ((log(.750)-log(.249))/3.98))
efficacy_mat <- 1 - exp(b)
hist(efficacy_mat)
mean(efficacy)

mat_cost_per_dose <- 1
discount_rate <- .03
age_at_death <- c(2/12, 8/12, 18/12, 42/12)
life_expectancy <- c(63.2, 63.2, 63.2, 63.2)
remaining_years <- life_expectancy - age_at_death
discounted_remaining_years <- (1 - exp(-discount_rate * remaining_years)) / discount_rate

h_vec_sim <- matrix(0, nrow=n_simulations, ncol=periods)
mu_gnr_sim <- matrix(0, nrow=n_simulations, ncol=periods)
#mu_ac_sim <- matrix(0, nrow=n_simulations, ncol=periods)

for (sim in 1:n_simulations) {
  for (age in 1:periods) {
    # Calculates number of hospitalizations, deaths due to GNR, and other cause deaths
    h_vec_sim[sim, age] <- rbeta(1, h_vec[age], h_n[age]-h_vec[age]) #For each age group, estimates a hospitalization probability based on the mean and standard deviation and 1 trial
    mu_gnr_sim[sim, age] <- rbeta(1, mu_gnr_vec[age], mu_gnr_n[age] - mu_gnr_vec[age])
    #mu_ac_sim[sim, age] <- rbeta(1, mu_ac_vec[age]
    }
}

ci <- function(x) quantile(x, probs=c(0.025, 0.975))

install.packages("writexl")

#library("writexl")
#write_xlsx(status_quo_data,"C:/Users/wstil/OneDrive/Desktop/Aim 3/GNR Model/status_quo_df.xlsx", col_name=TRUE, format_headers=FALSE)

