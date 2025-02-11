#Economic considerations

#ICER, incremental cost effectiveness ratio
#ratio between incremental DALYs averted by a scenario, and the incremental cost incurred

#DALYs averted
daly_diff_epi <- daly_status_quo - daly_epi
daly_diff_mat <- daly_status_quo - daly_mat
daly_diff_joint <- daly_status_quo - daly_joint

#Establish a function that outputs the DALYs averted, like the coverage plot
#Establish a function that outputs the incremental cost, given vaccine costs
mat_dose_cost <- 1
mat_doses_admin <- birth_cohort*vm 
cost_mat <- mat_doses_admin*mat_dose_cost

epi_dose_cost <- 1
epi_doses_admin <- birth_cohort*(1-mu_ac[1])*v1*3 #coverage applied to those who survive period 1 
cost_epi <- epi_doses_admin*epi_dose_cost

cost_joint <- cost_mat + cost_epi

icer_epi <- daly_diff_epi / cost_epi
icer_mat <- daly_diff_mat / cost_mat


age_at_death <- c(2/12, 8/12, 18/12, 42/12)
life_expectancy <- c(63.2, 63.2, 63.2, 63.2)
remaining_years <- life_expectancy - age_at_death
discounted_remaining_years <- (1 - exp(-discount_rate * remaining_years)) / discount_rate

discounted_YLL[sim, age] <- (1 - exp(-discount_rate*remaining_years[age]))

discount_rate <- .03



calculate_DALY_discounted <- function(deaths, age_at_death, life_expectancy, discount_rate) {
remaining_years <- life_expectancy - age_at_death
discounted_YLL <- (1 - exp(-discount_rate * remaining_years)) / discount_rate
return (deaths * discounted_YLL)
}

total_deaths_due_to_gnr <- c(8, 4, 2, 2)

calculate_DALY_undiscounted <- function(total_deaths_due_to_gnr, age_at_death, life_expectancy, age, periods) {
  for (age in 1:periods) {
  remaining_years[age] <- life_expectancy - age_at_death[age]
  daly <-(total_deaths_due_to_gnr[age] * remaining_years[age])
  return(daly[1] + daly[2])
 }
}

calculate_DALY_undiscounted(total_deaths_due_to_gnr, age_at_death, life_expectancy, age, periods)

for (age in 1:periods) {
  remaining_years <- life_expectancy - age_at_death[age]
  daly <-(total_deaths_due_to_gnr[age] * remaining_years)}


calculate_DALY(deaths=100, age_at_death[2], life_expectancy, discount_rate)

for 1 in n {
  ()
}


# Function to calculate Years of Life Lost (YLL) with 3% discounting
calculate_YLL <- function(deaths, age_at_death, life_expectancy = 63.2, discount_rate = 0.03) {
  remaining_years <- life_expectancy - age_at_death
  if (remaining_years <= 0) return(0)  # No YLL if past life expectancy
  
  # Apply discounting formula
  discounted_YLL <- (1 - exp(-discount_rate * remaining_years)) / discount_rate
  return(deaths * discounted_YLL)
}

# Function to calculate Years Lived with Disability (YLD) with 3% discounting
calculate_YLD <- function(incidence, disability_weight, duration, discount_rate = 0.03) {
  # Apply discounting formula
  discounted_YLD <- (1 - exp(-discount_rate * duration)) / discount_rate
  return(incidence * disability_weight * discounted_YLD)
}

# Function to calculate total DALYs
calculate_DALY <- function(deaths, age_at_death, incidence, disability_weight, duration, life_expectancy = 63.2, discount_rate = 0.03) {
  YLL <- calculate_YLL(deaths, age_at_death, life_expectancy, discount_rate)
  YLD <- calculate_YLD(incidence, disability_weight, duration, discount_rate)
  return(YLL + YLD)
}

# Example Usage
deaths <- 100      # Number of deaths
age_at_death <- 40 # Average age of death
incidence <- 500   # Number of incident cases
disability_weight <- 0.2  # Disability weight (ranges from 0 to 1)
duration <- 10     # Duration of disability in years

# Calculate DALYs
total_DALYs <- calculate_DALY(deaths, age_at_death, incidence, disability_weight, duration)
print(paste("Total DALYs:", round(total_DALYs, 2)))













