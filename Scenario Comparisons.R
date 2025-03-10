#Calculations to compare vaccine scenarios

# Maternal vs Status Quo
#Hospitalizations
hosp_diff_mat_sq <- sq_hosps_deaths_dalys[,1] - maternal_hosps_deaths_dalys[,1]
mean(hosp_diff_mat_sq) #Outputs Mean hosp averted
ci(hosp_diff_mat_sq) #Outputs CI hosp averted

death_diff_mat_sq <- sq_hosps_deaths_dalys[,2] - maternal_hosps_deaths_dalys[,2]
mean(death_diff_mat_sq) #Outputs Mean death averted
ci(death_diff_mat_sq) #Outputs CI death averted

daly_diff_mat_sq <- sq_hosps_deaths_dalys[,3] - maternal_hosps_deaths_dalys[,3]
mean(daly_diff_mat_sq) #Outputs Mean daly averted
ci(daly_diff_mat_sq) #Outputs CI daly averted

# EPI vs Status Quo
#Hospitalizations
hosp_diff_epi_sq <- sq_hosps_deaths_dalys[,1] - epi_hosps_deaths_dalys[,1]
mean(hosp_diff_epi_sq) #Outputs Mean hosp averted
ci(hosp_diff_epi_sq) #Outputs CI hosp averted

death_diff_epi_sq <- sq_hosps_deaths_dalys[,2] - epi_hosps_deaths_dalys[,2]
mean(death_diff_epi_sq) #Outputs Mean death averted
ci(death_diff_epi_sq) #Outputs CI death averted

daly_diff_epi_sq <- sq_hosps_deaths_dalys[,3] - epi_hosps_deaths_dalys[,3]
mean(daly_diff_epi_sq) #Outputs Mean daly averted
ci(daly_diff_epi_sq) #Outputs CI daly averted

# Joint vs Status Quo
#Hospitalizations
hosp_diff_joint_sq <- sq_hosps_deaths_dalys[,1] - joint_hosps_deaths_dalys[,1]
mean(hosp_diff_joint_sq) #Outputs Mean hosp averted
ci(hosp_diff_joint_sq) #Outputs CI hosp averted

death_diff_joint_sq <- sq_hosps_deaths_dalys[,2] - joint_hosps_deaths_dalys[,2]
mean(death_diff_joint_sq) #Outputs Mean death averted
ci(death_diff_joint_sq) #Outputs CI death averted

daly_diff_joint_sq <- sq_hosps_deaths_dalys[,3] - joint_hosps_deaths_dalys[,3]
mean(daly_diff_joint_sq) #Outputs Mean daly averted
ci(daly_diff_joint_sq) #Outputs CI daly averted

#All costs
#SQ
mean(sq_hosps_deaths_dalys[,1])*157.50 #Outputs expected medical costs, sq
ci(sq_hosps_deaths_dalys[,1])*157.50 #Outputs confidence interval of medical costs, sq

#Maternal
mean(maternal_hosps_deaths_dalys[,1])*157.50 #Outputs expected medical costs, maternal
ci(maternal_hosps_deaths_dalys[,1])*157.50 #Outputs CI of medical costs, maternal
coverage_mat*birth_cohort*1.47 #Outputs maternal vaccine cost
mean(maternal_hosps_deaths_dalys[,1])*157.50 + coverage_mat*birth_cohort*1.47 #Outputs expected total costs, maternal
ci(maternal_hosps_deaths_dalys[,1])*157.50 + c(coverage_mat*birth_cohort*1.47,coverage_mat*birth_cohort*1.47) #Outpus CI of total cost, maternal

#EPI
mean(epi_hosps_deaths_dalys[,1])*157.50 #Outputs expected medical costs, epi
ci(epi_hosps_deaths_dalys[,1])*157.50 #Outputs CI of medical costs, epi
birth_cohort*(1-mu_ac[1])*coverage_epi*2*1.47 #Outputs epi vaccine cost
mean(epi_hosps_deaths_dalys[,1])*157.50 + birth_cohort*(1-mu_ac[1])*coverage_epi*2*1.47 #Outputs expected total cost, epi
ci(epi_hosps_deaths_dalys[,1])*157.50 + c(birth_cohort*(1-mu_ac[1])*coverage_epi*2*1.47,birth_cohort*(1-mu_ac[1])*coverage_epi*2*1.47) #Outputs CI of total cost, epi

#Joint
mean(joint_hosps_deaths_dalys[,1])*157.50 #Outputs expected medical costs, joint
ci(joint_hosps_deaths_dalys[,1])*157.50 #Outputs CI of medical costs, joint
coverage_mat*birth_cohort*1.47 + birth_cohort*(1-mu_ac[1])*coverage_epi*2*1.47 #Outputs joint vaccine cost
mean(joint_hosps_deaths_dalys[,1])*157.50 + (coverage_mat*birth_cohort*1.47 + birth_cohort*(1-mu_ac[1])*coverage_epi*2*1.47) #Outputs expected total cost, epi
ci(joint_hosps_deaths_dalys[,1])*157.50 + c((coverage_mat*birth_cohort*1.47 + birth_cohort*(1-mu_ac[1])*coverage_epi*2*1.47),(coverage_mat*birth_cohort*1.47 + birth_cohort*(1-mu_ac[1])*coverage_epi*2*1.47)) #Outputs CI of total cost, epi

#Difference between joint and maternal costs
(mean(joint_hosps_deaths_dalys[,1])*157.50 + coverage_mat*birth_cohort*1.47 + birth_cohort*(1-mu_ac[1])*coverage_epi*2*1.47) - (mean(maternal_hosps_deaths_dalys[,1])*157.50 + coverage_mat*birth_cohort*1.47)
#Difference between DALYs averted, joint and maternal
mean(maternal_hosps_deaths_dalys[,3]) - mean(joint_hosps_deaths_dalys[,3])
#ICER estimate
((mean(joint_hosps_deaths_dalys[,1])*157.50 + coverage_mat*birth_cohort*1.47 + birth_cohort*(1-mu_ac[1])*coverage_epi*2*1.47) - (mean(maternal_hosps_deaths_dalys[,1])*157.50 + coverage_mat*birth_cohort*1.47))/(mean(maternal_hosps_deaths_dalys[,3]) - mean(joint_hosps_deaths_dalys[,3]))





reduced_hosp_epi <- status_quo_data[1,1] - epi_data[1,1]
std_dev_hosp_epi_sq <- sqrt(epi_data[1,2]^2 + status_quo_data[1,2]^2)

#Deaths
reduced_deaths_epi <- status_quo_data[1,3] - epi_data[1,3]
std_dev_deaths_epi_sq <- sqrt(epi_data[1,4]^2 + status_quo_data[1,4]^2)

# Maternal vs Status Quo
#Hospitalizations
reduced_hosp_mat <- status_quo_data[1,1] - maternal_data[1,1]
std_dev_hosp_mat_sq <- sqrt(maternal_data[1,2]^2 + status_quo_data[1,2]^2)

#Deaths
reduced_deaths_mat <- status_quo_data[1,3] - maternal_data[1,3]
std_dev_deaths_mat_sq <- sqrt(maternal_data[1,4]^2 + status_quo_data[1,4]^2)

# Maternal vs Status Quo
#Hospitalizations
reduced_hosp_joint <- status_quo_data[1,1] - joint_data[1,1]
std_dev_hosp_joint_sq <- sqrt(joint_data[1,2]^2 + status_quo_data[1,2]^2)

#Deaths
reduced_deaths_joint <- status_quo_data[1,3] - joint_data[1,3]
std_dev_deaths_joint_sq <- sqrt(joint_data[1,4]^2 + status_quo_data[1,4]^2)

#DALYs
#Calculates disability-adjusted life-years, sum of difference between life expectancy and duration of life for each death due to gnr
#Life expectancy: 63.2 years https://www.cia.gov/the-world-factbook/countries/mali/#:~:text=Life%20expectancy%20at%20birth,63.2%20years%20(2024%20est.)
#Age 1: median 2 months; Age 2: median 8 months; Age 3: median 18 months; Age 4: median 42 months
daly_status_quo <- status_quo_data[1, 9]*(63.2 - 2/12) + status_quo_data[1, 10]*(63.2 - 8/12) + status_quo_data[1, 11]*(63.2 - 18/12) + status_quo_data[1, 12]*(63.2 - 42/12)
daly_epi <- epi_data[1, 9]*(63.2 - 2/12) + epi_data[1, 10]*(63.2 - 8/12) + epi_data[1, 11]*(63.2 - 18/12) + epi_data[1, 12]*(63.2 - 42/12)
daly_mat <- maternal_data[1, 9]*(63.2 - 2/12) + maternal_data[1, 10]*(63.2 - 8/12) + maternal_data[1, 11]*(63.2 - 18/12) + maternal_data[1, 12]*(63.2 - 42/12)
daly_joint <- joint_data[1, 9]*(63.2 - 2/12) + joint_data[1, 10]*(63.2 - 8/12) + joint_data[1, 11]*(63.2 - 18/12) + joint_data[1, 12]*(63.2 - 42/12)

daly_diff_maternal <- status_quo_data - maternal_data
hist(daly_diff_maternal)
ci(daly_diff_maternal)

#Difference between maternal and status quo
mat_diff <- status_quo_hosps_deaths_dalys - maternal_hosps_deaths_dalys
mean(mat_diff[,1])
ci(mat_diff[,1])
mean(mat_diff[,2])
ci(mat_diff[,2])
mean(mat_diff[,3])
ci(mat_diff[,3])

#Difference between pediatric and status quo
epi_diff <- status_quo_hosps_deaths_dalys - epi_hosps_deaths_dalys
mean(epi_diff[,1])
ci(epi_diff[,1])
mean(epi_diff[,2])
ci(epi_diff[,2])
mean(epi_diff[,3])
ci(epi_diff[,3])

#Difference between joint and status quo
joint_diff <- status_quo_hosps_deaths_dalys - joint_hosps_deaths_dalys
mean(joint_diff[,1])
ci(joint_diff[,1])
mean(joint_diff[,2])
ci(joint_diff[,2])
mean(joint_diff[,3])
ci(joint_diff[,3])

simulate_maternal(birth_cohort, periods, h, mu_ac, mu_gnr, n_simulations, coverage_mat, efficacy_mat, remaining_years, mat_cost_per_dose)[,1]