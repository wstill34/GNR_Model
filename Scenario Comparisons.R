#Calculations to compare vaccine scenarios

# EPI vs Status Quo
#Hospitalizations
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


# Sample v1 and e1
v1_sample_epi[sim] <- rnorm(1, mean=v1, sd=v1_se)
e1_sample_epi[sim] <- rnorm(1, mean=e1, sd=e1_se)
#REMOVE VARIABILITY AROUND E1, EFFICACY FOR BOTH EPI AND MATERNAL VACCINES

# Ensure probabilities remain within [0, 1]
v1_sample_epi[sim] <- max(min(v1_sample_epi[sim], 1), 0)
e1_sample_epi[sim] <- max(min(e1_sample_epi[sim], 1), 0)