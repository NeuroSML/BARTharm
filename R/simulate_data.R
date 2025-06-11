
# This function generates synthetic cross-sectional data for a study with scanner-related batch effects,
# covariates (age, sex), biological features, and outcomes. It flexibly allows for correlated scanner assignment,
# linear or nonlinear relationships, and outcome contamination by image quality metrics (IQMs).

# Arguments:
# - n_subjects Number of simulated subjects (default = 1000)
# - linear_tau Logical. If TRUE, outcome depends linearly on covariates
# - linear_mu Logical. If TRUE, scanner effects on outcome are modeled linearly

simulate_data <- function(n_subjects = 1000,
                          linear_tau = TRUE,
                          linear_mu = TRUE) {
  
  # Step 1: Define scanner-level properties (resolution and noise_level)
  scanner_properties <- data.frame(scanner_id = 1:10)
  scanner_properties$scanner_resolution <- ifelse(scanner_properties$scanner_id <= 3, "low",
                                                  ifelse(scanner_properties$scanner_id <= 7, "medium", "high"))
  scanner_properties$noise_level <- ifelse(scanner_properties$scanner_id <= 3, 3,
                                           ifelse(scanner_properties$scanner_id <= 7, 2, 1))
  
  # Step 2: Simulate demographic variables
  simdata <- data.frame(
    subid = 1:n_subjects,
    age = sample(20:60, n_subjects, replace = TRUE),
    sex = sample(c(0, 1), n_subjects, replace = TRUE)
  )
  
  # Assign scanner_id based on correlation setting
  simdata$scanner_id <- sample(1:10, n_subjects, replace = TRUE)
  
  # Optionally add binary group variable
  simdata$group <- sample(c(0, 1), n_subjects, replace = TRUE)
  
  # Merge scanner resolution and noise_level into subject data
  simdata <- dplyr::left_join(simdata, scanner_properties, by = "scanner_id")
  simdata$scanner_res_numeric <- ifelse(simdata$scanner_resolution == "low", 1,
                                        ifelse(simdata$scanner_resolution == "medium", 2, 3))
  
  # Step 3: Simulate true biological features
  simdata$bio_feature1_true <- runif(n_subjects, 0, 25)
  simdata$bio_feature2_true <- runif(n_subjects, 0, 30)
  nonlinear_1 <- simdata$bio_feature1_true^2
  nonlinear_2 <- log(simdata$bio_feature2_true + 1)
  
  # Step 4: Generate baseline outcome (before IQMs)
  random_effects <- rnorm(n_subjects)
  outcome <- random_effects
  simdata$random_noise <- random_effects
  
  if (linear_tau) {
    # Linear outcome with additive effects
    sex_effect <- ifelse(simdata$sex == 0, -25, 25)
    outcome <- outcome + sex_effect
    outcome <- outcome - 4 * simdata$bio_feature1_true + 6 * simdata$bio_feature2_true +
                     0.1 * simdata$bio_feature1_true * simdata$sex +
                     0.1 * simdata$bio_feature1_true * simdata$bio_feature2_true

   outcome <- outcome + 0.01 * simdata$bio_feature2_true * simdata$age
  } else {
    # Nonlinear outcome with interactions
    outcome <- ifelse(simdata$sex == 0, -50, 50)
    outcome <- outcome - 0.1 * simdata$age^2 + 5 * sin(simdata$age)
    outcome <- outcome - 2 * nonlinear_1 + 3 * nonlinear_2 +
        0.5 * simdata$bio_feature1_true * (simdata$bio_feature2_true^2)/2
    outcome <- outcome + 0.1 * simdata$bio_feature1_true * simdata$sex +
        0.01 * simdata$bio_feature1_true * simdata$bio_feature2_true

    outcome <- outcome + 0.01 * simdata$bio_feature2_true * simdata$age
  }
  
  
  # Step 5: Add scanner effects (mu)
  outcome_simulated <- outcome

  if (linear_mu) {
    outcome_simulated <- outcome_simulated - 0.2 * simdata$scanner_res_numeric +
      0.15 * simdata$noise_level + 0.05 * simdata$noise_level * simdata$scanner_res_numeric
  } else {
    outcome_simulated <- outcome_simulated - 5 * simdata$scanner_res_numeric^2 +
      2 * exp(-simdata$noise_level) + 0.5 * simdata$scanner_res_numeric * simdata$noise_level^2
  }
  
  # Step 6: Add measurement noise to biological features
  simdata$bio_feature1_noise <- rnorm(n_subjects, 0, simdata$noise_level * 0.25)
  simdata$bio_feature2_noise <- rnorm(n_subjects, 0, simdata$noise_level * 0.1)
  
  simdata$bio_feature1_observed <- simdata$bio_feature1_true + simdata$bio_feature1_noise
  simdata$bio_feature2_observed <- simdata$bio_feature2_true + simdata$bio_feature2_noise

  # Step 7: Simulate SNR and CNR metrics
  simdata$SNR <- mapply(function(res, noise) {
    base <- ifelse(res == "high", runif(1, 50, 65),
                   ifelse(res == "medium", runif(1, 30, 49), runif(1, 10, 29)))
    base - (3 - noise) * 5
  }, simdata$scanner_resolution, simdata$noise_level)
  
  simdata$CNR <- mapply(function(res, noise) {
    base <- ifelse(res == "high", runif(1, 5, 14),
                   ifelse(res == "medium", runif(1, 15, 29), runif(1, 30, 40)))
    base + (3 - noise) * 2
  }, simdata$scanner_resolution, simdata$noise_level)
  
  # Step 9: Finalize outcome with IQM effects 
  outcome_simulated <- outcome_simulated + 1.75 * simdata$SNR - 1.55 * simdata$CNR - 0.5 * simdata$SNR * simdata$CNR
  
  # Save raw outcome values
  simdata$outcome_clean <- outcome  #this corresponds to the true, harmonized, outcome Y harm
  simdata$outcome_bio <- outcome - random_effects #this corresponds to noiseless biological contribution TAU
  simdata$outcome_iqm <- outcome_simulated - outcome #this corresponds to scanner contribution MU
  simdata$outcome_simulated <- outcome_simulated #this corresponds to observed outcome Y

  simdata$scanner_id <- as.factor(simdata$scanner_id)
  
  data_bio <-  simdata %>% dplyr::select(subid, outcome_simulated, age, sex, bio_feature1_observed, bio_feature2_observed, group )
  data_nonbio <-  simdata %>% dplyr::select(subid, scanner_id,  noise_level, scanner_res_numeric, SNR, CNR)
  
  return(list("simdata" = simdata, "data_bio" = data_bio, "data_iqm" = data_nonbio))
}

