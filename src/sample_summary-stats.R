# Function to simulate summary statistics (accuracy, mean RT, variance RT)
# based on the EZ-diffusion model parameters.
simSumStats <- function(ddm_parameters, n_trials){
  # Extract DDM parameters: drift rate (v), boundary separation (a), non-decision time (t)
  v <- ddm_parameters$drift
  a <- ddm_parameters$bound
  t <- ddm_parameters$nondt
  # Number of trials per condition
  N <- n_trials
  # Number of parameter sets (conditions)
  nP <- length(v)
  
  # Intermediate calculation used in EZ equations
  y <- exp(-a * v)
  
  # Calculate predicted summary statistics using the forward EZ-diffusion equations
  # Predicted probability of correct response (Accuracy Rate)
  PredAccuracyRate <- 1 / (1 + y)  # Equation 1
  # Predicted mean reaction time (MRT) for correct responses
  PredMean <- t + ((a / (2 * v)) * ((1 - y) / (1 + y)))  # Equation 2
  # Predicted variance of reaction time (VRT) for correct responses
  PredVariance <- (a / (2 * v^3)) * (( 1 - 2 * a * v * y - exp(-a*2*v)) / ((y + 1)^2))  # Equation 3
  
  # Simulate observed summary statistics by sampling from distributions
  # based on the predicted values and number of trials (N).
  
  # Simulate the total number of correct responses (Accuracy) from a binomial distribution.
  # The probability of success is the predicted accuracy rate.
  ObservedAccuracyTotal <- rbinom(nP, size = N, prob = PredAccuracyRate)
  
  # Simulate the observed mean reaction time (MRT) from a normal distribution.
  # The mean is the predicted MRT, and the standard deviation accounts for sampling variability (sqrt(PredVariance / N)).
  ObservedMean <- rnorm(nP, mean = PredMean, sd = sqrt(PredVariance / N))
  
  # Simulate the observed variance of reaction time (VRT) from a normal distribution.
  # The mean is the predicted VRT. The standard deviation is derived from the standard error of the variance estimate.
  # Note: This assumes the underlying RT distribution is approximately normal, which might be a simplification.
  ObservedVariance <- rnorm(nP, mean = PredVariance, sd = sqrt((2 * (PredVariance^2)) / (N - 1)))
  
  # Return a data frame containing the simulated observed summary statistics
  return(data.frame("A" = ObservedAccuracyTotal, # Simulated total correct responses
                    "Mrt" = ObservedMean,       # Simulated mean reaction time
                    "Vrt" = ObservedVariance))  # Simulated variance of reaction time
}