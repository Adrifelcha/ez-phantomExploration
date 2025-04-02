library(here)
library(knitr)
library(kableExtra) 
library(tidyverse)
source(here("src", "sample_summary-stats.R"))
library(R2jags)

K <- 2 # Number of factor levels
D <- 8 # Number of days (repetitions)
N <- 120 # Number of trials per condition

bound <- matrix(c(3,   4, 
                  3.05, 4.05), nrow = K, ncol = K, byrow = TRUE)
# Only Factor B has a main effect
drift <- matrix(c(0.7, 0.75, 
                  1.5, 1.45), nrow = K, ncol = K, byrow = TRUE)


layout(matrix(c(1,2), nrow = 2, ncol = 1))
par(mai = c(0.5,1,0.5,0.1), oma = c(0,0,0,0))
plot(c(0.5,1), bound[1,], type = "b", pch = 16, col = "red", ylim = c(0, max(bound)*1.2),
     ann = FALSE,axes = FALSE, xlim = c(0.4,1.1))
axis(1, at = c(0.5,1), labels = c("Factor A1", "Factor A2"))
axis(2, at = seq(0, max(bound)*1.2,length.out = 5), labels = round(seq(0, max(bound)*1.2,length.out = 5),1), las = 2)
points(c(0.5,1), bound[2,], type = "b", pch = 16, col = "blue")
legend("bottomright", legend = c("Factor B1", "Factor B2"), col = c("red", "blue"), pch = 16, bty = "n")
mtext("Only Factor A has a main effect", side = 3, line = 0, font = 2)
mtext(expression(bold(paste("Boundary separation ", alpha))), side = 2, line = 2.5, font = 2, cex = 1.2)

plot(c(0.5,1), drift[1,], type = "b", pch = 16, col = "red", ylim = c(0, max(drift)*1.2),
     ann = FALSE,axes = FALSE, xlim = c(0.4,1.1))
axis(1, at = c(0.5,1), labels = c("Factor A1", "Factor A2"))
axis(2, at = seq(0, max(drift)*1.2,length.out = 5), labels = round(seq(0, max(drift)*1.2,length.out = 5),1), las = 2)
points(c(0.5,1), drift[2,], type = "b", pch = 16, col = "blue")
mtext("Only Factor B has a main effect", side = 3, line = 0, font = 2)
mtext(expression(bold(paste("Drift rate ", nu))), side = 2, line = 2.5, font = 2, cex = 1.2)


ddm_parameters <- list(
  bound = rep(as.vector(bound), D),
  drift = rep(as.vector(drift), D),
  nondt = rep(0.2, D*K*K)  
)

FA <- rep(rep(paste("A",c(1,2), sep=""), each = K), D)
FB <- rep(paste("B",c(1,2), sep=""), K*D)

design <- data.frame("Drift" = ddm_parameters$drift, 
                     "Bound" = ddm_parameters$bound, 
                     "Factor A" = FA, "Factor B" = FB)

design$Condition <- 1
design$Condition[design$Factor.B == "B2" & design$Factor.A == "A1"] <- 2
design$Condition[design$Factor.B == "B1" & design$Factor.A == "A2"] <- 3
design$Condition[design$Factor.B == "B2" & design$Factor.A == "A2"] <- 4
head(design,8)

sum_stats <- simSumStats(ddm_parameters, N)
full_data <- cbind(design, sum_stats)

head(full_data)


modelFile <- here("output", "BUGS", "ezddm_model.bug")

model <- write("
    model {
              ####### Priors
              drift_mu ~ dnorm(0,1)                     # Baseline drift rate
              drift_lambda ~ dgamma(2,1)
              drift_sigma = pow(drift_lambda, -0.5)
              bound_mu ~ dnorm(2,1)T(0,)                     # Baseline boundary separation
              bound_lambda ~ dgamma(2,1)
              bound_sigma = pow(bound_lambda, -0.5)
              nondt ~ dexp(1)

              # Regression coefficients
              for(i in 1:3){                  
                  gamma[i] ~ dnorm(0,1)                # Drift rate
                  lambda[i] ~ dnorm(0,1)               # Boundary separation
              }      
              
              for(j in 1:4){
                  drift_pred[j] = drift_mu + gamma[1]*A[j]+gamma[2]*B[j]+gamma[3]*A[j]*B[j]
                  bound_pred[j] = bound_mu + lambda[1]*A[j]+lambda[2]*B[j]+lambda[3]*A[j]*B[j]
              }
              
              ####### Sampling model
              for (k in 1:length(meanRT)) {
                  # Person-by-condition parameters for DM parameters                                    
                  drift[k] ~ dnorm(drift_pred[cond[k]],drift_lambda)
                  bound[k] ~ dnorm(bound_pred[cond[k]],bound_lambda)                  
          
                  # Forward equations from EZ Diffusion
                  ey[k]  = exp(-bound[k] * drift[k])
                  Pc[k]  = 1 / (1 + ey[k])
                  PRT[k] = 2 * pow(drift[k], 3) / bound[k] * pow(ey[k] + 1, 2) / (2 * -bound[k] * 
                           drift[k] * ey[k] - ey[k] * ey[k] + 1)
                  MDT[k] = (bound[k] / (2 * drift[k])) * (1 - ey[k]) / (1 + ey[k])
                  MRT[k] = MDT[k] + nondt
          
                  # Sampling distributions for summary statistics
                  correct[k] ~ dbin(Pc[k], nTrials[k])
                  varRT[k]   ~ dnorm(1/PRT[k], 0.5*(nTrials[k]-1) * PRT[k] * PRT[k])
                  meanRT[k]  ~ dnorm(MRT[k], PRT[k] * nTrials[k])
                }
    }", modelFile)

# General setup
n.chains  <- 4;      n.iter    <- 5000
n.burnin  <- 250;    n.thin    <- 1

# Pass data to JAGS
data_toJAGS <- list("nTrials"  =  rep(N, length(full_data$Mrt)),
                    "meanRT"   =  full_data$Mrt,
                    "varRT"    =  full_data$Vrt,
                    "correct"  =  full_data$A,
                    "cond"     =  full_data$Condition,
                    "A"   =  as.integer(full_data$Factor.A)-1, 
                    "B"   =  as.integer(full_data$Factor.B)-1)

# Specify parameters to keep track of
parameters <- c('gamma', 'lambda', 'drift', 'bound', 'drift_pred', 'bound_pred',
                "Pc", "PRT", "MRT")

# Prepare initial values
myinits <- rep(list(list()), n.chains)
for(i in 1:n.chains){
    myinits[[i]] <- list(drift = rnorm(length(data_toJAGS$nTrials),1,0.1))
}


start <- Sys.time()
samples <- jags(data=data_toJAGS,
                parameters.to.save=parameters,
                model=modelFile,
                n.chains=n.chains,  n.iter=n.iter,
                n.burnin=n.burnin,  n.thin=n.thin,
                DIC=T, inits=myinits)
end <- Sys.time()

cat("Time taken:", round(difftime(end, start, units = "secs"), 2), "seconds")



source(here("src", "get_Rhat.R"))

rhats <- apply(samples$BUGSoutput$sims.array,3,getRhat)
rule <- 1.05
bad.Rhat <- which(rhats>rule)
test.rhat <- length(bad.Rhat) > 0
  if(test.rhat){
          par(mfrow=c(1,1))
          which.are.bad.Rhats <- names(bad.Rhat)
          hist(rhats, breaks = 50)
          abline(v=rule, col="red", lty=2)
          legend("top",paste("Rhat > ",rule," | ",
                             (round(nrow(bad.Rhat)/(length(as.vector(rhats))),5))*100,
                             "% of chains | ", length(which.are.bad.Rhats), " chains", sep=""), lty=2, col="red", cex=0.4)
          table(which.are.bad.Rhats)
  }else{      paste("No Rhat greater than ", rule, sep="")       }



  # Effects of instruction
gamma <- as.vector(samples$BUGSoutput$sims.list$gamma) # Main
lambda <- samples$BUGSoutput$sims.list$lambda # Interaction
# Fitted values / Predicted drift rates
drift_pred <- samples$BUGSoutput$sims.list$drift_pred
bound_pred <- samples$BUGSoutput$sims.list$bound_pred

# Load custom plotting functions
source(here("src", "plot_verticalHist.R"))
source(here("src", "plot_posteriors.R"))