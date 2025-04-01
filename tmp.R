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
sum_stats <- simSumStats(ddm_parameters, N)

full_data <- cbind(design, sum_stats)

head(full_data)

# --- Test Plotting Interactively ---


#pdf("tmp.pdf")

# Set up plot parameters (optional, adjust if needed)
layout(matrix(c(1), nrow = 1, ncol = 1))
par(mai = c(1, 1, 0.5, 0.1), oma = c(0, 0, 0, 0)) 

# Create an empty plot frame
plot(NULL, xlim = c(0.5, 2.5), 
     ylim = range(full_data$Mrt, na.rm = TRUE) * c(0.95, 1.05), 
     xlab = "Factor A",
     ylab = "Mean Reaction Time (Mrt)",
     main = "Simulated Mean RTs by Factor Combination (Interactive Test)",
     xaxt = "n", las = 1)     

# Add custom x-axis labels
axis(1, at = c(1, 2), labels = c("A1", "A2"))

# Define colors for Factor B levels
colors_b <- c("red", "blue")
names(colors_b) <- c("B1", "B2")

# Add jittered points for each combination
jitter_factor <- 0.1 

# Points for Factor B1
points(jitter(as.numeric(factor(full_data$Factor.A[full_data$Factor.B == "B1"])), factor = jitter_factor),
       full_data$Mrt[full_data$Factor.B == "B1"],
       pch = 16,
       col = colors_b["B1"])

# Points for Factor B2
points(jitter(as.numeric(factor(full_data$Factor.A[full_data$Factor.B == "B2"])), factor = jitter_factor),
       full_data$Mrt[full_data$Factor.B == "B2"],
       pch = 16,
       col = colors_b["B2"])

# Add a legend
legend("topright",
       legend = c("Factor B1", "Factor B2"),
       col = colors_b,
       pch = 16,
       bty = "n")

# --- End of Interactive Test ---
#dev.off()

# --- End of Saving Plot ---