
#####################################
# Title: Returns Based Style Analysis
# Author: Thomas Handscomb
#####################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Install packages and load into R session
#install.packages("quadprog")
#install.packages("matrix")

library(quadprog)
library(Matrix)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Specify the fund returns matrix over T = 10 time periods

T=10
R <- matrix(c(1.5, 1.6, 1.75, 1.70, 1.63, 1.59, 1.50, 1.59, 1.65, 1.70))

# n = 4 Index returns over T=10 time periods
n = 4

# The n style factor return streams
F1 <- c(0.76, 0.81, 0.89, 0.80, 0.81, 0.80, 0.79, 0.75, 0.69, 0.75)
F2 <- c(0.93, 0.45, 0.76, 0.45, 0.54, 0.64, 0.74, 0.56, 0.59, 0.45)
F3 <- c(0.15, 0.17, 0.21, 0.23, 0.26, 0.27, 0.28, 0.26, 0.23, 0.23)
F4 <- c(0.20, 0.25, 0.32, 0.35, 0.32, 0.30, 0.25, 0.18, 0.10, 0.09)

F <- cbind(F1, F2, F3, F4)
F

# Set up the symmetric, positive definite matrix D
e = matrix(1,T,1)

M <- diag(T) - ((e %*% t(e))/T)
M

eigen(M)$values

D <- (2/T)*t(F) %*% M %*% F

# Check that D is symmetric, positive semi-definite
D - (D+t(D))/2

## Set up the linear dvec coefficient
d <- 2*(((1/T) * t(R) %*% F) - (1/(T^2))*(t(e) %*% R %*% t(e) %*% F))
d

# Build up the constraint matrix A
L = c(rep(1,n))
L

A <- t(rbind(L, diag(n), -diag(n)))
A

# Define the constraint coefficient column vector b0
b0 <- matrix(c(1, rep(0,n), rep(-1,n)))
b0

# Solve the QP. Note the transpose taken for dvec and Amat
sol <- solve.QP(Dmat = D, dvec = d, Amat = A, bvec = b0, meq = 1)
sol$solution

# Double check that the sum of the beta coefficients = 1
runtot = 0
for (i in 1:n){
  runtot =+ runtot + sol$solution[i]
}
print (runtot)

# Construct the optimised solution from the solve.QP output
wvec  = {}
for (j in 1:n){
  wvec = rbind(wvec, sol$solution[j])
}
wvec

# Once the weight matrix has been constructed, the optimal solution can then be constructed
opsol <- F %*% Betavec
opsol

# View the return stream and the optimal solution
plot(R, ylim = c(0.0, 2.00))
points(opsol, col = 'blue')

# Calculate the fit of the optimal solution, recall the optimisation model is minimising var(R-F*wvec) 
var(R-opsol)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Try 1000 random solutions to check that you can't achieve a lower variance

# First create solution data frame with first row as the QP solved optimal solution
df_solutions = data.frame("Solution_ID" = character(0), "Variance" = numeric(0))
df_solutions = rbind(df_solutions, c(as.character("Optimal"), as.numeric(var(R-opsol))))
colnames(df_solutions) <- c("Solution_Id", "Variance")
#head(df_solutions)

# Then create 1000 random choices of weights that sum to 1 and 1000 corresponding random solutions
for (i in 1:1000){
  # Define random weights
  rancoef1 = runif(1, 0, 0.3)
  rancoef2 = runif(1, 0, 0.3)
  rancoef3 = runif(1, 0, 0.3)
  rancoef4 = 1- rancoef1 - rancoef2 - rancoef3
  random_wvec = c(rancoef1, rancoef2, rancoef3, rancoef4)
  
  # Construct the solution corresponding to the random weights
  randomsol = F %*% random_wvec
  
  # Create a temporary dataframe containing the iteration solution
  df_temp = data.frame("Solution_ID" = character(0), "Variance" = numeric(0))
  df_temp = rbind(df_temp, c(paste0("RanSol_",i), var(R-randomsol)))
  colnames(df_temp) <- c("Solution_Id", "Variance")
  
  # Append iteration solution to the total solution dataframe
  df_solutions = rbind(df_solutions, df_temp)
}
# The rbind has created Variance as a factor...
sapply(df_solutions, class)
#...so convert back to numeric
df_solutions$Variance = as.numeric(paste(df_solutions$Variance))
sapply(df_solutions, class)

# Sort the solution dataframe by descending Variance, show the lowest 10 only
head(df_solutions[order(df_solutions$Variance),], 10)

# Optimal solution is still the top!
