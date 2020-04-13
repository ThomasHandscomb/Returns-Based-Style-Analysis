
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

#eigen(M)$values

D <- (2/T)*t(F) %*% M %*% F

# Check that D is symmetric, positive semi-definite
D - (D+t(D))/2

## Set up dvec
d <- 2*(((1/T) * t(R) %*% F) - (1/(T^2))*(t(e) %*% R %*% t(e) %*% F))
d

# Build up the constraint matrix A and bvec, b0
A <- t(rbind(rep(1,n), diag(n), -diag(n)))
A

# Define the constraint coefficient column vector b0
b0 <- matrix(c(1, rep(0,n), rep(-1,n)))
b0

# Call the QP solver
sol <- solve.QP(Dmat = D, dvec = t(d), Amat = A, bvec = b0, meq = 1)

# The solution gives the vector w comprised of optimal weights
sol$solution

#> sol$solution
#[1] 0.6679001 0.0000000 0.3320999 0.0000000

# Double check that the sum of the w_i coefficients = 1
runtot = 0
for (i in 1:n){
  runtot =+ runtot + sol$solution[i]
}
print (runtot)
#[1] 1

# Construct the optimised solution from the solve.QP output
wvec  = {}
for (j in 1:n){
  wvec = rbind(wvec, sol$solution[j])
}
wvec

# Once the weight matrix has been constructed, the optimal solution can then be constructed
opsol <- F %*% wvec
opsol

# Calculate the fit of the optimal solution, recall the optimisation model is minimising var(R-F*wvec) 
var(R-opsol) #0.006430116

# View a plot of the return stream and the optimal solution
plot(R, ylim = c(0.0, 2.50), type = 'b'
     , main="Return stream vs. Optimal solution"
     , xlab="Time"
     , ylab="Return")
points(opsol, col = 'blue', type = 'b')
legend("topleft", legend=c("Investment", "Optimal Solution")
       , col=c("black", "blue")
       , lty=c(1,1))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Try 10000 random solutions to check that you can't achieve a lower variance

# First create solution data frame with first row as the QP solved optimal solution
df_solutions = data.frame("Solution_ID" = character(0), "Variance" = numeric(0))
df_solutions = rbind(df_solutions, c(as.character("Optimal"), as.numeric(var(R-opsol))))
colnames(df_solutions) <- c("Solution_Id", "Variance")
#head(df_solutions)

# Then create 10000 random choices of weights that sum to 1 and corresponding random solutions
for (i in 1:10000){
  # Define random weights
  R4 <- runif(4)
  R4 <- R4/sum(R4)
  rancoef1 = R4[1]
  rancoef2 = R4[2]
  rancoef3 = R4[3]
  rancoef4 = R4[4]
  random_wvec = as.matrix(c(rancoef1, rancoef2, rancoef3, rancoef4))
  
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

#Solution_Id    Variance
#1        Optimal 0.006430116
#8607 RanSol_8606 0.006632959
#4117 RanSol_4116 0.006643589
#509   RanSol_508 0.006652932
#8542 RanSol_8541 0.006659689
#6736 RanSol_6735 0.006747747
#7805 RanSol_7804 0.006753655
#4066 RanSol_4065 0.006755018
#9213 RanSol_9212 0.006757928
#7053 RanSol_7052 0.006819245

# Optimal solution is still the top!
