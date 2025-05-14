# Initialize the shared C library, R Scripts for drmdellasso
# This is assumed to be called from the directory root
# Changes to relative paths must be made if this is not the case

# Load .Renviron file to get shared lib name
readRenviron(".Renviron")
shared_lib = Sys.getenv("SHARED_LIB")

# Build the shared library
setwd("src")
lib_str = paste0(shared_lib, ".dll")
command_str = paste("R CMD SHLIB -o",lib_str,"drmdelLasso.c utilities.c basisFuncs.c")
system(command_str)

# Load the shared library
dyn.load(lib_str)

# Load the R wrappers from the R folder
source("..\\R\\drmdelLasso.R")

# Define a wrapper function to run simulations for the default distributions and parameter values
# These are gamma and normal distributions
runSimulation = function(distribution, n, d, model, lambdaVals, adaptive = FALSE, runs = 1000){
  # x ,n_total, n_samples, m, d, model, lambda_vals, adaptive = FALSE
  
  # Ensure distribution is normal or gamma
  if(!(distribution %in% c("normal", "gamma"))){
    stop('Invalid distribution string passed. Currently only "normal" and "gamma" are supported.')
  }
  
  # Define parameters for the distributions
  if(distribution == "normal"){
    mu_0 = 2
    mu_1 = 1.8
    mu_2 = 1.9
    sigma_0 = 2
    sigma_1 = 0.89
    sigma_2 = 0.875
  } else if (distribution == "gamma"){
    rate_0 = 2
    rate_1 = 1.4
    rate_2 = 1.2
    shape_0 = 1.8
    shape_1 = 1.2
    shape_2 = 1
  }
  
  # Set the seed based on the sample size
  if(n == 250){
    set.seed(20)
  } else if (n == 500){
    set.seed(21)
  } else if (n == 1000){
    set.seed(22)
  } else if (n == 2500){
    set.seed(23)
  } else if (n == 5000){
    set.seed(24)
  } else {
    stop("An invalid sample size (individual) was provided. Currently only 250, 500, 1000, 2500 and 5000 are supported.")
  }
  
  # Set m - always setting to two for these simulation
  m = 2
  
  # Set n_total and n_samples
  n_samples = rep(n, m+1)
  n_total = sum(n_samples)
  
  # Compute path length, noting that 0 will be added if not passed
  # It will be added in the solution path function if not
  if(0 %in% lambdaVals){
    pathLength = length(lambdaVals)
  } else{
    pathLength = length(lambdaVals)+1
  }
  
  # Create matrix for simulation results
  pathLength = length(lambdaVals)
  simulationResults = matrix(0, nrow = runs*pathLength, ncol = 2*(d+1) + 6)
  
  # Populate matrix by generating solution paths
  for(i in 1:runs){
    # Simulate data
    if(distribution == "normal"){
      x0_test = rnorm(n,mu_0,sigma_0)
      x1_test = rnorm(n,mu_1,sigma_1)
      x2_test = rnorm(n,mu_2,sigma_2)
    } else{ # The distribution is gamma
      x0_test = rgamma(n,shape=shape_0,rate=rate_0)
      x1_test = rgamma(n,shape=shape_1,rate=rate_1)
      x2_test = rgamma(n,shape=shape_2,rate=rate_2)
    }
    x_test = c(x0_test,x1_test,x2_test)
    
    # compute solution path
    solPath = solutionPath(x = x_test, n_total = n_total, n_samples = n_samples, 
                           m = m, d = d, model = model, lambdaVals, adaptive = adaptive)
    
    # Set the values of simulationResults
    simulationResults[(1+pathLength*(i-1)):(pathLength*i), 1] = i # Set first column to run iteration
    # Set the other values to solPath
    simulationResults[(1+pathLength*(i-1)):(pathLength*i), 2:(2*(d+1) + 6)] = solPath # essentially using m = 2
    
    # If we are on the first run, pull the column names from solPath
    colnames(simulationResults) = c("Run", colnames(solPath))
  }
  
  # Return simulationResults
  return(simulationResults)
}