# Run simulations

# Call init.R to setup shared library and all functions
# This should be run with the working directory as the project root
source("R\\init.R")

# Set parameter values
distribution = "normal"
n = 250
model = 12
lambdaVals = c(0, (n^(1/3)*(3))^(-4:4))
adaptive = TRUE
runs = 10

# Compute d based on model
if(model %in% 1:4){
  d = 1
}else if(model %in% 5:6){
  d = 2
}else if(model %in% 7:10){
  d = 3
}else if(model == 11){
  d = 4
}else if(model ==12){
  d = 5
}else{
  stop("Invalid model provided - adjust the program header.")
}

# Run simulation
simData = runSimulation(distribution, n, d, model, lambdaVals, adaptive = adaptive, runs = runs)

# Setwd back to data folder
setwd("..")
setwd("Data")

# Write data
simDataDf = as.data.frame(simData)
if(adaptive){
  adapStr = "adap"
} else{
  adapStr = "reg"
}
fileName = paste0(paste("sim_results", adapStr, n, sep = "_"), ".csv")
write.csv(simDataDf, fileName, row.names = FALSE)
