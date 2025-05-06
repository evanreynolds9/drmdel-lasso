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