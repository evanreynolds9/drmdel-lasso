The goal of this project is to build lasso related functionality into the existing `drmdel` package. It is inspired by the derivation of the Group Lasso Estimator for Density Ratio Models in my Masters thesis. Building this functionality directly into the C source code will greatly improve the efficiency of optimization functions.

Unlike the original project, this repository is not currently built as a standalone R package. Rather, it is a collection of programs that builds on functions from the original package to demonstrate the utility of the LASSO estimator in DRMs. As such, the documentation is sparse. Please view the original package, `drmdel`, for a more in-depth toolbox of DRM functions.

All credit for the original `drmdel` project goes to the author, Dr. Song Cai.

**Getting Started**

1. Ensure you have a C compiler installed on your machine.
2. Create a .Renviron file in root directory, and set the name you choose for your shared library to the following variable: `SHARED_LIB`
3. In your R session, from the project root, run the `init.R` file to build the shared library and load all functions. This can be done with `source("R\\init.R")`. 