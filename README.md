The goal of this project is to build lasso related functionality into the existing `drmdel` package. It is inspired by the derivation of the Group Lasso Estimator for Density Ratio Models in my Masters thesis. Building this functionality directly into the C source code will greatly improve the efficiency of optimization functions.

Unlike the original project, this repository is not currently built as a standalone R package. Rather, it is a collection of programs that builds on functions from the original package to demonstrate the utility of the LASSO estimator in DRMs. As such, the documentation is sparse. Please view the original package, `drmdel`, for a more in-depth toolbox of DRM functions.

All credit for the original `drmdel` project goes to the author, Dr. Song Cai.

**Getting Started**

1. Ensure you have a C compiler installed on your machine.
2. To compile the C functions to a shared library, run `R CMD SHLIB drmdelLasso.c utilities.c basisFuncs.c -o <name>.dll` from the `src` folder in the terminal.
3. Within R session, navigate to the `src` folder and run `dyn.load("<name>.dll")` to make the C functions accessible.
4. In the same R session, navigate to the R directory and run `source("drmdelLasso.R")` to import the available R wrappers and functions to your session.
