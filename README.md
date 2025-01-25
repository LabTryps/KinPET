# KinPET - Kinetic Parameter Estimation Toolkit

## Overview

This project aims at providing a simple, versatile and robust tool for users interested in obtaining kinetic parameters for enzymatic reaction models. Simple, as it provides unique and standardized interfaces to access complex functionalities, such as simulating, fitting a model to data and obtaining statistics on its parameters. Versatile, because it adapts to any model and collection of experiments that the user has, providing freedom of choice and unlimited experimentation. Robust, because it is able to deal with systems of extreme complexity, such as complete pathway models, and still provide useful results to inform the user about the uncertainty of each parameter and assist in the design of complementary experiments.

The following files are provided, each consisting of a method of the same name, which the user can import and call in their own code:
* **kinetics.jl** – simulates progress curves from initial conditions, a model, and its parameters
* **fitmodel.jl** – performs parameter adjustment, given a model and an experimental dataset. Minimizes mean squared error between simulated and experimental curves
* **bootstrap.jl** – gets the parameters for several sets of experiments, randomly selected from the available ones. On top of these values, mean and standard deviation can be calculated for each parameter, as well as generate histograms of their distribution

In addition, a *models.jl* file is provided, with a collection of example models, as well as validation scripts, exemplifying the use of the program in various situations and verifying its correct operation.

## Structure, interfaces and modules

The diagram below illustrates the program's modules and their connections:

![image](https://github.com/user-attachments/assets/7fe78e55-372c-4aa2-a8be-65c276f59547)


Above is presented the general structure of the tool: each blue block is a file in Julia containing a single method (bootstrap, fitmodel and kinetics). At the top are conceptually represented the inputs, and at the bottom the outputs of each method. The arrows indicate the order of execution of each step and calls from one method to another. The user can call each of the three illustrated methods independently, depending on the task to be performed (simulation of the model or adjustment of its parameters).

For a detailed description of each module, check out our wiki!
