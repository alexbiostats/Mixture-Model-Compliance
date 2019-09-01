Corresponding R code for "A Fully Bayesian Mixture Model Approach for Identifying Non-Compliance in a Regulatory Tobacco Clinical Trial" manuscript. Briefly, we propose estimating probability of compliance with a treatment given biomarker mixture models which either incorporate some biological relationship, estimate components without a proposed relationship, or consider a reversible-jump MCMC framework to choose between these two states per randomized dose-level/group.

GitHub_Functions_ComplianceMixture.R: a set of functions for implementing the mixture models, simulation studies, and figure creation

Simulation Results: a sub-folder that contains the simulation output used to create the tables and figures in the manuscript, includes README for folder

Simulation_ParallelComputing_Code.R: code to run the simulations in parallel cores

Simulation_Result_Code.R: code to take simulation results and format for paper

Example_Analysis_Code.R: code to simulate a sample data set with 5 dose-levels to implement the various models proposed in the manuscript
