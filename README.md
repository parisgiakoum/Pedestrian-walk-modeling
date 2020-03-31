# Study, modeling and simulation of pedestrian walk with regard to the improvement of stability and comfort on walkways
Diploma Thesis - Study, modeling and simulation of pedestrian walk with regard to the improvement of stability and comfort on walkways

This project was implemented in Department of Sciences and Methods for Engineering, University of Modena & Reggio Emilia (Italy) in cooperation with Technical University of Crete (Greece).

Project supervisor: Dr. Fabrizio Pancaldi, Assistant Professor at University of Modena and Reggio Emilia
Diploma thesis supervisor: Dr. Zervakis Michalis, Vice Rector and Professor at Technical University of Crete

The completion of the thesis is still remaining due to the COVID-19 pandemic situation in Italy.
A state of the art paper on the subject will soon be published.

The repository contains:
Impelementation_report file: Contains an extended implementation report - Read this for more information
Referenced Papers and Thesis projects file: References
Implemented code

MATLAB R2019b Update 4 was used for the development and testing of the implemented code. MATLAB is also required for the execution

# Outline of the implemented code:
The implemented code for all the functions described in the implementation report requiere the following files:
- main.m: The main file - it connects all the functions to implement all procedures and MATLAB figures included in this report.
- steps_database.mat: The database used for the modelisation.
- clearDb.m: The function clearDb that clears trash in database & fills force-time matrixes to match indices.
- retrieveAllVariables.m: The function retrieveAllVariables that retrieves all forces, times and coordinates from the database for each step.
- computeAllDesiredVariables.m: The function computeAllDesiredVariables  that computes the interarrival time between one step and the next, mean force induced, length and angle for each step, as well as clears any wrong values.
- fitGMMtoData.m: The function fitGMMtoData that fits a gaussian mixture model to each subject's data and performs the two-sample Kolmogorov-Smirnov tests to determine if the distributions describe the data.
- mu_weight_statDescription.m: The function mu_weight_statDescription that creates a statistical description for mu (mean values) and mixing proportions of mixture components.
- sigmaStatDescription.m: The function sigmaStatDescription that create a statistical description for Sigma (differentiated covariance values)
- generateParameters.m: The function generateParameters that extracts a random set of parameters to use in the GMM of the simulator

All those files contain extensive comments in order to be easily understandable. A run of the main file will produce all the tables acquired by the database and the fitted GMMs for each person, the GMMs of the parameters describing the initial GMMs, the final GMM and finally, the figures included in the report.

The table randomWalk will contain all the important variables of all n steps of the simulated random walk.
