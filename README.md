# Study and modeling of pedestrian walk with regard to the improvement of stability and comfort on walkways
Diploma Thesis - Study and modeling of pedestrian walk with regard to the improvement of stability and comfort on walkways

A part of this project was implemented in Department of Sciences and Methods for Engineering, University of Modena & Reggio Emilia (Italy) in cooperation with Technical University of Crete (Greece).

# Committee:
- Project supervisor: Assistant Professor Fabrizio Pancaldi, Assistant Professor at University of Modena and Reggio Emilia (UNIMORE)

- Diploma thesis supervisor: Professor Zervakis Michalis, Vice Rector and Professor at Technical University of Crete

- Member of the committee: Professor George Karystinos, Professor at Technical University of Crete

# Abstract:
The static stability of footbridges or pedestrian walkways can be effectively assessed through several approaches developed in the fields of mechanical and civil engineering. On the other hand, the dynamic stability of pedestrian walkways represents an underexplored field and only in the last 2 years, the comfort of such structures has been investigated. A walkway under the tendency to oscillate, provokes panic and insecurity of the users and needs to be appropriately addressed in order to guarantee the safety of pedestrians.

In this thesis, we introduce an innovative algorithm for modeling and simulation of human walk using Gaussian Mixture Models. Our model satisfies the requirements of simplicity, ease of use by engineers and is suitable to accurately assess the dynamic stability of walkways. Furthermore, we implement a simulator that can be used to provide reliable prediction and assessment of floor vibrations under human actions. Evaluation results are promising, showing that our simulator is capable of supplementing the experimental procedure in future research.

# The repository contains:
- Pedestrian_walk_modeling.pdf: The diploma thesis report - Read this for more information.

- Implemented code and important files.

# Note:
MATLAB R2019b Update 4 was used for the development and testing of the implemented code. MATLAB is also required for the execution.

The final version of this thesis was submitted at 9 Jan. 2021.

A state of the art paper on the subject is in work.

# Outline of the implemented code:
The implemented code for all the functions described in the implementation report requiere the following files:
- main.m: The main file - it connects all the functions to implement all procedures and MATLAB figures included in the report.
- steps_database.mat: The database used for the modelisation.
- clearDb.m: The function clearDb that clears trash in database & fills force-time matrixes to match indices.
- retrieveAllVariables.m: The function retrieveAllVariables that retrieves all forces, times and coordinates from the database for each step.
- computeAllDesiredVariables.m: The function computeAllDesiredVariables that computes the interarrival time between one step and the next, mean force induced, length and angle for each step, as well as clears any wrong values.
- fitGMMtoData.m: The function fitGMMtoData that fits a gaussian mixture model to each subject's data.
- mu_weight_statDescription.m: The function mu_weight_statDescription that creates a statistical description for mu (mean values) and mixing proportions of mixture components.
- sigmaStatDescription.m: The function sigmaStatDescription that create a statistical description for Sigma (differentiated covariance values).
- generateParameters.m: The function generateParameters that extracts a random set of parameters to use in the GMM of the simulator.
- findNumberOfGMMComponents.m - Independent MATLAB file used to estimate the number of components to use on the GMM fittings.

All files contain extensive comments in order to be easily understandable. A run of the main file will produce all the tables acquired by the database and the variable models, the parameter models, the simulation and finally, the figures included in the report.

The table randomWalk describes the simulated random walk.
