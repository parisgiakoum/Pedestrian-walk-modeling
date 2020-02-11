%% 25 - Plots of GMMs fitted on statistical description of Sigma for angles
%% Clear junk, retrieve force-time-position measurements and find meanF-Dt-length-angle for each step of each subject

database=load('steps_database').database_passi;

database = clearDb(database);

[time,force, x_coord, y_coord] = retrieveAllVariables(database);

[X, Dt,meanF, len, angle] = computeAllDesiredVariables(force, time, x_coord, y_coord);

%% GMMs for each subject
[GMMAngle, h_angle] = fitGMMtoData(X, 3, 'angle');

%% GMMs for Sigma
[GMModelSigmaAngle, SigmaValuesAngle] = sigmaStatDescription(GMMAngle, 'angle');

%% 25 - Sigma of angle
figure;
histogram (SigmaValuesAngle, 'BinWidth', 0.002, 'BinLimits',[0,0.04], 'normalization' , 'pdf' );
title(strcat('GMM fitted on Sigma - angle variable'));
xlabel('Sigma data');
ylabel('Density');
xgrid = linspace(0,0.04,1000)';
hold on; plot(xgrid,pdf(GMModelSigmaAngle,xgrid),'r-'); hold off

