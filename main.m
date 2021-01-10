%% Study and modeling of pedestrian walk with regard to the improvement of stability and comfort on walkways
% Pavlos Paris Giakoumakis, Technical University of Crete
% Implemented as a thesis project at University of Modena & Reggio Emilia
% Final version : 9-Jan-2021
%%%
% The implemented code has been tested on MATLAB R2019b only and
% requires the following files:
%%%
% main.m
% steps_database.mat
% clearDb.m
% retrieveAllVairables.m
% computeAllDesiredVariables.m
% fitGMMtoData.m
% mu_weight_statDescription.m
% sigmaStatDescription.m
% generateParameters.m
%
% findNumberOfGMMComponents.m -- Independent file used to estimate the
% number of components to use on GMMs
%%%
%% Abstract
% The static stability of footbridges or pedestrian walkways can be
% effectively assessed through several approaches developed in the fields
% of mechanical and civil engineering. On the other hand, the dynamic
% stability of pedestrian walkways represents an underexplored field and
% only in the last 2 years, the comfort of such structures has been
% investigated. A walkway under the tendency to oscillate, provokes panic
% and insecurity of the users and needs to be appropriately addressed in
% order to guarantee the safety of pedestrians.
%
% In this thesis, we introduce an innovative algorithm for modeling and
% simulation of human walk using Gaussian Mixture Models. Our model
% satisfies the requirements of simplicity, ease of use by engineers and
% is suitable to accurately assess the dynamic stability of walkways.
% Furthermore, we implement a simulator that can be used to provide
% reliable prediction and assessment of floor vibrations under human
% actions. Evaluation results are promising, showing that our simulator is
% capable of supplementing the experimental procedure in future research.
%%%
%% Initialisation

clear;
close all;
clc;

%% Basic functions
%% Clear junk, retrieve force-time-position measurements and find meanF-Dt-length-angle for each step of each subject

database=load('steps_database').database_passi; % load database

database = clearDb(database);   % Clear database

[time,force, x_coord, y_coord] = retrieveAllVariables(database);    %  Retrieve forces, times and coordinates

[X, Dt, meanF, len, angle] = computeAllDesiredVariables(force, time, x_coord, y_coord);  % Extract Dt, meanf, len, angle


%% Modeling of data
% Fit GMMs on each person
GMModel = fitGMMtoData(X, 3, 'variables'); % Dt, meanF, len
GMMAngle = fitGMMtoData(X, 2, 'angle');  % angle

% Statistical description of parameters
% Statistical description of mu and componentProportion
[nd_s_mu_table, nd_s_weight_table, s_mu_table{1}, s_weight_table] = mu_weight_statDescription(GMModel, 1); % mu for Dt and comonentProportion
[nd_s_mu_table_angle, nd_s_weight_table_angle, s_mu_table_angle, s_weight_table_angle] = mu_weight_statDescription(GMMAngle, 4);    % mu and componentProportion for angle
% mu for meanF, len
for variable=2:3
    [nd_s_mu_table(variable, :), ~, s_mu_table{variable}, ~] = mu_weight_statDescription(GMModel, variable);
end

% Statistical description of Sigma
[nd_Sigma, SigmaValues] = sigmaStatDescription(GMModel, 'variables');   % Dt, meanF, len
[nd_SigmaAngle, SigmaValuesAngle] = sigmaStatDescription(GMMAngle, 'angle');    % angle

%% Simulator
% Generate random parameters for final GMMs
[randomWeight, randomMu, randomSigma, randomWeightAngle, randomMuAngle, randomSigmaAngle] = generateParameters(nd_s_weight_table, nd_s_mu_table, nd_Sigma, nd_s_weight_table_angle, nd_s_mu_table_angle, nd_SigmaAngle);

% Extract a random walk
% Fit a GMM to the generated parameters
simulatedGMM = gmdistribution(randomMu, randomSigma, randomWeight); % Dt, meanF, len
simulatedGMMAngle = gmdistribution(randomMuAngle, randomSigmaAngle, randomWeightAngle); % angle
    
% Simulate a random walk
n = 30 ; % number of steps
randomWalk = random(simulatedGMM, n);   % Simulation of Dt, meanF, len
randomWalk(:,4) = random(simulatedGMMAngle, n) % Simulation of angle

% Calculate speed and mean speed
speed = randomWalk(:,3)./randomWalk(:,1);
meanSpeed = mean(speed).*(3.6) % km/h

%% Figures
%% Plot force-time measurements of first 10 steps of of a random subject i
figure;
i=randi([1,215]);   % random subject

plot(time(:,1:10,i),force(:,1:10,i))
xlabel('time');
ylabel('force');


%% GMM fitted on the interarrival time data of a random subject i
i=randi([1,215]);   % random subject

tempgm = gmdistribution(GMModel{i}.mu(:,1), GMModel{i}.Sigma(1), GMModel{i}.PComponents); % fit a GMM with parameters corresponding to Dt
figure;
temp = X(:,1,i);    % Dt data
histogram (temp(~isnan(temp)), 'BinWidth', 0.02, 'BinLimits',[0,1.6], 'normalization' , 'pdf' );
title('GMM fitted on Dt')
xlabel('Dt Data')
ylabel('Density')

% Curve of GMM
xgrid = linspace(0,1.6,1000)';
hold on; plot(xgrid,pdf(tempgm,xgrid),'r-'); hold off

%% GMM fitted on the mean force data of a random subject i
%i=randi([1,215]);   % random subject

tempgm = gmdistribution(GMModel{i}.mu(:,2), GMModel{i}.Sigma(2,2), GMModel{i}.PComponents); % fit a GMM with parameters corresponding to meanF
figure;
temp = X(:,2,i);    % meanF data
histogram (temp(~isnan(temp)), 'BinWidth', 2, 'BinLimits',[0 ,120], 'normalization' , 'pdf' );
title('GMM fitted on mean force')
xlabel('meanF Data')
ylabel('Density')

% Curve of GMM
xgrid = linspace(0,120,1000)';
hold on; plot(xgrid,pdf(tempgm,xgrid),'r-'); hold off

%% GMM fitted on the length of step data of a random subject
%i=randi([1,215]);   % random subject

tempgm = gmdistribution(GMModel{i}.mu(:,3), GMModel{i}.Sigma(3,3), GMModel{i}.PComponents); % fit a GMM with parameters corresponding to len
figure;
temp = X(:,3,i);    % len data
histogram (temp(~isnan(temp)), 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
title('GMM fitted on length of step')
xlabel('Length Data')
ylabel('Density')

% Curve of GMM
xgrid = linspace(0,1.7,1000)';
hold on; plot(xgrid,pdf(tempgm,xgrid),'r-'); hold off

%% GMM fitted on the step angle data of a random subject
% i=randi([1,215]);  % random subject

tempgm = gmdistribution(GMMAngle{i}.mu, GMMAngle{i}.Sigma, GMMAngle{i}.PComponents); % fit a GMM with parameters corresponding to len
figure;
temp = X(:,4,i);    % len data
histogram (temp(~isnan(temp)), 'BinWidth', 0.03, 'normalization' , 'pdf' );
title('GMM fitted on angle of step')
xlabel('Angle Data')
ylabel('Density')

% Curve of GMM
xgrid = linspace(-0.5,0.7,1000)';
hold on; plot(xgrid,pdf(tempgm,xgrid),'r-'); hold off

%% Gaussians fitted on mu of Dt, meanF, length
% Plot gmm and data
for variable = 1:3
    figure;
    for i=1:size(s_mu_table{variable},2)  % all components
        subplot(3,1,i)
        if (variable == 2)  % meanF
            temp_mu_table = s_mu_table{variable};
            histogram (temp_mu_table(:,i), 'BinWidth', 2, 'BinLimits',[0 ,120], 'normalization' , 'pdf' );
            title(strcat('GMM fitted on mu for component ', {' '}, num2str(i),' - mean force variable'));
            xlabel('mu data');
            ylabel('Density');
            
            % Curve of GMM
            xgrid = linspace(0,120,1000)';
            hold on; plot(xgrid,pdf(nd_s_mu_table{variable,i},xgrid),'r-'); hold off
        elseif  (variable == 3) % len
            temp_mu_table = s_mu_table{variable};
            histogram (temp_mu_table(:,i), 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
            title(strcat('GMM fitted on mu for component ', {' '} ,num2str(i),' - Length variable'));
            xlabel('mu data');
            ylabel('Density');
            xgrid = linspace(0,1.8,1000)';
            hold on; plot(xgrid,pdf(nd_s_mu_table{variable,i},xgrid),'r-'); hold off
        else    % Dt
            temp_mu_table = s_mu_table{variable};
            histogram (temp_mu_table(:,i), 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
            title(strcat('GMM fitted on mu for component ', {' '} ,num2str(i),' - Dt variable'));
            xlabel('mu data');
            ylabel('Density');
            xgrid = linspace(0,1.8,1000)';
            hold on; plot(xgrid,pdf(nd_s_mu_table{variable,i},xgrid),'r-'); hold off
        end
    end
end

%% Gaussians fitted on mixing probability (component proportion coefficients) of Dt, meanF, length
figure;
for i=1:size(s_weight_table,2)  % all components
    subplot(3,1,i)  
    histogram (s_weight_table(:,i), 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
    title(strcat('GMM fitted on mixing probability for component', {' '}, num2str(i)));
    xlabel('mixing probability');
    ylabel('Density');
    
    % Curve of GMM
    xgrid = linspace(0,1.8,1000)';
    hold on; plot(xgrid,pdf(nd_s_weight_table{i},xgrid),'r-'); hold off
end

%% Gaussians fitted on Sigma's critical values of Dt, meanF, length

figure;
for i=1:size(SigmaValues, 2) % all differentiated values

    subplot(3,2,i)
    % The seperation in the if clause is for aesthetic purposes on diagrams
    if i==1 | i==3 | i==5 
        histogram (SigmaValues(:,i), 'normalization' , 'pdf'  );
        title(strcat('GMM fitted on differentiated value', {' '}, num2str(i)))
        xlabel('Data');
        ylabel('Density');
        xgrid = linspace(-0.05,0.16,1000)';
        hold on; plot(xgrid,pdf(nd_Sigma{i},xgrid),'r-'); hold off
    elseif i==2
        histogram (SigmaValues(:,i), 'normalization' , 'pdf'  );
        title(strcat('GMM fitted on differentiated value', {' '}, num2str(i)))
        xlabel('Data');
        ylabel('Density');
        xgrid = linspace(-20,100,1000)';
        hold on; plot(xgrid,pdf(nd_Sigma{i},xgrid),'r-'); hold off
    else
        histogram (SigmaValues(:,i), 'normalization' , 'pdf'  );
        title(strcat('GMM fitted on differentiated value', {' '}, num2str(i)))
        xlabel('Data');
        ylabel('Density');
        xgrid = linspace(-2,2,1000)';
        hold on; plot(xgrid,pdf(nd_Sigma{i},xgrid),'r-'); hold off
    end    
end

%% Gaussians fitted on the angle data of a random subject
i=randi([1,215]);  % random subject

figure;
temp = X(:,4,i);    % angle data
histogram (temp(~isnan(temp)), 'BinWidth', 0.03, 'BinLimits',[-0.25,0.25], 'normalization' , 'pdf' );
title('GMM fitted on angle of step')
xlabel('Angle Data')
ylabel('Density')
ylim([0 11])

% Curve of GMM
xgrid = linspace(-0.25,0.25,1000)';
hold on; plot(xgrid,pdf(GMMAngle{i},xgrid),'r-'); hold off

%% Gaussians fitted on mu of angle
figure;
for i=1:size(s_mu_table_angle,2)  % all components
    subplot(2,1,i)
    histogram (s_mu_table_angle(:,i), 'BinWidth', 0.03, 'BinLimits',[-1.4,1.6], 'normalization' , 'pdf' );
    title(strcat('GMM fitted on mu for component', {' '} ,num2str(i),' - angle variable'));
    xlabel('mu data');
    ylabel('Density');
    xgrid = linspace(-1.4,1.6,1000)';
    hold on; plot(xgrid,pdf(nd_s_mu_table_angle{i},xgrid),'r-'); hold off
end

%% Gaussians fitted on mixing probability (component proportion coefficients) of angle
figure;
for i=1:size(s_weight_table_angle,2)  % all components
    subplot(2,1,i)
    histogram (s_weight_table_angle(:,i), 'BinWidth', 0.03, 'BinLimits',[0,1], 'normalization' , 'pdf' );
    title(strcat('GMM fitted on mixing probabilities for component', {' '} ,num2str(i),' - angle variable'));
    xlabel('mu data');
    ylabel('Density');
    xgrid = linspace(0,1,1000)';
    hold on; plot(xgrid,pdf(nd_s_weight_table_angle{i},xgrid),'r-'); hold off
end

%% Gaussians fitted on Sigma of angle
figure;
histogram (SigmaValuesAngle, 'BinWidth', 0.002, 'BinLimits',[0,0.04], 'normalization' , 'pdf' );
title(strcat('GMM fitted on Sigma - angle variable'));
xlabel('Sigma data');
ylabel('Density');
xgrid = linspace(-0,0.04,1000)';
hold on; plot(xgrid,pdf(nd_SigmaAngle,xgrid),'r-'); hold off

%% fig 34: Scattering a random walk
% Plot data acquired by simulation
dx = randomWalk(:,3).*cos(randomWalk(1,4));
dy = randomWalk(:,3).*sin(randomWalk(:,4));
for i=2:size(randomWalk , 1)    % all simulated steps
    dx(i) = dx(i-1)+ randomWalk(i,3).*cos(randomWalk(i,4));
end

figure;
plot(dx,dy,'--*','MarkerEdgeColor','k')
title('Scattering random walk');
xlabel('x\_coord');
ylabel('y\_coord');
%camroll(90)    % a different view
for k = 1: length (dx)
    text (dx (k)+0.07, dy (k)+0.01, num2str (k), 'Color','r')
end

%% fig. 35: Scatter first 10 steps of a subject from database
% Plot data from database to compare
clear xpos ypos
subj = 5; % subject of figure 9, 35 in the implementation report

for i=1:10 % number of steps to scatter
    if ~isempty(database{subj,i})
        xpos(i) = database{subj,i}.x_step;
        ypos(i) = database{subj,i}.y_step;
    end
end

figure;
plot(xpos,ypos,'--*','MarkerEdgeColor','k')
title('Scattering walk');
xlabel('x\_coord');
ylabel('y\_coord');

for k = 1: length (xpos)
    text (xpos (k), ypos (k)+0.01, num2str (k), 'Color','r')
end
