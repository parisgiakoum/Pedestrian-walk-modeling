%% Find the optimal number of GMM components for the modelisation of experimental walk
%% Initialisation
%%
clear;
close all;
clc;

% Basic functions
% Clear junk, retrieve force-time-position measurements and find meanF-Dt-length-angle for each step of each subject

database=load('steps_database').database_passi; % load database

database = clearDb(database);   % Clear database

[time,force, x_coord, y_coord] = retrieveAllVariables(database);    %  Retrieve forces, times and coordinates

[X, ~, ~, ~, ~] = computeAllDesiredVariables(force, time, x_coord, y_coord);  % Extract Dt, meanf, len, angle

clear database time force x_coord y_coord;

%% Plot the normal AIC & BIC for a different number of components
%% Dt, meanF, length
% The problem is identified at: https://towardsdatascience.com/gaussian-mixture-model-clusterization-how-to-select-the-number-of-components-clusters-553bef45f6e4
subj = randi(125);    % subject
maxComp = 15;   % maximum components
tries = 20;

Xtemp = X(:,1:3,subj);
Xtemp = Xtemp(any(~isnan(Xtemp), 2), :);
figure;
t = tiledlayout(5,4,'TileSpacing','Compact');
for i=1:tries
    for comp=1:maxComp % Components of GMM
        GMModel = fitgmdist (Xtemp, comp, 'Options', statset('MaxIter', 1500), 'SharedCovariance', true);
        AIC(i, comp)= GMModel.AIC;
        BIC(i, comp)= GMModel.BIC;
    end
    
    nexttile;
    plot(BIC(i,:), '.-', 'DisplayName','BIC')
    hold on;
    plot(AIC(i,:), 'o-', 'DisplayName','AIC')
    hold off;
    lg = legend;
end
title(t,'GMMs fitted on a random pedestrian')
xlabel(t, 'Components')
ylabel(t, 'AIC/BIC')

meanBIC = sum(BIC,1)/tries;
meanAIC = sum(AIC,1)/tries;
figure;
plot(meanBIC, '.-', 'DisplayName','BIC')
title('Fitting a GMM on a random pedestrian - Mean histogram')
xlabel('Components')
ylabel('AIC/BIC')
hold on;
plot(meanAIC, 'o-', 'DisplayName','AIC')
hold off;
lgd = legend;

dBIC = diff(meanBIC);
dAIC = diff(meanAIC);
figure;
plot([dBIC(1) dBIC], '.-', 'DisplayName','BIC')
title('Fitting a GMM on a random pedestrian – Score intervals')
xlabel('Components')
ylabel('AIC/BIC interval (current to previous)')
hold on;
plot([dAIC(1) dAIC], 'o-', 'DisplayName','AIC')
hold off;
lgd = legend;

%% 2nd option: Put all steps from many experiments into one

% Xtotal is a matrix with all steps performed
temp = X(:,1:3,1);
Xtotal = temp(any(~isnan(temp), 2), :);
for i=2:size(X,3)
    temp = X(:,1:3,i);
    Xtotal = [Xtotal; temp(any(~isnan(temp), 2), :)];
end

Xtotal=zscore(Xtotal);

maxComp = 15;   % maximum components
tries = 12; % number of modelisation tries

figure;
t = tiledlayout(3,4,'TileSpacing','Compact');
for i=1:tries
    for comp=1:maxComp % Components of GMM
        GMModel = fitgmdist (Xtotal, comp, 'Options', statset('MaxIter', 1500), 'SharedCovariance', true);
        AIC(i,comp)= GMModel.AIC;
        BIC(i,comp)= GMModel.BIC;
    end
    
    nexttile;
    plot(BIC(i,:), '.-', 'DisplayName','BIC')
    hold on;
    plot(AIC(i,:), '.-', 'DisplayName','AIC')
    hold off;
    lg = legend;
end
title(t,'GMMs fitted on the unified table')
xlabel(t, 'Components')
ylabel(t, 'AIC/BIC')

meanBIC = sum(BIC,1)/tries;
meanAIC = sum(AIC,1)/tries;
figure;
plot(meanBIC, '.-', 'DisplayName','BIC')
title('Fitting a GMM on the unified table - Mean histogram')
xlabel('Components')
ylabel('AIC/BIC')
hold on;
plot(meanAIC, 'o-', 'DisplayName','AIC')
hold off;
lgd = legend;

dBIC = diff(meanBIC);
dAIC = diff(meanAIC);
figure;
plot([dBIC(1) dBIC], '.-', 'DisplayName','BIC')
title('Fitting a GMM on the unified table – Score intervals')
xlabel('Components')
ylabel('AIC/BIC interval (current to previous)')
hold on;
plot([dAIC(1) dAIC], 'o-', 'DisplayName','AIC')
hold off;
lgd = legend;

%% Angle
minSubj = 1;
maxSubj = 11;
maxComp = 8;
        
for comp=1:maxComp % Components of GMM
    [GMModel, h] = fitGMMtoData(X(:,:,minSubj:maxSubj), comp, 'angle'); % Dt, meanF, len
    for subj=1:maxSubj - minSubj+1
        AIC(comp, subj)= GMModel{subj}.AIC;
        BIC(comp, subj)= GMModel{subj}.BIC;
    end
end

for j=1:maxSubj-minSubj+1
    [minAIC(j), numComponentsAIC(j)] = min(AIC(:, j));
    [minBIC(j), numComponentsBIC(j)] = min(BIC(:, j));
end

for comp=1:maxComp % Components of GMM
    totalAIC(comp) = sum(numComponentsAIC == comp);
    totalBIC(comp) = sum(numComponentsBIC == comp);
end

% Plot
figure;
for subj = 1:maxSubj-minSubj+1
    nexttile
    plot(BIC(:, subj), '.-', 'DisplayName','BIC')
    hold on;
    plot(AIC(:, subj), '.-', 'DisplayName','AIC')
    hold off;
    lgd = legend;
end

% Gradient plot
figure;
for subj = 1:maxSubj-minSubj+1
    nexttile
    plot(gradient(BIC(:, subj)), '.-', 'DisplayName','BIC')
    hold on;
    plot(gradient(AIC(:, subj)), '.-', 'DisplayName','AIC')
    hold off;
    lgd = legend;
end



