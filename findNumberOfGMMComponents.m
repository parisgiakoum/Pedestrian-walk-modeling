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
minSubj = 1;    % starting subject
maxSubj = 11;   % ending subject
maxComp = 6;   % maximum components

for comp=1:maxComp % Components of GMM
    [GMModel, h] = fitGMMtoData(X(:,:, minSubj:maxSubj), comp, 'variables'); % Dt, meanF, len
    for subj=1:maxSubj - minSubj+1
        AIC(comp, subj)= GMModel{subj}.AIC;
        BIC(comp, subj)= GMModel{subj}.BIC;
    end
end

for j=1:maxSubj - minSubj+1
    [minAIC(j), numComponentsAIC(j)] = min(AIC(:, j));
    [minBIC(j), numComponentsBIC(j)] = min(BIC(:, j));
end

for comp=1:maxComp % Components of GMM
    totalAIC(comp) = sum(numComponentsAIC == comp);
    totalBIC(comp) = sum(numComponentsBIC == comp);
end

% Plot
figure;
for subj = 1:maxSubj - minSubj+1
    nexttile
    plot(BIC(:, subj), '.-', 'DisplayName','BIC')
    xlabel('components')
    ylabel('AIC/BIC')
    hold on;
    plot(AIC(:, subj), '.-', 'DisplayName','AIC')
    hold off;
    lgd = legend;
end

% Gradient Plot
figure;
for subj = 1:maxSubj - minSubj+1 
    nexttile
    plot(gradient(BIC(:, subj)), '.-', 'DisplayName','BIC')
    hold on;
    plot(gradient(AIC(:, subj)), '.-', 'DisplayName','AIC')
    hold off;
    lgd = legend;
end

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

%% Plot the AIC & BIC for a different number of components after normalisation of histograms
%% Dt, meanF, length

minSubj = 1;    % starting subject
maxSubj = 11;   % ending subject
maxComp = 6;   % maximum components

for comp=1:maxComp % Components of GMM
    for i=minSubj:maxSubj
        temp = X(:,1:3,i); %%%%%%
        GMModel{i} = fitgmdist (temp(any(~isnan(temp), 2), :), comp, 'Options', statset('MaxIter', 1500), 'SharedCovariance', true);
    end
    for subj=1:maxSubj - minSubj+1
        AIC(comp, subj)= GMModel{subj}.AIC;
        BIC(comp, subj)= GMModel{subj}.BIC;
    end
end

for j=1:maxSubj - minSubj+1
    [minAIC(j), numComponentsAIC(j)] = min(AIC(:, j));
    [minBIC(j), numComponentsBIC(j)] = min(BIC(:, j));
end

for comp=1:maxComp % Components of GMM
    totalAIC(comp) = sum(numComponentsAIC == comp);
    totalBIC(comp) = sum(numComponentsBIC == comp);
end

% Plot
figure;
for subj = 1:maxSubj - minSubj+1
    nexttile
    plot(BIC(:, subj), '.-', 'DisplayName','BIC')
    xlabel('components')
    ylabel('AIC/BIC')
    hold on;
    plot(AIC(:, subj), '.-', 'DisplayName','AIC')
    hold off;
    lgd = legend;
end

%% 2nd option: Put all steps from many experiments into one

% Xtotal is a matrix with all steps performed
temp = X(:,:,1);
Xtotal = temp(any(~isnan(temp), 2), :);

for i=2:size(X,3)
    temp = X(:,:,i);
    Xtotal = [Xtotal; temp(any(~isnan(temp), 2), :)];
end

maxComp = 6;   % maximum components

for comp=1:maxComp % Components of GMM
    GMModel = fitgmdist (Xtotal, comp, 'Options', statset('MaxIter', 1500), 'SharedCovariance', true);
    AIC(comp)= GMModel.AIC;
    BIC(comp)= GMModel.BIC;
end

[minAIC, numComponentsAIC] = min(AIC);
[minBIC, numComponentsBIC] = min(BIC);

figure;
plot(BIC, '.-', 'DisplayName','BIC')
xlabel('components')
ylabel('AIC/BIC')
hold on;
plot(AIC, '.-', 'DisplayName','AIC')
hold off;
lgd = legend;
