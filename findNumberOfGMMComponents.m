%%
% Find the optimal number of GMM components for the modelisation of experimental walk
%
% The lowest pointy in the score differential diagram indicates the
% selected number of components
%
% The problem is identified at: https://towardsdatascience.com/gaussian-mixture-model-clusterization-how-to-select-the-number-of-components-clusters-553bef45f6e4
%%%
%% Initialisation
%%
clear;
close all;
clc;

database=load('steps_database').database_passi; % load database

database = clearDb(database);   % Clear database

[time,force, x_coord, y_coord] = retrieveAllVariables(database);    %  Retrieve forces, times and coordinates

[X, ~, ~, ~, ~] = computeAllDesiredVariables(force, time, x_coord, y_coord);  % Extract Dt, meanf, len, angle

clear database time force x_coord y_coord;

%% Multivariate GMM

%% Single subject
subj = randi(125);    % random subject
maxComp = 15;   % maximum components
tries = 20;

Xtemp = X(:,1:3,subj);
Xtemp = Xtemp(any(~isnan(Xtemp), 2), :);

% Plot AIC/BIC for a variety of components and tries
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

% Plot mean index evaluation diagram
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

% Plot socre differentials
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

%% 2nd option: Put all steps into a unified walk table

% Xtotal is a matrix with all steps performed
temp = X(:,1:3,1);
Xtotal = temp(any(~isnan(temp), 2), :);
for i=2:size(X,3)
    temp = X(:,1:3,i);
    Xtotal = [Xtotal; temp(any(~isnan(temp), 2), :)];
end

maxComp = 15;   % maximum components
tries = 12; % number of modelisation tries

% Plot AIC/BIC for a variety of components and tries
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

% Plot mean index evaluation diagram
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

% Plot socre differentials
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
% Single subject
subj = randi(125);    % subject
maxComp = 15;   % maximum components
tries = 20;

% Plot AIC/BIC for a variety of components and tries
Xtemp = X(:,4,subj);
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

% Plot mean index evaluation diagram
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
