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

%% run seperately after normalization
for subj=1:215
    temp = X(:,1:3,subj);
    temp = temp(any(~isnan(temp), 2), :);
    z=zscore(temp);

    maxComp = 15;   % maximum components
    for comp=1:maxComp % Components of GMM
        GMModel{subj} = fitgmdist (z, comp, 'Options', statset('MaxIter', 1500), 'SharedCovariance',true);
        AIC(subj,comp)= GMModel{subj}.AIC;
        BIC(subj,comp)= GMModel{subj}.BIC;
    end
    [minAIC(subj), numComponentsAIC(subj)] = min(AIC(subj, :));
    [minBIC(subj), numComponentsBIC(subj)] = min(BIC(subj,:));
end

%% Option 1
% THELEI ZSCORE I OXI???

temp = X(:,1:3,1);
temp = temp(any(~isnan(temp), 2), :);
Xadded = temp(1:18,:);
for i=2:215
   temp =X(:,1:3,i);
   temp = temp(any(~isnan(temp), 2), :);
   Xadded = Xadded + temp(1:18,:); 
end

Xadded = zscore(Xadded);
Xadded = Xadded/215;

maxComp = 8;   % maximum components
for comp=1:maxComp % Components of GMM
    GMModel = fitgmdist (Xadded, comp, 'Options', statset('MaxIter', 1500), 'SharedCovariance',true);
    AIC(comp)= GMModel.AIC;
    BIC(comp)= GMModel.BIC;
end
[minAIC, numComponentsAIC] = min(AIC)
[minBIC, numComponentsBIC] = min(BIC)

%% Option 2: Put all steps from all experiments into one
% Xtotal is a matrix with all steps performed
temp = X(:,1:3,1);
Xtotal = temp(any(~isnan(temp), 2), :);
for i=2:size(X,3)
    temp = X(:,1:3,i);
    Xtotal = [Xtotal; temp(any(~isnan(temp), 2), :)];
end

z=zscore(Xtotal);
n=normalize(Xtotal,'norm');
maxComp = 15;   % maximum components
for comp=1:maxComp % Components of GMM
    GMModel = fitgmdist (n, comp, 'Options', statset('MaxIter', 1500), 'SharedCovariance',true);
    AIC(comp)= GMModel.AIC;
    BIC(comp)= GMModel.BIC;
end
[minAIC, numComponentsAIC] = min(AIC)
[minBIC, numComponentsBIC] = min(BIC)


[r,~] = find(z>3)

[r,~] = find(z>3)
z(r,:) = []
[r,~] = find(z<-1.5)
z(r,:) = []

plot(z(:,1))
figure;
plot(z(:,2))
figure;
plot(z(:,3))



