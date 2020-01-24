%% Init
clear all;
close all;
clc;

%%
%% Clear junk, retrieve force-time-position measurements and find meanF-Dt-length-angle for each step of each subject

database=load('steps_database').database_passi;

database = clearDb(database);

[time,force, x_coord, y_coord] = retrieveAllVariables(database);

[X, Dt,meanF, len, angle] = computeAllDesiredVariables(force, time, x_coord, y_coord);

%% GMMs for each subject
[GMModel, h] = fitGMMtoData(X, 5, 'variables');
[GMMAngle, h_angle] = fitGMMtoData(X, 3, 'angle');

%% Simulator
%% Generate random parameters for Dt, meanF, length
%%% Random Sigma
[GMModelSigma, SigmaValues] = sigmaStatDescription(GMModel, 'variables');

for i=1:length(GMModelSigma)    % 1-6
    randomSigmaValues(i) = random(GMModelSigma{i},1);
end

% 1-3 is the diagonal, 4 is 1-2, 5 is 1-3 and 6 is 2-3
randomSigma = diag(randomSigmaValues(1:3));
randomSigma(1,2) = randomSigmaValues(4); randomSigma(2,1) = randomSigmaValues(4);
randomSigma(1,3) = randomSigmaValues(5); randomSigma(3,1) = randomSigmaValues(5);
randomSigma(2,3) = randomSigmaValues(6); randomSigma(3,2) = randomSigmaValues(6);


%%% Random w, mu
[GM_s_mu_mat, GM_s_weight_mat, s_mu_mat, s_weight_mat] = mu_weight_statDescription(GMModel, 1);

for i=1:length(GM_s_weight_mat) % 1-5
    randomWeight(i) = random(GM_s_weight_mat{i},1);
    randomMuValues(i) = random(GM_s_mu_mat{i},1);
end

% Weight needs to sum to 1, so the values are normalised
randomWeight = randomWeight./sum(randomWeight);

%%% Random mu
randomMu(:, 1) = randomMuValues.';
for variable=2:3
    [GM_s_mu_mat, GM_s_weight_mat, s_mu_mat, s_weight_mat] = mu_weight_statDescription(GMModel, variable);
    for i=1:length(GM_s_mu_mat) % 1-5
        randomMuValues(i) = random(GM_s_mu_mat{i},1);
    end
    randomMu(:, variable) = randomMuValues.';
end

%% Generate random parameters for angle
%%% Random Sigma
[GMModelSigmaAngle, SigmaValuesAngle] = sigmaStatDescription(GMMAngle, 'angle');
randomSigmaAngle = random(GMModelSigmaAngle, 1);

%%% Random w, mu
[GM_s_mu_mat_angle, GM_s_weight_mat_angle, s_mu_mat_angle, s_weight_mat_angle] = mu_weight_statDescription(GMMAngle, 4);
for i=1:length(GM_s_weight_mat_angle) % 1-3
    randomWeightAngle(i) = random(GM_s_weight_mat_angle{i},1);
    randomMuAngle(i) = random(GM_s_mu_mat_angle{i},1);
end
randomWeightAngle = randomWeightAngle./sum(randomWeightAngle);
randomMuAngle = randomMuAngle.';

%% Extract a random walk
% Fit a GMM to the generated parameters
simulatedGMM = gmdistribution(randomMu, randomSigma, randomWeight);
simulatedGMMAngle = gmdistribution(randomMuAngle, randomSigmaAngle, randomWeightAngle);

% Simulate a random walk
RandomWalk = random(simulatedGMM, 50);
RandomWalk(:,4) = random(simulatedGMMAngle, 50);

% Calculate velocity and mean velocity
velocity = RandomWalk(:,3)./RandomWalk(:,1);
meanVel = mean(velocity).*(3.6); % km/h

%% Scatter the random walk
x0 = 0;
y0 = 0;

scatter(x0,y0,'.')

%% Plotting GMM and data
%% Time
i=1;

tempgm = gmdistribution(GMModel{i}.mu(:,1), GMModel{i}.Sigma(1), GMModel{i}.PComponents);
figure;
temp = X(:,1,i);
histogram (temp(~isnan(temp)), 'BinWidth', 0.02, 'BinLimits',[0,1.6], 'normalization' , 'pdf' );
xgrid = linspace(0,1.6,1000)';
hold on; plot(xgrid,pdf(tempgm,xgrid),'r-'); hold off

%% Force
i=1;

tempgm = gmdistribution(GMModel{i}.mu(:,2), GMModel{i}.Sigma(2,2), GMModel{i}.PComponents);
figure;
temp = X(:,2,i);
histogram (temp(~isnan(temp)), 'BinWidth', 2, 'BinLimits',[0 ,120], 'normalization' , 'pdf' );
xgrid = linspace(0,120,1000)';
hold on; plot(xgrid,pdf(tempgm,xgrid),'r-'); hold off

%% Length
i=5;

tempgm = gmdistribution(GMModel{i}.mu(:,3), GMModel{i}.Sigma(3,3), GMModel{i}.PComponents);
figure;
temp = X(:,3,i);
histogram (temp(~isnan(temp)), 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
xgrid = linspace(0,1.7,1000)';
hold on; plot(xgrid,pdf(tempgm,xgrid),'r-'); hold off


%% Statistical description of mu, w, sigma & Plot GMM and data
%% Sigma
% 1-3 is the diagonal, 4 is 1-2, 5 is 1-3 and 6 is 2-3
[GMModelSigma, SigmaValues] = sigmaStatDescription(GMModel);

for i=1:size(SigmaValues, 2) % 1-6
%   rand_sigma_val(:,i) = random(GMModelSigma{i},415);
% 	h_sigma(i) = kstest2(SigmaValues(:,i),rand_sigma_val(:,i))

    figure;
    if i==1 | i==3 | i==5
        histogram (SigmaValues(:,i), 'normalization' , 'pdf'  );
        xgrid = linspace(-0.1,0.08,1000)';
        hold on; plot(xgrid,pdf(GMModelSigma{i},xgrid),'r-'); hold off
    elseif i==2
        histogram (SigmaValues(:,i), 'normalization' , 'pdf'  );
        xgrid = linspace(-0.1,20,1000)';
        hold on; plot(xgrid,pdf(GMModelSigma{i},xgrid),'r-'); hold off
    else
        histogram (SigmaValues(:,i), 'normalization' , 'pdf'  );
        xgrid = linspace(-0.35,0.5,1000)';
        hold on; plot(xgrid,pdf(GMModelSigma{i},xgrid),'r-'); hold off
    end
    
end

%% mu & weights
% 1st approach

variable = 1;
% Create vectors
all_mu = GMModel{1}.mu(:,variable);
all_weights =  GMModel{1}.ComponentProportion.';
for i=2:length(GMModel) %2-215
    all_mu = [all_mu ; GMModel{i}.mu(:,variable)];
    all_weights = [all_weights ; GMModel{i}.ComponentProportion.'];
end

% Fit GMM to mu and weight arrays
GM_all_mu = fitgmdist (all_mu, 3);
% rand_all_mu = random(GM_all_mu,415);
% h_all_mu = kstest2(all_mu,rand_all_mu)

GM_all_weights = fitgmdist (all_weights,3);
% rand_all_weights = random(GM_all_weights,415);
% h_all_weights = kstest2(all_weights,rand_all_weights)

% Plot gmm and data
figure;
if (variable == 2)
    histogram (all_mu, 'BinWidth', 2, 'BinLimits',[0 ,120], 'normalization' , 'pdf'  );
    xgrid = linspace(0,120,1000)';
    hold on; plot(xgrid,pdf(GM_all_mu,xgrid),'r-'); hold off
    figure;
    histogram (all_weights, 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf'  );
    xgrid = linspace(-0.1,1.8,1000)';
    hold on; plot(xgrid,pdf(GM_all_weights,xgrid),'r-'); hold off
else
    histogram (all_mu, 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
    xgrid = linspace(0,1.8,1000)';
    hold on; plot(xgrid,pdf(GM_all_mu,xgrid),'r-'); hold off
    figure;
    histogram (all_weights, 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
    xgrid = linspace(-0.1,1.8,1000)';
    hold on; plot(xgrid,pdf(GM_all_weights,xgrid),'r-'); hold off
end

%% 2nd approach
variable = 1;
[GM_s_mu_mat, GM_s_weight_mat, s_mu_mat, s_weight_mat] = mu_weight_statDescription(GMModel, variable);

% Plot gmm and data
for i=1:size(s_mu_mat,2)  %1-5
    figure;
    if (variable == 2)
        histogram (s_mu_mat(:,i), 'BinWidth', 2, 'BinLimits',[0 ,120], 'normalization' , 'pdf' );
        xgrid = linspace(0,120,1000)';
        hold on; plot(xgrid,pdf(GM_s_mu_mat{i},xgrid),'r-'); hold off
        figure;
        histogram (s_weight_mat(:,i), 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
        xgrid = linspace(0,1.8,1000)';
        hold on; plot(xgrid,pdf(GM_s_weight_mat{i},xgrid),'r-'); hold off
    else
        histogram (s_mu_mat(:,i), 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
        xgrid = linspace(0,1.8,1000)';
        hold on; plot(xgrid,pdf(GM_s_mu_mat{i},xgrid),'r-'); hold off
        figure;
        histogram (s_weight_mat(:, i), 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
        xgrid = linspace(0,1.8,1000)';
        hold on; plot(xgrid,pdf(GM_s_weight_mat{i},xgrid),'r-'); hold off
    end
end


%% First plots of data
%% Plots
% plot force-time measurements of a subject (j) in one plot
figure;
i=110;
plot(time(:,:,i),force(:,:,i))
%plot(time(:,10:15,i),force(:,10:15,i))
xlabel('time');
ylabel('force');

% plot all force-time measurements in one plot
figure;
for i=1:size(force,3) % All subjects
    plot(time(:,:,i),force(:,:,i))
    hold on
end
xlabel('time');
ylabel('force');
hold off

% Scatter diagram of meanF-Dt
figure;
for i=1
plot(X(:,3,i),X(:,2,i),'.')
hold on
end
hold off