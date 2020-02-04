%% Init
clear;
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

% Generate random Sigma matrix until it is symmetric positive definite
while 1
    for i=1:length(GMModelSigma)    % 1-6
        randomSigmaValues(i) = random(GMModelSigma{i},1);
    end

    % 1-3 is the diagonal, 4 is 1-2, 5 is 1-3 and 6 is 2-3
    randomSigma = diag(randomSigmaValues(1:3));
    randomSigma(1,2) = randomSigmaValues(4); randomSigma(2,1) = randomSigmaValues(4);
    randomSigma(1,3) = randomSigmaValues(5); randomSigma(3,1) = randomSigmaValues(5);
    randomSigma(2,3) = randomSigmaValues(6); randomSigma(3,2) = randomSigmaValues(6);

    [~,posdef] = chol(randomSigma); % posdef checks if randomSigma is a symmetric positive definite matrix
    if posdef == 0
        break;
    end
end

%%% Random w, mu
[GM_s_mu_table, GM_s_weight_table, s_mu_table, s_weight_table] = mu_weight_statDescription(GMModel, 1);

for i=1:length(GM_s_weight_table) % 1-5
    randomWeight(i) = random(GM_s_weight_table{i},1);
    randomMuValues(i) = random(GM_s_mu_table{i},1);
end

% Weight needs to sum to 1, so the values are normalised
randomWeight = randomWeight./sum(randomWeight);

%%% Random mu
randomMu(:, 1) = randomMuValues.';
for variable=2:3
    [GM_s_mu_table, GM_s_weight_table, s_mu_table, s_weight_table] = mu_weight_statDescription(GMModel, variable);
    for i=1:length(GM_s_mu_table) % 1-5
        randomMuValues(i) = random(GM_s_mu_table{i},1);
    end
    randomMu(:, variable) = randomMuValues.';
end

%% Generate random parameters for angle
%%% Random Sigma
[GMModelSigmaAngle, SigmaValuesAngle] = sigmaStatDescription(GMMAngle, 'angle');
randomSigmaAngle = random(GMModelSigmaAngle, 1);

%%% Random w, mu
[GM_s_mu_table_angle, GM_s_weight_table_angle, s_mu_table_angle, s_weight_table_angle] = mu_weight_statDescription(GMMAngle, 4);
for i=1:length(GM_s_weight_table_angle) % 1-3
    randomWeightAngle(i) = random(GM_s_weight_table_angle{i},1);
    randomMuAngle(i) = random(GM_s_mu_table_angle{i},1);
end
randomWeightAngle = randomWeightAngle./sum(randomWeightAngle);
randomMuAngle = randomMuAngle.';

%% Extract a random walk
% Fit a GMM to the generated parameters
simulatedGMM = gmdistribution(randomMu, randomSigma, randomWeight);
simulatedGMMAngle = gmdistribution(randomMuAngle, randomSigmaAngle, randomWeightAngle);

% Simulate a random walk
randomWalk = random(simulatedGMM, 25);
randomWalk(:,4) = random(simulatedGMMAngle, 25);

% Calculate velocity and mean velocity
velocity = randomWalk(:,3)./randomWalk(:,1);
meanVel = mean(velocity).*(3.6); % km/h


%%
%% PLOTS
%% Scatter the random walk positions
% Plot data acquired by simulation
dx = randomWalk(1,3).*cos(randomWalk(1,4));
dy = randomWalk(:,3).*sin(randomWalk(:,4));
for i=2:size(randomWalk , 1)
    dx(i) = dx(i-1)+ randomWalk(i,3).*cos(randomWalk(i,4));
end

figure;
plot(dx,dy,'--*','MarkerEdgeColor','k')
%camroll(90)
for k = 1: length (dx)
    text (dx (k)+0.07, dy (k)+0.01, num2str (k), 'Color','r')
end

%% Scatter walk from database
% Plot data from database to compare
clear xpos ypos
for i=1:3
    if ~isempty(database{2,i})
        xpos(i) = database{2,i}.x_step;
        ypos(i) = database{2,i}.y_step;
    end
end

figure;
plot(xpos,ypos,'-or','MarkerEdgeColor','k')

ylim([-0.2 0.2])
xlim([0 2])
xlabel('x_coord')
ylabel('y\_coord')
title('Calculation of length and angle')

text (xpos (1), ypos (1)+0.04, num2str (1), 'Color','r')
for k = 2: length (xpos)
    line([xpos(k-1),xpos(k-1)], [ypos(k-1),ypos(k)], 'Color', 'm', 'LineWidth', 1, 'LineStyle', '--');
    text (xpos(k-1)+0.01,ypos (k-1) + (ypos(k)-ypos(k-1))/2, strcat('x_',num2str(k-1),'_-_',num2str(k)), 'Color','m')
    line([xpos(k-1),xpos(k)], [ypos(k-1),ypos(k-1)], 'Color', 'b', 'LineWidth', 1,  'LineStyle', '--');
    text (xpos(k-1)+(xpos(k)-xpos(k-1))/2, ypos (k-1)+0.012, strcat('x_',num2str(k-1),'_-_',num2str(k)), 'Color','b')
    
    if ypos(k) < ypos(k-1)
        text (xpos (k), ypos (k)-0.04, num2str (k), 'Color','r')
        text (xpos(k-1)+(xpos(k)-xpos(k-1))/2, ypos (k-1)+ (ypos(k)-ypos(k-1))/2 + 0.015, strcat('length_',num2str(k-1),'_-_',num2str(k)), 'Color','r')
    else
        text (xpos (k), ypos (k)+0.04, num2str (k), 'Color','r')
        text (xpos(k-1)+(xpos(k)-xpos(k-1))/2, ypos (k-1)+ (ypos(k)-ypos(k-1))/2 - 0.015, strcat('length_',num2str(k-1),'_-_',num2str(k)), 'Color','r')
    end
end



%% Plotting GMM and data
%% 8 - Dt
i=31;

tempgm = gmdistribution(GMModel{i}.mu(:,1), GMModel{i}.Sigma(1), GMModel{i}.PComponents);
figure;
temp = X(:,1,i);
histogram (temp(~isnan(temp)), 'BinWidth', 0.02, 'BinLimits',[0,1.6], 'normalization' , 'pdf' );
title('GMM fitted on Dt')
xlabel('Dt Data')
ylabel('Density')

xgrid = linspace(0,1.6,1000)';
hold on; plot(xgrid,pdf(tempgm,xgrid),'r-'); hold off

%% 9 - Mean Force
i=1;

tempgm = gmdistribution(GMModel{i}.mu(:,2), GMModel{i}.Sigma(2,2), GMModel{i}.PComponents);
figure;
temp = X(:,2,i);
histogram (temp(~isnan(temp)), 'BinWidth', 2, 'BinLimits',[0 ,120], 'normalization' , 'pdf' );
title('GMM fitted on mean force')
xlabel('meanF Data')
ylabel('Density')

xgrid = linspace(0,120,1000)';
hold on; plot(xgrid,pdf(tempgm,xgrid),'r-'); hold off

%% 10 - Length
i=20;

tempgm = gmdistribution(GMModel{i}.mu(:,3), GMModel{i}.Sigma(3,3), GMModel{i}.PComponents);
figure;
temp = X(:,3,i);
histogram (temp(~isnan(temp)), 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
title('GMM fitted on length of step')
xlabel('Length Data')
ylabel('Density')

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
% 14-
variable = 1;
[GM_s_mu_table, GM_s_weight_table, s_mu_table, s_weight_table] = mu_weight_statDescription(GMModel, variable);

% Plot gmm and data
for i=1:size(s_mu_table,2)  %1-5
    figure;
    if (variable == 2)
        histogram (s_mu_table(:,i), 'BinWidth', 2, 'BinLimits',[0 ,120], 'normalization' , 'pdf' );
        title(strcat('GMM fitted on mu for component', {' '}, num2str(i),' - mean force variable'));
        xlabel('mu data');
        ylabel('Density');
        xgrid = linspace(0,120,1000)';
        hold on; plot(xgrid,pdf(GM_s_mu_table{i},xgrid),'r-'); hold off
        
        figure;
        histogram (s_weight_table(:,i), 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
        title(strcat('GMM fitted on mixing probability for component', {' '}, num2str(i)));
        xlabel('mixing probability');
        ylabel('Density');
        xgrid = linspace(0,1.8,1000)';
        hold on; plot(xgrid,pdf(GM_s_weight_table{i},xgrid),'r-'); hold off
    else
        histogram (s_mu_table(:,i), 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
        title(strcat('GMM fitted on mu for component', {' '} ,num2str(i),' - Length variable'));
        xlabel('mu data');
        ylabel('Density');
        xgrid = linspace(0,1.8,1000)';
        hold on; plot(xgrid,pdf(GM_s_mu_table{i},xgrid),'r-'); hold off
        
        figure;
        histogram (s_weight_table(:, i), 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
        title(strcat('GMM fitted on mixing probability for component', {' '}, num2str(i)));
        xlabel('mixing probability');
        ylabel('Density');
        xgrid = linspace(0,1.8,1000)';
        hold on; plot(xgrid,pdf(GM_s_weight_table{i},xgrid),'r-'); hold off
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