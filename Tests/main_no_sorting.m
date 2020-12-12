clear;
close all;
clc;

%% Basic functions
%% Clear junk, retrieve force-time-position measurements and find meanF-Dt-length-angle for each step of each subject

database=load('steps_database').database_passi; % load database

database = clearDb(database);   % Clear database

[time,force, x_coord, y_coord] = retrieveAllVariables(database);    %  Retrieve forces, times and coordinates

[X, Dt, meanF, len, angle] = computeAllDesiredVariables(force, time, x_coord, y_coord);  % Extract Dt, meanf, len, angle


%% Modelisation of data
% Fit GMMs on each person
[GMModel, h] = fitGMMtoData(X, 3, 'variables'); % Dt, meanF, len
[GMMAngle, h_angle] = fitGMMtoData(X, 2, 'angle');  % angle

% Statistical description of parameters
% Statistical description of mu and componentProportion
[GM_mu_table, GM_weight_table, mu_table{1}, weight_table] = mu_weight_statDescription_no_sorting(GMModel, 1); % mu for Dt and comonentProportion
[GM_mu_table_angle, GM_weight_table_angle, mu_table_angle, weight_table_angle] = mu_weight_statDescription_no_sorting(GMMAngle, 4);    % mu and componentProportion for angle
% mu for meanF, len
for variable=2:3
    [GM_mu_table(variable, :), ~, mu_table{variable}, ~] = mu_weight_statDescription_no_sorting(GMModel, variable);
end

% Statistical description of Sigma
[GMModelSigma, SigmaValues] = sigmaStatDescription(GMModel, 'variables');   % Dt, meanF, len
[GMModelSigmaAngle, SigmaValuesAngle] = sigmaStatDescription(GMMAngle, 'angle');    % angle


%%
for variable = 1:3
    figure;
    for i=1:size(mu_table,2)  % all components
            subplot(3,1,i)
            tempgm = gmdistribution(GM_mu_table_test.mu(i), GM_mu_table_test.Sigma(i,i), GM_mu_table_test.PComponents); % fit a GMM with parameters corresponding to meanF
            histogram (mu_table(:,i), 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
            title(strcat('GMM fitted on mu for component ', {' '} ,num2str(i),' - Dt variable'));
            xlabel('mu data');
            ylabel('Density');
            xgrid = linspace(0,1.8,1000)';
            hold on; plot(xgrid,pdf(tempgm,xgrid),'r-'); hold off
    end
    end

y = normpdf(x,m,s);
pd = fitdist(mu_table(:,i),'Normal')

    figure;
    for i=1:size(mu_table,2)  % all components
             subplot(3,1,i)
          %  y = normpdf(mu_table(:,i),m(i),s(i));
            pd = fitdist(mu_table(:,i),'Normal')
            histogram (mu_table(:,i), 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
            title(strcat('GMM fitted on mu for component ', {' '} ,num2str(i),' - Dt variable'));
            xlabel('mu data');
            ylabel('Density');
            xgrid = linspace(0,1.8,1000)';
            hold on; plot(xgrid,pdf(pd,xgrid),'r-'); hold off
    end
    
    normrnd(pd.mu,pd.sigma)
