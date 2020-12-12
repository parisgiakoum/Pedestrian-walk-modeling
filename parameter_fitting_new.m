
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

%% UNSORTED
% Modelisation of mu 
%% GAUSSIAN METHOD
for variable = 1:4
    if variable == 4    % angle
        % Create matrixes
        for i=1:length(GMModel) %1-215
            mu_table_angle(i,:) = GMMAngle{i}.mu;
            weight_table_angle(i, :) = GMMAngle{i}.ComponentProportion;
        end
    else
        % Create matrixes
        for i=1:length(GMModel) %1-215
            mu_table(i,:) = GMModel{i}.mu(:,variable);
            weight_table (i, :) = GMModel{i}.ComponentProportion;
        end
    end

    % PLOT    
    if variable == 4    % angle
        figure;
        for i=1:size(mu_table_angle,2)  % all components
            subplot(2,1,i)
            pd = fitdist(mu_table_angle(:,i),'Normal')
            histogram (mu_table_angle(:,i), 'BinWidth', 0.03, 'BinLimits',[-1.4,1.6], 'normalization' , 'pdf' );
            title(strcat('Gaussian fitted on unsorted mu for component ', {' '} ,num2str(i),' - angle variable'));
            xlabel('mu data');
            ylabel('Density');
            xgrid = linspace(-1.4,1.6,1000)';
            hold on; plot(xgrid,pdf(pd,xgrid),'r-'); hold off
        end
        
        figure;
        for i=1:size(weight_table_angle,2)  % all components
            subplot(2,1,i)
            pd = fitdist(weight_table_angle(:,i),'Normal')
            histogram (weight_table_angle(:,i), 'BinWidth', 0.03, 'BinLimits',[0,1], 'normalization' , 'pdf' );
            title(strcat('GMM fitted on unsorted mixing probability for component ', {' '} ,num2str(i),' - angle variable'));
            xlabel('mu data');
            ylabel('Density');
            xgrid = linspace(0,1,1000)';
            hold on; plot(xgrid,pdf(pd,xgrid),'r-'); hold off
        end
else
        figure;
        for i=1:size(mu_table,2)  % all components
            if (variable == 2)  % meanF
                subplot(3,1,i)
                pd = fitdist(mu_table(:,i),'Normal')
                histogram (mu_table(:,i), 'BinWidth', 2, 'BinLimits',[0 ,120], 'normalization' , 'pdf' );
                title(strcat('Gaussian fitted on unsorted mu for component ', {' '}, num2str(i),' - mean force variable'));
                xlabel('mu data');
                ylabel('Density');
                % Curve of GMM
                xgrid = linspace(0,120,1000)';
                hold on; plot(xgrid,pdf(pd,xgrid),'r-'); hold off
            elseif  (variable == 3) % len
                subplot(3,1,i)
                pd = fitdist(mu_table(:,i),'Normal')
                histogram (mu_table(:,i), 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
                title(strcat('Gaussian fitted on unsorted mu for component ', {' '} ,num2str(i),' - Length variable'));
                xlabel('mu data');
                ylabel('Density');
                xgrid = linspace(0,1.8,1000)';
                hold on; plot(xgrid,pdf(pd,xgrid),'r-'); hold off
            else    % Dt
                subplot(3,1,i)
                pd = fitdist(mu_table(:,i),'Normal')
                histogram (mu_table(:,i), 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
                title(strcat('Gaussian fitted on unsorted mu for component ', {' '} ,num2str(i),' - Dt variable'));
                xlabel('mu data');
                ylabel('Density');
                xgrid = linspace(0,1.8,1000)';
                hold on; plot(xgrid,pdf(pd,xgrid),'r-'); hold off 
            end
        end
        
        figure;
        for i=1:size(weight_table,2)  % all components
            subplot(3,1,i)  
            pdw = fitdist(weight_table(:,i),'Normal')
            histogram (weight_table(:,i), 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
            title(strcat('GMM fitted on unsorted mixing probability for component', {' '}, num2str(i), ' - Dt, F, L'));
            xlabel('mixing probability');
            ylabel('Density');
            xgrid = linspace(0,1.8,1000)';
            hold on; plot(xgrid,pdf(pdw,xgrid),'r-'); hold off
        end
    end
end
%% SORTED
% Modelisation of mu 
%% GAUSSIAN METHOD
variable = 4
if (variable == 4)    % angle
    % Create matrixes
    for i=1:length(GMModel) %1-215
        mu_table(i,:) = GMMAngle{i}.mu;
        weight_table (i, :) = GMMAngle{i}.ComponentProportion;
    end
    [s_mu_table_angle s_mu_idx] = sort(mu_table,2);
    for i=1:size(weight_table, 1) %1-215
        for j=1:size(weight_table, 2)
            s_weight_table_angle(i,j) = weight_table(i,s_mu_idx(i,j));
        end
    end
elseif variable == 1
    % Create matrixes
    for i=1:length(GMModel) %1-215
        mu_table(i,:) = GMModel{i}.mu(:,variable);
        weight_table (i, :) = GMModel{i}.ComponentProportion;
    end
    [s_mu_table s_mu_idx] = sort(mu_table,2);
    for i=1:size(weight_table, 1) %1-215
        for j=1:size(weight_table, 2)
            s_weight_table(i,j) = weight_table(i,s_mu_idx(i,j));
        end
    end
else
    for i=1:length(GMModel) %1-215
        mu_table(i,:) = GMModel{i}.mu(:,variable);
        weight_table (i, :) = GMModel{i}.ComponentProportion;
    end
end
    % Sort mu_table and fix index pairs for weight_table
%    [s_mu_table s_mu_idx] = sort(mu_table,2);
    % Keep same index for weights
    % Keep same index for weights
%     for i=1:size(weight_table, 1) %1-215
%          for j=1:size(weight_table, 2)
%              s_weight_table(i,j) = weight_table(i,s_mu_idx(i,j));
%          end
%      end


% PLOT
if variable==4
    figure;
    for i=1:size(s_mu_table_angle,2)  % all components
        subplot(2,1,i)
        pd = fitdist(s_mu_table_angle(:,i),'Normal')
        histogram (s_mu_table_angle(:,i), 'BinWidth', 0.03, 'BinLimits',[-1.4,1.6], 'normalization' , 'pdf' );
        title(strcat('Gaussian fitted on sorted mu for component ', {' '} ,num2str(i),' - angle variable'));
        xlabel('mu data');
        ylabel('Density');
        xgrid = linspace(-1.4,1.6,1000)';
        hold on; plot(xgrid,pdf(pd,xgrid),'r-'); hold off
    end
    figure;
    for i=1:size(s_weight_table_angle,2)  % all components
        subplot(2,1,i)
        pd = fitdist(s_weight_table_angle(:,i),'Normal')
        histogram (s_weight_table_angle(:,i), 'BinWidth', 0.03, 'BinLimits',[0,1], 'normalization' , 'pdf' );
        title(strcat('GMM fitted on mixing probability for component ', {' '} ,num2str(i),' - angle variable'));
        xlabel('mu data');
        ylabel('Density');
        xgrid = linspace(0,1,1000)';
        hold on; plot(xgrid,pdf(pd,xgrid),'r-'); hold off
    end
else
    figure;
    for i=1:size(s_mu_table,2)  % all components
        subplot(3,1,i)
        pd = fitdist(s_mu_table(:,i),'Normal')
        if (variable == 2)  % meanF
            histogram (s_mu_table(:,i), 'BinWidth', 2, 'BinLimits',[0 ,120], 'normalization' , 'pdf' );
            title(strcat('Gaussian fitted on sorted mu for component ', {' '}, num2str(i),' - mean force variable'));
            xlabel('mu data');
            ylabel('Density');
            % Curve of GMM
            xgrid = linspace(0,120,1000)';
            hold on; plot(xgrid,pdf(pd,xgrid),'r-'); hold off
        elseif  (variable == 3) % len
            pd = fitdist(s_mu_table(:,i),'Normal')
            histogram (s_mu_table(:,i), 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
            title(strcat('Gaussian fitted on sorted mu for component ', {' '} ,num2str(i),' - Length variable'));
            xlabel('mu data');
            ylabel('Density');
            xgrid = linspace(0,1.8,1000)';
            hold on; plot(xgrid,pdf(pd,xgrid),'r-'); hold off
        else    % Dt
            pd = fitdist(s_mu_table(:,i),'Normal')
            histogram (s_mu_table(:,i), 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
            title(strcat('Gaussian fitted on sorted mu for component ', {' '} ,num2str(i),' - Dt variable'));
            xlabel('mu data');
            ylabel('Density');
            xgrid = linspace(0,1.8,1000)';
            hold on; plot(xgrid,pdf(pd,xgrid),'r-'); hold off
        end

    end

    figure;
    for i=1:size(s_weight_table,2)  % all components
        subplot(3,1,i)  
        pdw = fitdist(s_weight_table(:,i),'Normal')
        histogram (s_weight_table(:,i), 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
        title(strcat('GMM fitted on mixing probability for component', {' '}, num2str(i)));
        xlabel('mixing probability');
        ylabel('Density');
        xgrid = linspace(0,1.8,1000)';
        hold on; plot(xgrid,pdf(pdw,xgrid),'r-'); hold off
    end
end