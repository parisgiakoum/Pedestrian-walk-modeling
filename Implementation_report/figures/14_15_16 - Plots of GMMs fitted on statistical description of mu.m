%% 14-15-16 - Plot of GMMs fitted on statistical descripiton of mu
%% Clear junk, retrieve force-time-position measurements and find meanF-Dt-length-angle for each step of each subject

database=load('steps_database').database_passi;

database = clearDb(database);

[time,force, x_coord, y_coord] = retrieveAllVariables(database);

[X, Dt,meanF, len, angle] = computeAllDesiredVariables(force, time, x_coord, y_coord);

%% GMMs for each subject
[GMModel, h] = fitGMMtoData(X, 5, 'variables');

%% 14-15-16
for variable = 1:3
    [GM_s_mu_table, GM_s_weight_table, s_mu_table, s_weight_table] = mu_weight_statDescription(GMModel, variable);
    figure;
    for i=1:size(s_mu_table,2)  %1-5
        subplot(3,2,i)
        if (variable == 2)
            histogram (s_mu_table(:,i), 'BinWidth', 2, 'BinLimits',[0 ,120], 'normalization' , 'pdf' );
            title(strcat('GMM fitted on mu for component', {' '}, num2str(i),' - mean force variable'));
            xlabel('mu data');
            ylabel('Density');
            xgrid = linspace(0,120,1000)';
            hold on; plot(xgrid,pdf(GM_s_mu_table{i},xgrid),'r-'); hold off
        elseif  (variable == 3)
            histogram (s_mu_table(:,i), 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
            title(strcat('GMM fitted on mu for component', {' '} ,num2str(i),' - Length variable'));
            xlabel('mu data');
            ylabel('Density');
            xgrid = linspace(0,1.8,1000)';
            hold on; plot(xgrid,pdf(GM_s_mu_table{i},xgrid),'r-'); hold off
        else
            histogram (s_mu_table(:,i), 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
            title(strcat('GMM fitted on mu for component', {' '} ,num2str(i),' - Dt variable'));
            xlabel('mu data');
            ylabel('Density');
            xgrid = linspace(0,1.8,1000)';
            hold on; plot(xgrid,pdf(GM_s_mu_table{i},xgrid),'r-'); hold off
        end
    end
end

