%% 23, 24 - Plots of GMMs fitted on statistical description of mu and proportions for angles
%% Clear junk, retrieve force-time-position measurements and find meanF-Dt-length-angle for each step of each subject

database=load('steps_database').database_passi;

database = clearDb(database);

[time,force, x_coord, y_coord] = retrieveAllVariables(database);

[X, Dt,meanF, len, angle] = computeAllDesiredVariables(force, time, x_coord, y_coord);

%% GMMs for each subject
[GMMAngle, h_angle] = fitGMMtoData(X, 3, 'angle');

%% GMMs for mu and component proportion
[GM_s_mu_table_angle, GM_s_weight_table_angle, s_mu_table_angle, s_weight_table_angle] = mu_weight_statDescription(GMMAngle, 4);

%% 23 - mu of angle
figure;
for i=1:size(s_mu_table_angle,2)  %1-3
    subplot(2,2,i)
    histogram (s_mu_table_angle(:,i), 'BinWidth', 0.03, 'BinLimits',[-1.4,1.6], 'normalization' , 'pdf' );
    title(strcat('GMM fitted on mu for component', {' '} ,num2str(i),' - angle variable'));
    xlabel('mu data');
    ylabel('Density');
    xgrid = linspace(-1.4,1.6,1000)';
    hold on; plot(xgrid,pdf(GM_s_mu_table_angle{i},xgrid),'r-'); hold off
end

%% 24 - component proportion of angle
figure;
for i=1:size(s_weight_table_angle,2)  %1-3
    subplot(2,2,i)
    histogram (s_weight_table_angle(:,i), 'BinWidth', 0.03, 'BinLimits',[0,1], 'normalization' , 'pdf' );
    title(strcat('GMM fitted on mu for component', {' '} ,num2str(i),' - angle variable'));
    xlabel('mu data');
    ylabel('Density');
    xgrid = linspace(0,1,1000)';
    hold on; plot(xgrid,pdf(GM_s_weight_table_angle{i},xgrid),'r-'); hold off
end
