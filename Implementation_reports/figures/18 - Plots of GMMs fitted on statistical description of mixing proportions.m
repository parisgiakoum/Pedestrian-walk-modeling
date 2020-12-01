%% 18 - Plot of GMMs fitted on statistical descripiton of mixing proportion
%% Clear junk, retrieve force-time-position measurements and find meanF-Dt-length-angle for each step of each subject

database=load('steps_database').database_passi;

database = clearDb(database);

[time,force, x_coord, y_coord] = retrieveAllVariables(database);

[X, Dt,meanF, len, angle] = computeAllDesiredVariables(force, time, x_coord, y_coord);

%% GMMs for each subject
[GMModel, h] = fitGMMtoData(X, 5, 'variables');

%% 18
figure;
for i=1:size(s_mu_table,2)  %1-5
    subplot(3,2,i)  
    histogram (s_weight_table(:,i), 'BinWidth', 0.03, 'BinLimits',[0,1.7], 'normalization' , 'pdf' );
    title(strcat('GMM fitted on mixing probability for component', {' '}, num2str(i)));
    xlabel('mixing probability');
    ylabel('Density');
    xgrid = linspace(0,1.8,1000)';
    hold on; plot(xgrid,pdf(GM_s_weight_table{i},xgrid),'r-'); hold off
end