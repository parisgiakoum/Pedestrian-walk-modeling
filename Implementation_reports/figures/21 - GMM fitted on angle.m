%% 21 - Plotting GMMs fitted on angle
%% Clear junk, retrieve force-time-position measurements and find meanF-Dt-length-angle for each step of each subject

database=load('steps_database').database_passi;

database = clearDb(database);

[time,force, x_coord, y_coord] = retrieveAllVariables(database);

[X, Dt,meanF, len, angle] = computeAllDesiredVariables(force, time, x_coord, y_coord);

%% GMMs for each subject fitted on angle
[GMMAngle, h_angle] = fitGMMtoData(X, 3, 'angle');


%% 21 - angle
i=20;

figure;
temp = X(:,4,i);
histogram (temp(~isnan(temp)), 'BinWidth', 0.03, 'BinLimits',[-0.25,0.25], 'normalization' , 'pdf' );
title('GMM fitted on angle of step')
xlabel('Angle Data')
ylabel('Density')
ylim([0 11])

xgrid = linspace(-0.25,0.25,1000)';
hold on; plot(xgrid,pdf(GMMAngle{i},xgrid),'r-'); hold off