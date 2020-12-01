%% 8, 9, 10 - Plotting GMMs fitted in data
%% Clear junk, retrieve force-time-position measurements and find meanF-Dt-length-angle for each step of each subject

database=load('steps_database').database_passi;

database = clearDb(database);

[time,force, x_coord, y_coord] = retrieveAllVariables(database);

[X, Dt,meanF, len, angle] = computeAllDesiredVariables(force, time, x_coord, y_coord);

%% GMMs for each subject
[GMModel, h] = fitGMMtoData(X, 5, 'variables');

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