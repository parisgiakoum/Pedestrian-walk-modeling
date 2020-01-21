%% Init
clear all;
close all;
clc;

%% Retrieve force-time measurements and find meanF-Dt for each step of each subject

database=load('steps_database').database_passi;

database = clearDb(database);

[time,force, x_coord, y_coord] = retrieveAllVariables(database);

[X, Dt,meanF, len, phi] = computeAllDesiredVariables(force, time, x_coord, y_coord);

%% GM model
for i=1:size(X,3)   %1-215
    temp = X(:,1:3,i);

    GMModel{i} = fitgmdist (temp(any(~isnan(temp), 2), :), 5, 'SharedCovariance',true);
  
    Y = random(GMModel{i},415);

    T = X(:,1,i);
    F = X(:,2,i);
    L = X(:,3,i);
    TSim=Y(:,1);             
    FSim=Y(:,2);
    LSim=Y(:,3);

    h1(i)=kstest2(T,TSim);
    h2(i)=kstest2(F,FSim);
    h3(i)=kstest2(L,LSim);
end

h=[h1; h2; h3];
% sum(h==1,'all')
clear h1 h2 h3 T F L TSim FSim LSim

%% Try plotting gmm and data
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