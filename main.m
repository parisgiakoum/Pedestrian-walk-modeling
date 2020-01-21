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