%% Init
clear all;
close all;
clc;

%% Gain Force-time measurements and find meanF-Dt for each step of a subject
database=load('steps_database.mat').database_passi;
subj = 2;

% identify max force-time indices
max_length1 = 0;
for i=1:size(database,2) %1-32
    if ~isempty(database{subj,i})
        if length(database{subj,i}.force) > max_length1
            max_length1=length(database{subj,i}.force);
        end
    end
end

% fill with nan all forces and times up to max indices
for i=1:size(database,2) %1-32
    if ~isempty(database{subj,i})
        database{subj,i}.force(end+1:max_length1)=nan;
        database{subj,i}.time(end+1:max_length1)=nan;
    end
end

% Find forces-time for each step
for i=1:size(database,2) %% 1-max steps
    if ~isempty(database{subj,i})
        force(:,i)=database{subj,i}.force;
        time(:,i)=database{subj,i}.time;
    else
        force(:,i)=nan;
        time(:,i)=nan;
    end
end

% plot all force-time measurements in one plot for subject 1
figure;
plot(time,force)

% Find mean force of each step
for i=1:size(force,2) %1-max steps
    meanF(:,i) = nanmean(nonzeros(force(:,i)));
    Dt(:,i) = max(time(:,i))-time(1,i);
end

% Scatter diagram of meanF-Dt for all steps
figure;
scatter(Dt,meanF,'.')
