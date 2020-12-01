%% Clear junk, retrieve force-time-position measurements and find meanF-Dt-length-angle for each step of each subject

database=load('steps_database').database_passi; % load database

database = clearDb(database);   % Clear database

[time,force, x_coord, y_coord] = retrieveAllVariables(database);    %  Retrieve forces, times and coordinates

[X, Dt, meanF, len, angle] = computeAllDesiredVariables(force, time, x_coord, y_coord);  % Extract Dt, meanf, len, angle

%%

Xtotal = X(:,:,1);
for i=2:size(X, 3)
    Xtotal = [Xtotal ; X(:,:,i)];
end

Xtotal2 = Xtotal(any(~isnan(Xtotal), 2), :);

%normal
counts(2,:) = zeros(1,15);
for itter=1:5
    for i=1:size(X,3)   % All subjects
        for comp=1:15
            temp = X(:,1:3,i);
            % 0.1
            GMM = fitgmdist (temp(any(~isnan(temp), 2), :), comp, 'RegularizationValue', 0.01, 'Options', statset('MaxIter', 1500), 'SharedCovariance',true);
            BIC(i,comp) = GMM.BIC;
        end
    end
    [mint, indexmin] = min(BIC,[],2);
    for i=1:15  % How many minBics in each component number
        counts(2,i) = counts(2,i) + sum(indexmin==i);
    end
end

figure;
for subj = 1:5
    nexttile
    plot(BIC(subj, :), '.-', 'DisplayName','BIC')
    lgd = legend;
end


%%
% First approach: regularization of GMMs and spot mean BIC and AIC to find
% the best fitting
for j=1:10 % j tries
    for i=1:size(X,3)   % All subjects
        for comp=1:10
            temp = X(:,1:3,i);
            GMM = fitgmdist (temp(any(~isnan(temp), 2), :), comp, 'RegularizationValue', 0.01, 'Options', statset('MaxIter', 1500), 'SharedCovariance',true);
            BIC(i,comp) = GMM.BIC;
            AIC(i,comp) = GMM.AIC;
        end 
    end
    meanBIC(j,:) = mean(BIC);
    meanAIC(j,:) = mean(AIC);
end

[val,idx] = min(meanBIC,[],2) % MINBIC IS IN 4 COMPONENTS
[val,idx] = min(meanAIC,[],2) % MINAIC IS IN 5 COMPONENTS