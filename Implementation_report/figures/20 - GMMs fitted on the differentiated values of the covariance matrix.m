%% 20 - Plot of GMMs fitted on the differentiated values of the covariance matrix
%% Clear junk, retrieve force-time-position measurements and find meanF-Dt-length-angle for each step of each subject

database=load('steps_database').database_passi;

database = clearDb(database);

[time,force, x_coord, y_coord] = retrieveAllVariables(database);

[X, Dt,meanF, len, angle] = computeAllDesiredVariables(force, time, x_coord, y_coord);

%% GMMs for each subject
[GMModel, h] = fitGMMtoData(X, 5, 'variables');

%% 20
% 1-3 is the diagonal, 4 is 1-2, 5 is 1-3 and 6 is 2-3
[GMModelSigma, SigmaValues] = sigmaStatDescription(GMModel, 'variables');

figure;
for i=1:size(SigmaValues, 2) % 1-6

    subplot(3,2,i)
    if i==1 | i==3 | i==5
        histogram (SigmaValues(:,i), 'normalization' , 'pdf'  );
        title(strcat('GMM fitted on differentiated value', {' '}, num2str(i)))
        xlabel('Data');
        ylabel('Density');
        xgrid = linspace(-0.1,0.08,1000)';
        hold on; plot(xgrid,pdf(GMModelSigma{i},xgrid),'r-'); hold off
    elseif i==2
        histogram (SigmaValues(:,i), 'normalization' , 'pdf'  );
        title(strcat('GMM fitted on differentiated value', {' '}, num2str(i)))
        xlabel('Data');
        ylabel('Density');
        xgrid = linspace(-0.1,20,1000)';
        hold on; plot(xgrid,pdf(GMModelSigma{i},xgrid),'r-'); hold off
    else
        histogram (SigmaValues(:,i), 'normalization' , 'pdf'  );
        title(strcat('GMM fitted on differentiated value', {' '}, num2str(i)))
        xlabel('Data');
        ylabel('Density');
        xgrid = linspace(-0.35,0.5,1000)';
        hold on; plot(xgrid,pdf(GMModelSigma{i},xgrid),'r-'); hold off
    end
    
end