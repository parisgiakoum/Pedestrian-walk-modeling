%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- fitGMMtoData(X, components)
% Fit GMMs to all data
%
% Fit a gaussian mixture model to each subject's data and perform
% Two-sample Kolmogorov-Smirnov tests to determine if the distributions
% describe the data
%
%%% Returns %%%
%%%
% GMModel: 1x215 (subjects)
% GMM for each subject fitted on the data for Dt-meanF-length (3 varibles -
% components determined by the variable components)
%%%
% h : 3x215 (variables x subjects)
% All results of Two-sample Kolmogorov-Smirnov tests for each variable
% (row). If the result is 0, then the test is a success, thus the GMM
% fitted to data describes the data accurately (5% significance level is
% applied)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [GMModel, h] = fitGMMtoData(X, components)

    for i=1:size(X,3)   % All subjects
        % Use only Dt-meanF-length
        temp = X(:,1:3,i);
   
        % GM model fitting on the 3 variables for each subject
        GMModel{i} = fitgmdist (temp(any(~isnan(temp), 2), :), components, 'SharedCovariance',true);

        % Kolmogorov-Smirnov tests
        % Generate random sample based on the GMM
        Y = random(GMModel{i},415);

        T = X(:,1,i);
        F = X(:,2,i);
        L = X(:,3,i);
        TSim=Y(:,1);             
        FSim=Y(:,2);
        LSim=Y(:,3);

        % Perform the test using the randomly generated data to determine
        % if they belong to the same distribution as the original
        h1(i)=kstest2(T,TSim);
        h2(i)=kstest2(F,FSim);
        h3(i)=kstest2(L,LSim);
    end

    % h will hold the results for all 3 tests
    h=[h1; h2; h3];
    
end