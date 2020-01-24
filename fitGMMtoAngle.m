%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- fitGMMtoAngle(angle, components)
% Fit GMMs to angle
%
% Fit a gaussian mixture model to each subject's angle data and perform
% Two-sample Kolmogorov-Smirnov tests to determine if the distribution
% describes the data
%
%%% Returns %%%
%%%
% GMMAngle: 1x215 (subjects)
% GMM for each subject fitted on the data for angle (1 varible -
% components determined by the variable components)
%%%
% h_angle : 1x215 (subjects)
% Results of Two-sample Kolmogorov-Smirnov tests for angle. If the results
% are 0, then the test is a success, thus the GMM  fitted to data describes
% the data accurately (5% significance level is applied)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% angle is univariate
function [GMMAngle, h_angle] = fitGMMtoAngle(angle, components)

    for i=1:size(angle, 1)   % All subjects

        temp = angle(i,:);

        % GM model fitting on the angle for each subject
        GMMAngle{i} = fitgmdist (temp(~isnan(temp)).', components, 'SharedCovariance',true);

        % Kolmogorov-Smirnov tests
        % Generate random sample based on the GMM
        Y = random(GMMAngle{i},415);

        % Perform the test using the randomly generated data to determine
        % if they belong to the same distribution as the original
        h_angle(i)=kstest2(temp(~isnan(temp)),Y);
    end
end