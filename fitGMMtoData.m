 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- fitGMMtoData(X, components, mode)
% Fit GMMs to data
%
% Fit a gaussian mixture model to each subject's data
%
%%% Returns %%%
%%%
% GMModel: subjects
% GMM for each subject fitted on the data for Dt-meanF-length (3 varibles -
% components determined by the variable components)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GMModel = fitGMMtoData(X, components, mode)

    if strcmp(mode, 'variables')
        for i=1:size(X,3)   % All subjects
            % Use only Dt-meanF-length
            temp = X(:,1:3,i);

            % GM model fitting on the 3 variables for each subject
            GMModel{i} = fitgmdist (temp(any(~isnan(temp), 2), :), components, 'Options', statset('MaxIter', 1500), 'SharedCovariance',true);

        end

    elseif strcmp(mode, 'angle')

        for i=1:size(X,3)   % All subjects
            temp = X(:,4,i);

            % GM model fitting on the angle for each subject
            GMModel{i} = fitgmdist (temp(~isnan(temp)), components, 'Options', statset('MaxIter', 1500), 'SharedCovariance',true);

        end
        
    else
        
        error('Wrong mode input! Put either variables or angle');

    end
    
end