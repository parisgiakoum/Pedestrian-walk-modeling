%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- mu_weight_statDescription(GMModel, variable)
% Create a statistical description for mu (mean values) and weight (mixing
% proportions of mixture components)
% 
% Approach: Create matrixes with all mu and weights for each subject (215x5
% as they contain all values for each subject and each component for a
% specific variable) and do a column sort on mu (while also fixing index
% pairs on weights). Then a GMM is fitted for each of the 5
% components-columns (10 in total - 5 for mu and 5 for weight)
%
%%% Returns %%%
%%%
% GM_s_mu_table: 1x5 (components)
% Cell array containing 5 GMMs to describe mu, 1 for each component (1
% variable - 3 components each)
%%%
% GM_s_weight_table : 1x5 (components)
% Cell array containing 5 GMMs to describe the mixing proportions of
% mixture components, 1 for each component (1 variable - 3 components each)
%%%
% s_mu_table : 215x5 (subjects x components)
% All mu values gathered by the data
%%%
% s_weight_table : 215x5 (subjects x components)
% All mixing proportions values gathered by the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nd_s_mu_table, nd_s_weight_table, s_mu_table, s_weight_table] = mu_weight_statDescription(GMModel, variable)

    if variable == 4    % angle
        % Create matrixes
        for i=1:length(GMModel) %1-215
            mu_table(i,:) = GMModel{i}.mu;
            weight_table (i, :) = GMModel{i}.ComponentProportion;
        end
    else
        % Create matrixes
        for i=1:length(GMModel) %1-215
            mu_table(i,:) = GMModel{i}.mu(:,variable);
            weight_table (i, :) = GMModel{i}.ComponentProportion;
        end
    end

    % Sort mu_table and fix index pairs for weight_table
    [s_mu_table s_mu_idx] = sort(mu_table,2);
    % Keep same index for weights
    for i=1:size(weight_table, 1) %1-215
        for j=1:size(weight_table, 2)
            s_weight_table(i,j) = weight_table(i,s_mu_idx(i,j));
        end
    end

    % Fit GMM to mu and weight matrixes
    for i=1:size(s_mu_table,2)   
        nd_s_mu_table{i} = fitdist(s_mu_table(:,i),'Normal');
        nd_s_weight_table{i} = fitdist(s_weight_table(:,i),'Normal')
    end
    
end