%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- mu_weight_statDescription(GMModel, variable)
% Create a statistical description for mu (mean values) and weight (mixing
% proportions of mixture components)
% 
% Approach: Create matrixes with all mu and weights for each subject
% and do a column sort on mu (while also fixing index pairs on weights).
% Then a GMM is fitted for each component-column
%
%%% Returns %%%
%%%
% nd_s_mu_table: components
% Cell array containing the Gaussians describing mu, 1 for each component
%%%
% nd_s_weight_table : components
% Cell array containing the Gaussians describing the mixing proportions of
% mixture components, 1 for each component
%%%
% s_mu_table : subjects x components
% All mu values gathered by the data
%%%
% s_weight_table : subjects x components
% All mixing proportions values gathered by the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nd_s_mu_table, nd_s_weight_table, s_mu_table, s_weight_table] = mu_weight_statDescription(GMModel, variable)

    if variable == 4    % angle
        % Create matrixes
        for i=1:length(GMModel) % all subjects
            mu_table(i,:) = GMModel{i}.mu;
            weight_table (i, :) = GMModel{i}.ComponentProportion;
        end
    else
        % Create matrixes
        for i=1:length(GMModel)
            mu_table(i,:) = GMModel{i}.mu(:,variable);
            weight_table (i, :) = GMModel{i}.ComponentProportion;
        end
    end

    % Sort mu_table and fix index pairs for weight_table
    [s_mu_table s_mu_idx] = sort(mu_table,2);
    % Keep same index for weights
    for i=1:size(weight_table, 1)
        for j=1:size(weight_table, 2)
            s_weight_table(i,j) = weight_table(i,s_mu_idx(i,j));
        end
    end

    % Fit Gaussians to mu and weight matrixes
    for i=1:size(s_mu_table,2)   
        nd_s_mu_table{i} = fitdist(s_mu_table(:,i),'Normal');
        nd_s_weight_table{i} = fitdist(s_weight_table(:,i),'Normal');
    end
    
end