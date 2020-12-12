function [GM_mu_table, GM_weight_table, mu_table, weight_table] = mu_weight_statDescription_no_sorting(GMModel, variable)

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
  %  [s_mu_table s_mu_idx] = sort(mu_table,2);
    % Keep same index for weights
%    for i=1:size(weight_table, 1) %1-215
%         for j=1:size(weight_table, 2)
%             s_weight_table(i,j) = weight_table(i,s_mu_idx(i,j));
%         end
%     end
        GM_mu_table_test = fitgmdist (mu_table , 1, 'Options', statset('MaxIter', 1500), 'SharedCovariance',true);

    for i=1:size(mu_table,2)   
       [m(i),s(i)] = normfit(mu_table(:,i));
    end
       
    % Fit GMM to mu and weight matrixes
    for i=1:size(mu_table,2)   
        GM_mu_table{i} = fitgmdist (mu_table(:,i) , 3, 'Options', statset('MaxIter', 1500), 'SharedCovariance',true);
        GM_weight_table{i} = fitgmdist (weight_table(:,i) , 3, 'Options', statset('MaxIter', 1500), 'SharedCovariance',true);
    end
    
end