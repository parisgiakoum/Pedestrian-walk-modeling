%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- sigmaStatDescription(GMModel, mode)
% Create a statistical description for Sigma (covariances/variances)
%
% The covariances matrix is symmetriv matrix with critical values e.g.,
% the Sigma of a 3-component GMM is a 3x3 symmetric matrix with 6 critical
% values. In this function, a GMM is fitted on each critical value
%
%%% Returns %%%
%%%
% GMModelSigma: critical values
% Cell array containing Gaussians, 1 for each significant value
%%%
% SigmaValues : subjects x critical values
% All critical values contained in the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nd_Sigma, SigmaValues] = sigmaStatDescription(GMModel, mode)

    if strcmp(mode, 'variables')
        
        % 1-3 is the diagonal, 4 is 1-2, 5 is 1-3 and 6 is 2-3
        SigmaValues = [diag(GMModel{1}.Sigma).'  GMModel{1}.Sigma(1,2)  GMModel{1}.Sigma(1,3)  GMModel{1}.Sigma(2,3)];
        for i=2:length(GMModel) % 1-215
           SigmaValues = [SigmaValues ; diag(GMModel{i}.Sigma).' GMModel{i}.Sigma(1,2)  GMModel{i}.Sigma(1,3)  GMModel{i}.Sigma(2,3)];
        end

        for i=1:size(SigmaValues, 2) % all critical values
            nd_Sigma{i} = fitdist(SigmaValues(:,i),'Normal');
        end
        
    elseif strcmp(mode, 'angle')
        
        for i=1:length(GMModel) % all subjects
            SigmaValues(i) = GMModel{i}.Sigma;
        end
        SigmaValues = SigmaValues.';
        
        nd_Sigma = fitdist(SigmaValues,'Normal');
        
    else
        
        error('Wrong input! Put either variables or angle');
        
    end
    
end