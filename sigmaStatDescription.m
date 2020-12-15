%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- sigmaStatDescription(GMModel, mode)
% Create a statistical description for Sigma (covariances)
%
% The covariances matrix is a 3x3 (as 3 are the variables) matrix with 6 
% significant values (the matrix is symetrical so the rest values are equal
% ie. (1,2)=(2,1)). In this function, a GMM is fitted on each of the 6
% values of significance (6 GMMs)
%
%%% Returns %%%
%%%
% GMModelSigma: 1x6 (significant values)
% Cell array containing 6 GMMs, 1 for each significant value (1 variable -
% 4 components each)
%%%
% SigmaValues : 215x6 (subjects x significant values)
% All significant values gathered by the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nd_Sigma, SigmaValues] = sigmaStatDescription(GMModel, mode)

    if strcmp(mode, 'variables')
        
        % 1-3 is the diagonal, 4 is 1-2, 5 is 1-3 and 6 is 2-3
        SigmaValues = [diag(GMModel{1}.Sigma).'  GMModel{1}.Sigma(1,2)  GMModel{1}.Sigma(1,3)  GMModel{1}.Sigma(2,3)];
        for i=2:length(GMModel) % 1-215
           SigmaValues = [SigmaValues ; diag(GMModel{i}.Sigma).' GMModel{i}.Sigma(1,2)  GMModel{i}.Sigma(1,3)  GMModel{i}.Sigma(2,3)];
        end

        for i=1:size(SigmaValues, 2) % 1-6
            nd_Sigma{i} = fitdist(SigmaValues(:,i),'Normal');
        end
        
    elseif strcmp(mode, 'angle')
        
        for i=1:length(GMModel) % 1-215
            SigmaValues(i) = GMModel{i}.Sigma;
        end
        SigmaValues = SigmaValues.';
        
        nd_Sigma = fitdist(SigmaValues,'Normal');
        
    else
        
        error('Wrong input! Put either variables or angle');
        
    end
    
end