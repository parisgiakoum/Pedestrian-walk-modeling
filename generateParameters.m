%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- generateParameters(nd_s_weight_table, nd_s_mu_table, nd_Sigma, nd_s_weight_table_angle, nd_s_mu_table_angle, nd_SigmaAngle)
% Extract a random set of parameters to use in simulator's GMM
%
% Takes the Gaussians of each parameter as input and estimates the
% parameters of the final GMM in the output
%
%%% Returns %%%
%%%
% randomWeight : components of final GMM for Dt, meanF, len
% Mixing probabilities Component Proportion Coefficients of each component 
% in the final GMM for Dt, meanF, len
%%%
% randomMu : components of final GMM for Dt, meanF, len x variables (Dt, meanF, len)
% Mean table of the final GMM for Dt, meanF, len
%%%
% randomSigma :  variables (Dt, meanF, len) x variables (Dt, meanF, len)
% Covariance matrix (Sigma) of the final GMM for Dt, meanF, len
%%%
% randomWeightAngle : components of final GMM for angle
% Mixing probabilities (componentProportion) of each component in the final
% GMM for angle
%%%
% randomMuAngle : components of final GMM for angle x variables (angle)
% Mean table of the final GMM for angle
%%%
% randomSigmaAngle :  variable (angle)
% Variance (Sigma) of the final GMM for angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [randomWeight, randomMu, randomSigma, randomWeightAngle, randomMuAngle, randomSigmaAngle] = generateParameters(nd_s_weight_table, nd_s_mu_table, nd_Sigma, nd_s_weight_table_angle, nd_s_mu_table_angle, nd_SigmaAngle)

    % Generate random parameters for Dt, meanF, length
    %%% Random componentProportion, mu
    for i=1:length(nd_s_weight_table) % 1-5
        while 1
            tempWeight = random(nd_s_weight_table{i},1);
            if tempWeight < 0
                tempWeight = random(nd_s_weight_table{i},1);
            else
                randomWeight(i) = tempWeight;
                break;
            end
        end
        randomMuValues(i) = random(nd_s_mu_table{1,i},1);
    end

    % Weight needs to sum to 1, so the values are normalised
    randomWeight = randomWeight./sum(randomWeight);

    %%% Random mu
    randomMu(:, 1) = randomMuValues.';
    for variable=2:3
        
        for i=1:length(nd_s_mu_table) % 1-5
            randomMuValues(i) = random(nd_s_mu_table{variable,i},1);
        end
        randomMu(:, variable) = randomMuValues.';
    end

    %%% Random Sigma
    % Generate random Sigma matrix until it is symmetric positive definite
    while 1
        for i=1:length(nd_Sigma)    % 1-6
            randomSigmaValues(i) = random(nd_Sigma{i},1);
        end

        % 1-3 is the diagonal variances, 4 is 1-2 covariance, 5 is 1-3
        % covariance and 6 is 2-3 covariance
        randomSigma = diag(randomSigmaValues(1:3));
        randomSigma(1,2) = randomSigmaValues(4); randomSigma(2,1) = randomSigmaValues(4);
        randomSigma(1,3) = randomSigmaValues(5); randomSigma(3,1) = randomSigmaValues(5);
        randomSigma(2,3) = randomSigmaValues(6); randomSigma(3,2) = randomSigmaValues(6);

        [~,posdef] = chol(randomSigma); % posdef checks if randomSigma is a symmetric positive definite matrix
        if posdef == 0
            break;
        end
    end

    % Generate random parameters for angle
    %%% Random w, mu
    for i=1:length(nd_s_weight_table_angle) % 1-3
        while 1
            tempWeight = random(nd_s_weight_table_angle{i},1);
            if tempWeight < 0
                tempWeight = random(nd_s_weight_table_angle{i},1);
            else
                randomWeightAngle(i) = tempWeight;
                break;
            end
        end
        randomMuAngle(i) = random(nd_s_mu_table_angle{i},1);
    end
    randomWeightAngle = randomWeightAngle./sum(randomWeightAngle);
    randomMuAngle = randomMuAngle.';

    %%% Random Sigma
    while 1
        randomSigmaAngle = random(nd_SigmaAngle, 1);
        [~,posdef] = chol(randomSigmaAngle); % posdef checks if randomSigma is a symmetric positive definite matrix
        if posdef == 0
            break;
        end
    end
end