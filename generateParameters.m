%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- generateParameters(GM_s_weight_table, GM_s_mu_table, GMModelSigma, GM_s_weight_table_angle, GM_s_mu_table_angle, GMModelSigmaAngle)
% Extract a random set of parameters to use in the GMM of the simulator
%
% Takes the GMMs of each parameter as input (for angle and all the other
% variables seperately) and gives the parameters of the final GMM as output
%
%%% Returns %%%
%%%
% randomWeight : 1x5 (components of final GMM for Dt, meanF, len)
% Mixing probabilities (componentProportion) of each component in the final
% GMM for Dt, meanF, len
%%%
% randomMu : 5x3 (components of final GMM for Dt, meanF, len x variables (Dt, meanF, len))
% Mean table of the final GMM for Dt, meanF, len
%%%
% randomSigma :  3x3 (variables (Dt, meanF, len) x variables (Dt, meanF, len))
% Covariance matrix (Sigma) of the final GMM for Dt, meanF, len
%%%
% randomWeightAngle : 1x3 (components of final GMM for angle)
% Mixing probabilities (componentProportion) of each component in the final
% GMM for angle
%%%
% randomMuAngle : 3x1 (components of final GMM for angle x variables (angle))
% Mean table of the final GMM for angle
%%%
% randomSigmaAngle :  val (variable (angle))
% Variance (Sigma) of the final GMM for angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [randomWeight, randomMu, randomSigma, randomWeightAngle, randomMuAngle, randomSigmaAngle] = generateParameters(GM_s_weight_table, GM_s_mu_table, GMModelSigma, GM_s_weight_table_angle, GM_s_mu_table_angle, GMModelSigmaAngle)

    % Generate random parameters for Dt, meanF, length
    %%% Random componentProportion, mu
    for i=1:length(GM_s_weight_table) % 1-5
        while 1
            tempWeight = random(GM_s_weight_table{i},1);
            if tempWeight < 0
                tempWeight = random(GM_s_weight_table{i},1);
            else
                randomWeight(i) = tempWeight;
                break;
            end
        end
        randomMuValues(i) = random(GM_s_mu_table{1,i},1);
    end

    % Weight needs to sum to 1, so the values are normalised
    randomWeight = randomWeight./sum(randomWeight);

    %%% Random mu
    randomMu(:, 1) = randomMuValues.';
    for variable=2:3
        
        for i=1:length(GM_s_mu_table) % 1-5
            randomMuValues(i) = random(GM_s_mu_table{variable,i},1);
        end
        randomMu(:, variable) = randomMuValues.';
    end

    %%% Random Sigma
    % Generate random Sigma matrix until it is symmetric positive definite
    while 1
        for i=1:length(GMModelSigma)    % 1-6
            randomSigmaValues(i) = random(GMModelSigma{i},1);
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
    for i=1:length(GM_s_weight_table_angle) % 1-3
        while 1
            tempWeight = random(GM_s_weight_table_angle{i},1);
            if tempWeight < 0
                tempWeight = random(GM_s_weight_table_angle{i},1);
            else
                randomWeightAngle(i) = tempWeight;
                break;
            end
        end
        randomMuAngle(i) = random(GM_s_mu_table_angle{i},1);
    end
    randomWeightAngle = randomWeightAngle./sum(randomWeightAngle);
    randomMuAngle = randomMuAngle.';

    %%% Random Sigma
    while 1
        randomSigmaAngle = random(GMModelSigmaAngle, 1);
        [~,posdef] = chol(randomSigmaAngle); % posdef checks if randomSigma is a symmetric positive definite matrix
        if posdef == 0
            break;
        end
    end
end