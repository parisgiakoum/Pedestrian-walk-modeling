%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- computeAllDesiredVariables(force, time, x_coord, y_coord)
% Compute interarrival time between one step and the next, mean vertical force
% induced, length and angle for each step - Clear wrong values
%
%%% Returns %%%
%%%
% X: steps x Dt-MeanF-length-angle x subjects
% Mean force, interarrival time, length and angle of each step for each
% person in one matrix. For each person(represented in the 3rd dimension),
% 1st column represents the Dt, 2nd column is the meanF, 3rd column is the
% length and 4th is the angle for the specific step (row)
%%%
% meanF : subjects x steps
% Mean force induced on each step (col) for each subject (row)
%%%
% Dt: subjects x steps
% Interarrival time between each step (col) for each subject (row).
%%%
% len: subjects x steps
% The euclidian distance between both hills (two points) on each step
%%%
% angle: subjects x steps
% The angle of each step along the gait horizontal direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X, Dt,meanF, len, angle] = computeAllDesiredVariables(force, time, x_coord, y_coord)
    
     for i=1:size(time,3) % All subjects
        meanTime(i,:) = nanmean(time(:,:,i));   % Compute mean time for each step excluding nan values
        
        % Calculate meanF
        % meanF is the mean force for each step excluding nan values
        meanF(i,:) = nanmean(force(:,:,i));
        
        % Init Dt, length, angle
        Dt(i,1) = nan;  % Dt1 cannot be computed
        len(i,1) = nan;   % length1 cannot be computed
        angle(i,1) = nan; % angle1 cannot be computed
        
        % Calculate Dt, length, angle
        for j=2:size(meanTime,2) % All steps
            
            % Clear first measurement when new measurement of same subject begins
            % If current x_coord is less than the previous, it means that
            % the subject started from the beggining, so measurements for
            % Dt, length and angle are not valid
            if x_coord(i,j) < x_coord(i,j-1)
               Dt(i,j)=nan;
               len(i,j)=nan;
               angle(i,j)=nan;
            else
                % Calculate Dt
                % t(j) : mean time of j step
                % t(j-1): mean time of j-1 step
                % Dt(j) = t(j)-t(j-1) : the interarrival time from j-1 to jth step
                Dt(i,j) = meanTime(i,j) - meanTime(i,j-1);

                % Calculate length and angle
                % dx, dy are the intervals of distance in x and y axis
                % length is calculated using euclidian distance
                % angle is calculated as the arctan of dy/dx
                dx = x_coord(i,j)-x_coord(i,j-1);
                dy = y_coord(i,j)-y_coord(i,j-1);
                len(i,j) = sqrt(dx.^2 + dy.^2 );
                angle(i,j) = atan2(dy,dx);
            end
        end
     end
    
    % Clear meanF measurements if there is no pair in time
    % When there is no pair in time, there is also none in length and angle
    meanF(isnan(Dt))=nan;
    
    % Shift right all nan values
    [~, I] = sort (isnan (Dt), 2);
    for i=1:size(Dt,1)
        Dt (i, :) = Dt (i, I(i, :));
        meanF(i,:) = meanF(i, I(i,:));
        len(i,:) = len(i, I(i,:));
        angle(i,:) = angle(i, I(i,:));
    end
       
    % Compute the X matrixes
    for i=1:size(Dt,1) % All subjects
        for j =1:size(Dt,2) % All steps
            X(j,1,i) = Dt(i,j);
            X(j,2,i) = meanF(i,j);
            X(j,3,i) = len(i,j);
            X(j,4,i) = angle(i,j);
        end
    end

end
