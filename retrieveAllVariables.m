%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- retrieveAllVariables(database)
% Retrieve all forces, times and coordinates from the database for each step
%
%%% Returns %%%
%%%
% time : 67x32x215 (max measurements x steps x subjects
% Times in which measurements occured(1st dim) for each step (2nd dim) and
% each subject (3rd dim).
%%%
% force: 67x32x215 (max measurements x steps x subjects)
% Force induced on each moment (1st dim) for each step (2nd dim) and each
% subject (3rd dim).
%%%
% x_coord: 215x32 (subjects x steps)
% The step coordinates in x axis for all subjects(rows) and steps (col)
%%%
% y_coord: 215x32 (subjects x steps)
% The step coordinates in y axis for all subjects(rows) and steps (col)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [time,force, x_coord, y_coord] = retrieveAllVariables(database)

    for i=1:size(database,1) % All subjects
        for j=1:size(database,2) % All Steps
            
            % Retrieve measurements if not empty, else put nan
            if ~isempty(database{i,j})
                force(:,j,i)=database{i,j}.force;
                time(:,j,i)=database{i,j}.time;
                x_coord(i,j)=database{i,j}.x_step;
                y_coord(i,j)=database{i,j}.y_step;
            else
                force(:,j,i)=nan;
                time(:,j,i)=nan;
                x_coord(i,j)=nan;
                y_coord(i,j)=nan;
            end
            
        end
    end
    
end