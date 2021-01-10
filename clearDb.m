%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- clearDb(database)
% Clear trash in database & fill force-time matrixes to match indices
%
%%% Returns %%%
%%%
% database : subjects x steps
% The new database - clear of junk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function database = clearDb(database)

    % clear empty spaces in database
    for i=1:size(database,1) % All subjects
        if isempty(database{i})
            database(i,:)=[];
        end
    end

    % Clear negative and zero force values
    for i=1:size(database,1) % All subjects
        for j=1:size(database,2) % All Steps
            if ~isempty(database{i,j})
                database{i,j}.time(database{i,j}.force <= 0) = nan;
                database{i,j}.force(database{i,j}.force <= 0) = nan;
            end
        end
    end
    
    % identify max force-time indices
    max_length = 0;
    for i=1:size(database,1) % All subjects
        for j=1:size(database,2) % All Steps
            if ~isempty(database{i,j})
                if length(database{i,j}.force) > max_length
                    max_length=length(database{i,j}.force);
                end
            end
        end
    end

    % fill with nan all forces and times up to max indices
    for i=1:size(database,1) % All subjects
        for j=1:size(database,2) % All Steps
            if ~isempty(database{i,j})
                database{i,j}.force(end+1:max_length)=nan;
                database{i,j}.time(end+1:max_length)=nan;
            end
        end
    end
    
end