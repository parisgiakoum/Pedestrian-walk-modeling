%% Plotting some steps from database to showcase how length and angle are calculated
%% Clear junk, retrieve force-time-position measurements and find meanF-Dt-length-angle for each step of each subject

database=load('steps_database').database_passi;

database = clearDb(database);

[time,force, x_coord, y_coord] = retrieveAllVariables(database);

%%
clear xpos ypos
for i=1:3
    if ~isempty(database{2,i})
        xpos(i) = database{2,i}.x_step;
        ypos(i) = database{2,i}.y_step;
    end
end

figure;
plot(xpos,ypos,'-or','MarkerEdgeColor','k')

ylim([-0.2 0.2])
xlim([0 2])
xlabel('x\_coord')
ylabel('y\_coord')
title('Calculation of length and angle')

text (xpos (1), ypos (1)+0.015, num2str (1), 'Color','r')
for k = 2: length (xpos)
    line([xpos(k-1),xpos(k-1)], [ypos(k-1),ypos(k)], 'Color', 'm', 'LineWidth', 1, 'LineStyle', '--');
    text (xpos(k-1)+0.01,ypos (k-1) + (ypos(k)-ypos(k-1))/2, strcat('y_',num2str(k-1),'_-_',num2str(k)), 'Color','m')
    line([xpos(k-1),xpos(k)], [ypos(k-1),ypos(k-1)], 'Color', 'b', 'LineWidth', 1,  'LineStyle', '--');
    text (xpos(k-1)+(xpos(k)-xpos(k-1))/2, ypos (k-1)+0.012, strcat('x_',num2str(k-1),'_-_',num2str(k)), 'Color','b')
    
    if ypos(k) < ypos(k-1)
        text (xpos (k), ypos (k)-0.015, num2str (k), 'Color','r')
        text (xpos(k-1)+(xpos(k)-xpos(k-1))/2, ypos (k-1)+ (ypos(k)-ypos(k-1))/2 + 0.015, strcat('length_',num2str(k-1),'_-_',num2str(k)), 'Color','r')
    else
        text (xpos (k), ypos (k)+0.015, num2str (k), 'Color','r')
        text (xpos(k-1)+(xpos(k)-xpos(k-1))/2, ypos (k-1)+ (ypos(k)-ypos(k-1))/2 - 0.015, strcat('length_',num2str(k-1),'_-_',num2str(k)), 'Color','r')
    end
end
