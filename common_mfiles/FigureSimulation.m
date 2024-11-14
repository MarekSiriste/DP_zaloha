function FigureSimulation(trajectories, estimates, nSim)
    % Function will figure chosen simulation
    % Inputs:
    %           - nSim - number of simulation to figure
    %           - trajectories - struct of trajectories
    %           - estimates - struct of estimates
    
    
    x_r = trajectories{nSim};
    est = estimates{nSim}.x;
    
    % Models are 9 dimentional or 6 dimentional
    if length(est(:,1)) == 9
        pos= [1 4 7];
    else
        pos = [1 3 5];
    end

    figure;
    plot3(x_r(:,1), x_r(:,2), x_r(:,3),'b--', 'LineWidth', 2);
    hold on;
    grid on;
    xlabel('X [m]');
    ylabel('Y [m]');
    zlabel('Z [m]');
    title('3D Position Plot');
    % Highlight starting and ending points
    plot3(x_r(1,1), x_r(1,2), x_r(1,3), 'ro', 'MarkerSize', 10, 'LineWidth', 2); % starting point
    plot3(x_r(end,1), x_r(end,2), x_r(end,3), 'go', 'MarkerSize', 10, 'LineWidth', 2); % ending point
    plot3(est(pos(1),:), est(pos(2),:), est(pos(3),:),'r', 'LineWidth', 2);
    hold off;
    legend('Pos', 'Starting Point', 'Ending Point','IMM - Kalman', 'FontSize', 8, 'Location', 'best');
    view(3); 
    %saveas(gcf,'fig7_IMM','epsc')
end

