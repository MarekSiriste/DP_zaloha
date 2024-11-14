function FigureEstimateQuality(trajectories, estimates, nSim, T, zoom)
    % Function will figure chosen simulation/Monte carlo sample
    % Inputs:
    %           - nSim - number of simulation to figure
    %           - trajectories - struct of trajectories
    %           - estimates - struct of estimates
    %           - T - sampling period
    %           - zoom = [min max] - interval of zoomed figure 
    % Print: comparison of estimates and truth for each axis, RMSE for each axis and anees
    
    x_r = trajectories{nSim};
    est = real(estimates{nSim}.x');
    
    t = (0:size(est,1)-1) * T;
    
    % Models are 9 dimensional or 6 dimensional
    if size(est, 2) == 9
        pos_index = [1 4 7];    
    else
        pos_index = [1 3 5];
    end
    
    %rmse_values = sqrt(mean(errors.^2, 2));  % RMSE at each time step
    N = size(est,1);
    errors = zeros(N, 3);
    rmse_values = zeros(N,3);
    if isfield(estimates{nSim}, 'P')
        anees_values = zeros(size(N, 1), 1);
        for i = 1:N
            errors(i, :) = est(i, pos_index) - x_r(i, 1:3);
            P_inv = inv(estimates{nSim}.P(pos_index,pos_index,i));
            anees_values(i) = (errors(i, :) * P_inv * errors(i, :)') / ((length(pos_index)));
            rmse_values(i,:) = sqrt((errors(i, :)).^2)./N;
        end
    else
        anees_values = nan(size(errors, 1), 1); % No covariance data available
    end
    
    anees_values = real(anees_values);
    rmse_values = real(rmse_values);
    
    % Figure for subplots
    figure;
    subplot(3,1,1); % Subplot for x-axis
    plot(t, x_r(:,1), 'b-', 'LineWidth', 2);
    hold on;
    grid on;
    plot(t,est(:,pos_index(1)), 'r--', 'LineWidth', 2);
    ylabel('x Position [m]');
    title('Comparison of True Positions and Estimated positions');
    legend('True', 'Estimate');
    hold off;
    xlim([1 t(end)]);
    % Zoom na detail
    ax = gca;
    ax_pos = ax.Position;
    zoom_ax = axes('Position', [ax_pos(1)+0.05 ax_pos(2)+0.05 ax_pos(3)*0.35 ax_pos(4)*0.25]);
    plot(t, x_r(:,1), 'b-', 'LineWidth', 2);
    hold on;
    grid on;
    plot(t, est(:,pos_index(1)), 'r--', 'LineWidth', 2);
    xlim([zoom(1) zoom(2)])

    subplot(3,1,2); % Subplot for y-axis
    plot(t, x_r(:,2), 'b-', 'LineWidth', 2);
    hold on;
    grid on
    plot(t, est(:,pos_index(2)), 'r--', 'LineWidth', 2);
    ylabel('y Position [m]');
    legend('True', 'Estimate');
    xlim([1 t(end)]);
    hold off;
    % Zoom na detail
    ax = gca;
    ax_pos = ax.Position;
    zoom_ax = axes('Position', [ax_pos(1)+0.05 ax_pos(2)+0.05 ax_pos(3)*0.35 ax_pos(4)*0.25]);
    plot(t, x_r(:,2), 'b-', 'LineWidth', 2);
    hold on;
    grid on
    plot(t, est(:,pos_index(2)), 'r--', 'LineWidth', 2);
    xlim([zoom(1) zoom(2)])

    subplot(3,1,3); % Subplot for z-axis
    plot(t, x_r(:,3), 'b-', 'LineWidth', 2);
    hold on;
    grid on;
    plot(t, est(:,pos_index(3)), 'r--', 'LineWidth', 2);
    ylabel('z Position [m]');
    legend('True', 'Estimate');
    xlabel('t [s]');
    xlim([1 t(end)]);
    hold off;
    % Zoom na detail
    ax = gca;
    ax_pos = ax.Position;
    zoom_ax = axes('Position', [ax_pos(1)+0.05 ax_pos(2)+0.05 ax_pos(3)*0.35 ax_pos(4)*0.25]);
    plot(t, x_r(:,3), 'b-', 'LineWidth', 2);
    hold on;
    grid on
    plot(t, est(:,pos_index(3)), 'r--', 'LineWidth', 2);
    xlim([zoom(1) zoom(2)])
    saveas(gcf,'fig10_compare_est_true','epsc')
    
    % Plotting RMSE of each axis
    figure;
    subplot(3,1,1);
    plot(t, rmse_values(:,1), 'b-', 'LineWidth', 2);
    title('Root Mean Square Error Over Time');
    ylabel('x [m]');
    xlim([1 t(end)]);
    grid on;

    subplot(3,1,2);
    plot(t, rmse_values(:,2), 'b-', 'LineWidth', 2);
    title('Root Mean Square Error Over Time');
    ylabel('y [m]');
    xlim([1 t(end)]);
    grid on;
    
    subplot(3,1,3);
    plot(t, rmse_values(:,3), 'b-', 'LineWidth', 2);
    title('Root Mean Square Error Over Time');
    ylabel('z [m]');
    xlabel('t [s]');
    grid on;
    xlim([1 t(end)]);
    saveas(gcf,'fig11_rmse_each_axis','epsc')
    
    
    figure;
    grid on;
    title('(A)NEES Over Time');
    xlabel('t [s]');
    ylabel('(A)NEES');
    hold on;

    red = [1, 0, 0];
    blue = [0, 0, 1];
    green = [0, 1, 0];

    % Initialize color array
    colors = zeros(N, 3); 

    % Assign colors based on ANEES values
    for i = 1:N
        if anees_values(i) > 1.1
            colors(i, :) = blue;
        elseif anees_values(i) >= 0.9
            colors(i, :) = green;
        else
            colors(i, :) = red; 
        end
    end

    scatter(t, anees_values, 36, colors, 'filled', 'HandleVisibility', 'off');
%     textLocation = [t(round(N/2)) max(anees_values)]; 
%     textBorder = {'BackgroundColor', 'white', 'Margin', 3, 'EdgeColor', 'black'};
%     text(textLocation(1), textLocation(2), sprintf('ANEES: %.2f', anees), ...
%      'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', textBorder{:});

    % Create a custom legend
    hold on;
    plot(t, anees_values, '-', 'Color', 'black', 'LineWidth', 1,'HandleVisibility', 'off')
    one_vec = ones(N,1);
    plot(t, one_vec, '--', 'Color', '#EDB120', 'LineWidth', 2,'HandleVisibility', 'off')
    scatter([], [], 36, blue, 'filled', 'DisplayName', 'Overconfident');
    scatter([], [], 36, green, 'filled', 'DisplayName', 'Confident');
    scatter([], [], 36, red, 'filled', 'DisplayName', 'Underconfident');
    legend('show', 'Location', 'best');
    hold off;
    xlim([1 t(end)]);
    set(gca, 'YScale', 'log');
    saveas(gcf,'fig12_anees','epsc')
end