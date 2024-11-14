function [x_kk,P_kk]  = calcEstimate(x_traj, models, sensors, surf, T, initialProbabilities, transitionMatrix, Meas_type, MeasInterval)
% Function to calculate estimate via IMM algorithm
% Inputs:
%           - x_traj - simulated trajectory (true values)
%           - models - models which will be used for estimation process
%           - sensors - sensors config (scenario also - number of TDOA/AOA)
%           - surf - surface config
%           - T - sampling period
%           - initialProbabilities: init probab. of each model
%           - transitionMatrix: Markov matrix of transition between models
%           - MeasType - Linear/NonLinear (positions/TDOA + AOA)
%           - MeasInterval - interval between measurements that the filters consider
% Outputs:
%           - x_kk - estimated states
%           - P_kk - estimated covariance matrices

    if strcmp(Meas_type, 'Linear') || strcmp(Meas_type, 'linear')
       % Linear scenario % 
       
       % Linear measurements
       % Calc linear measurements
       [Linear_meas, R] = calcLinearMeasurements(sensors, surf, x_traj); 
       
        for r=1:numel(models) % Init for each model
            [x_init(:,r), P_init(:,:,r)] = calcInit(models{r}, Linear_meas, R, T);
        end
        
        % IMM 
        [x_kk(:,1),P_kk(:,:,1),modeProb(1,:),x_kk_j(:,:,1),P_kk_j(:,:,:,1)] = IMM(initialProbabilities, transitionMatrix, x_init, P_init, models, sensors, Linear_meas(1,:), R(:,:,1), Meas_type, 1, MeasInterval);
         n = length(Linear_meas);
        for i = 2:n
            [x_kk(:,i),P_kk(:,:,i),modeProb(i,:),x_kk_j(:,:,i),P_kk_j(:,:,:,i)] = IMM(modeProb(i-1,:),transitionMatrix, x_kk_j(:,:,i-1),P_kk_j(:,:,:,i-1), models, sensors, Linear_meas(i,:), R(:,:,i), Meas_type, i, MeasInterval);
        end
    else
        % NonLinear scenario %
        
        % Calc linear measurements for init 
        %       - only first three/two measurements needed
        [Linear_meas, R_lin] = calcLinearMeasurements(sensors, surf, x_traj(1:3,:)); 

        for r=1:numel(models) % Init for each model
            [x_init(:,r), P_init(:,:,r)] = calcInit(models{r}, Linear_meas, R_lin, T);
        end

        % Calc Non-linear measurements
        [DeltaT, R] = calcNonLinearMeasurements(sensors, x_traj);

        % IMM 
        [x_kk(:,1),P_kk(:,:,1),modeProb(1,:),x_kk_j(:,:,1),P_kk_j(:,:,:,1)] = IMM(initialProbabilities, transitionMatrix, x_init, P_init, models, sensors, DeltaT(1,:), R, Meas_type, 1, MeasInterval);
        n = length(DeltaT);
        for i = 2:n
            [x_kk(:,i),P_kk(:,:,i),modeProb(i,:),x_kk_j(:,:,i),P_kk_j(:,:,:,i)] = IMM(modeProb(i-1,:),transitionMatrix, x_kk_j(:,:,i-1),P_kk_j(:,:,:,i-1), models, sensors, DeltaT(i,:), R, Meas_type, i, MeasInterval);
        end
        
    end
    %plotComparisons(x_traj, Linear_meas); % Just for figures to thesis

    return
end

%{
% Just for testing
function plotComparisons(x_traj, Linear_meas)
    % This function is just for generating figures
    x_r = x_traj;
    est = Linear_meas';
    
    % Calculate squared differences
    squaredDifferences = (x_r(:,1:3) - est').^2;

    figure;
    subplot(3,1,1);
    plot(squaredDifferences(:,1), 'b-', 'LineWidth', 2);
    ylabel('MSE x');
    title('MSE Between True Positions and Estimates');
    legend('MSE x'); 
    hold off;

    subplot(3,1,2);
    plot(squaredDifferences(:,2), 'g-', 'LineWidth', 2);
    ylabel('MSE y');
    legend('MSE y'); 
    hold off;

    subplot(3,1,3);
    plot(squaredDifferences(:,3), 'r-', 'LineWidth', 2);
    ylabel('MSE z');
    xlabel('Sample Index');
    legend('MSE z'); 
    hold off;
    saveas(gcf,'fig5_linear_est1_MSE','epsc')
        

    % Figure for subplots
    figure;
    subplot(3,1,1); % Subplot for x-axis
    plot(x_r(:,1), 'b-', 'LineWidth', 2);
    hold on;
    plot(est(1,:), 'r--', 'LineWidth', 2);
    ylabel('X Position');
    title('Comparison of True Positions and Linear Measurements');
    legend('True Position', 'Linear Measurement');
    hold off;

    subplot(3,1,2); % Subplot for y-axis
    plot(x_r(:,2), 'b-', 'LineWidth', 2);
    hold on;
    plot(est(2,:), 'r--', 'LineWidth', 2);
    ylabel('Y Position');
    legend('True Position', 'Linear Measurement');
    hold off;

    subplot(3,1,3); % Subplot for z-axis
    plot(x_r(:,3), 'b-', 'LineWidth', 2);
    hold on;
    plot(est(3,:), 'r--', 'LineWidth', 2);
    ylabel('Z Position');
    legend('True Position', 'Linear Measurement');
    xlabel('Sample Index');
    hold off;
    saveas(gcf,'fig5_linear_est1','epsc')
    
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
    plot3(est(1,:), est(2,:), est(3,:),'r', 'LineWidth', 2);
    hold off;
    legend('True emitter position', 'Starting Point', 'Ending Point','Linear Measurements', 'FontSize', 8, 'Location', 'best');
    view(3); 
    saveas(gcf,'fig6_linear_est2','epsc')

end
%}
