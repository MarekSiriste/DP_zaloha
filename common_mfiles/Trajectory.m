% -------- Generate trajectory ----------%
function x_r = Trajectory(models, initialProbabilities, transitionMatrix, T, simTime, init)
    % Generate trajectory using two models - HMM algorithm
    % inputs:
    %        - models - models used in HMM to generate traj
    %        - initialProbabilities - init probabilities of models 
    %        - transitionMatrix - Matrix of trasitions used for HMM algorithm
    %        - T - sampling period
    %        - simTime - time of simulations in seconds [s]
    %        - init: init state
    
    x = zeros(ceil((simTime/T)+1), size(models{1}.mu,2));
    
    % Init from input
    x(1,:) = init;
    
    % Choose init model to start
    currentModel = randsample(1:2, 1, true, initialProbabilities);
    model_act(1) = currentModel;
    
    % Generate init based on model
    if currentModel == 1 % First model is active
      curr_mod = models{1};
    else % second is active
      curr_mod = models{2};
    end
    x(2,:) = SimStep(curr_mod, x(1,:), T);
    model_act(2) = currentModel;
    
    scaled_index = 3;
    for t = T:T:(simTime)
         if currentModel == 1 % First model is active
             curr_mod = models{1};
         else % second is active
             curr_mod = models{2};
         end
         model_act(scaled_index) = currentModel;
         x(scaled_index,:) = SimStep(curr_mod, x(scaled_index-1,:), T);
         scaled_index = scaled_index + 1;
         % Choose new model based on transition matrix
         currentModel = randsample(1:2, 1, true, transitionMatrix(currentModel, :));
    end
    
    % Rearrange Vector:
    x_r = rearrangeVec(x);
%     plotData(x_r, simTime, T, model_act);
end


%-----------------------------------------------------%
%--------------- Simulate next Step ------------------%
%-----------------------------------------------------%
function x = SimStep(model, x, T)
    % model - current model to use
    % x - state vector
    % T = sampling period 

    mu_x = model.mu;
    F = model.F;
    Q = model.Q;
    [U,S,~] = svd(Q);
    sqrtP = U*S^(1/2);
    
    noise = sqrtP*randn(length(mu_x), 1) + mu_x'; % OR: mvnrnd(mu_x, Q);
    x = F * x(end,:)' + noise;
end

%-----------------------------------------------------%
%--------------- vector rearrangement ----------------%
%-----------------------------------------------------%
function x = rearrangeVec(a)
% The input vector a is [x_pos, x_vel, "x_acc", y_pos, y_vel, "y_acc", z_pos, z_vel, "z_acc"]
% For clarity it is better to reorganize the vector to [x_pos, y_pos, z_pos, x_vel,.....]
    if size(a,2) > 6
        % The acceleration data is also available
        % The vector is [x, df(x), ddf(x), y, df(y), ddf(y), z, df(z), ddf(z)]
        % Indeces of positions, velocities and accelerations
        ipos = [1, 4, 7];
        ivel = [2, 5, 8];
        iacc = [3, 6, 9];
        % Extracting from vector
        x(:, 1:3) = a(:, ipos);
        x(:, 4:6) = a(:, ivel);
        x(:, 7:9) = a(:, iacc);
    else
        % The vector is [x, df(x), y, df(y), z, df(z)]
        % Indeces of positions, velocities
        ipos = [1, 3, 5];
        ivel = [2, 4, 6];
        % Extracting position from vector
        x(:, 1:3) = a(:, ipos);
        x(:, 4:6) = a(:, ivel);
    end
end

%{
% Todo: just for testing
function plotData(x_r, simTime, T, modelActive)
% Function only for generating plots for thesis
time = 0:T:(simTime+T);
positions = x_r(:, 1:3);
velocities = x_r(:, 4:6);
accelerations = x_r(:, 7:9);

figure;
subplot(3,1,1);
plot(time, positions);
title('Positions');
xlabel('Time [s]');
ylabel('$[m]$');
legend('x', 'y', 'z');

subplot(3,1,2);
plot(time, velocities);
title('Velocities');
xlabel('Time [s]');
ylabel('$[m/s]$');
legend('x', 'y', 'z');

subplot(3,1,3);
plot(time, accelerations);
title('Accelerations');
xlabel('Time [s]');
ylabel('$[m/s^2]$');
legend('x', 'y', 'z');
saveas(gcf,'fig3_generated_trajectory','epsc')


% Figure for active model
figure;
stairs(time, modelActive, 'LineWidth', 2);
ylim([0.8 2.2]);
xlim([0 time(end)]);
yticks([1 2]);
yticklabels({'Singer', 'DWPA'});
title('Active model over time');
xlabel('Time [s]');
ylabel('Model');
grid on;

% Colors for each model
hold on;
for i = 1:length(time)-1
    if modelActive(i) == 1
        stairs(time(i:i+1), modelActive(i:i+1), 'b', 'LineWidth', 2);
    else
        stairs(time(i:i+1), modelActive(i:i+1), 'r', 'LineWidth', 2);
    end
end
saveas(gcf,'fig4_HMM','epsc')
end
%}