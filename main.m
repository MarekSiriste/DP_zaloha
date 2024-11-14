%% --------------------------------------------------------------- %%
%  --------------------- Marek Siřiště ---------------------------  %
%  ---------------------------------------------------------------  %
clc
clear all
close all

set(0,'defaulttextInterpreter','latex') %latex axis labels
set(0,'defaultaxesFontSize',12)
set(0, 'DefaultFigureRenderer', 'painters');
format shortG 
warning('off','all')

% Add path to the folder with common mfiles
common_files_path = 'common_mfiles';
addpath(genpath(common_files_path));


%  ---------------------------------------------------------------- %
%  ----------------- 1.) Configuration ---------------------------  %
%  ---------------------------------------------------------------  %

% a) Config area surface and receivers

% for configuration it is necessary to modify ConfigLoader.m
%   - in this m-file new function with config can be added
%   - This configuration is only for surface

% Load the configuration
% TODO - more configuration can be added
%      - LoadConfig1 = configuration for testing
%      - LoadCofig2 = testing (faster)
%      - LoadCofig3 = era
config = ConfigLoader(2);
% Accessing the configuration data
surfaceLimits = [config.surface.xMin, config.surface.xMax, config.surface.yMin, config.surface.yMax, config.surface.zMin, config.surface.zMax];
quantization = [config.surface.xQuantization, config.surface.yQuantization, config.surface.zQuantization];
sensors = config.sensors;
surf = config.surface;
init_states = config.initStates;


% Print info about configuration of surface and receivers
printInfoConfig(config);

% b) Config models and simulations

% Configuration for Models
%   - To modify models - ConfigModels.m needs to be changed
[NCV, DWNA, DWPA, SINGER, initialProbabilities, transitionMatrix, T] = ConfigModels();

% Time of each simulation in seconds
simTime = 10;

% Number of Monte Carlo simulations
nMC = length(init_states(:,1:9));

% Denote the frequency at which measurements are sampled and used
measInt = 1; % Each 10-th measurement

bChooseMdl = 1; % 1 for [DWPA, SINGER], 0 for [NCV, DWNA]

if bChooseMdl
    models = {SINGER, DWPA};
    start_pos = init_states;
    pos= [1 4 7];
else
    models = {NCV, DWNA};
    % These models do not have accelerations
    start_pos = init_states(:, [1, 2, 4, 5, 7, 8]);
    pos = [1 3 5];
end

%  ---------------------------------------------------------------  %
%  ---------------------------------------------------------------  %

%%
%  ---------------------------------------------------------------- %
%  --------------- 2.) Generating Trajectory ---------------------  %
% In this section:
%   - trajectory is generated using Trajectory.m
%   - each trajectory is stored for later estimates

% Initialize arrays to store
trajectories = cell(1, nMC);

for j = 1:nMC
    x_traj = Trajectory(models, initialProbabilities, transitionMatrix, T, simTime, start_pos(j,:));
    trajectories{j} = x_traj;
end

%  ---------------------------------------------------------------  %
%  ---------------------------------------------------------------  %

%%
% Main
% {
% Main simulation
%  --------------- 3.) Estimate - KF, EKF ------------------------  %
%  ---------------------------------------------------------------  %

% Info:
%   - Several scenarios included
%   - scenario can be modulated by changing:
%                   - sensors config (TDOA/AOA availability)
%                   - surf config
%                   - models config
%                   - MeasType - Linear/Nonlinear

% Main 6 scenarios
%  ---------------------------------------------------------------- %
%  -----------Scenario 1 - KF with linear meas (TDOA)-------------  %
% Initialize array to store
estimates_1 = cell(1, nMC);

% Config scenario:
% Turn ON all TOA
sensors(1).TOA = true;
sensors(2).TOA = true;
sensors(3).TOA = true;
sensors(4).TOA = true;
% Turn off all AOA
sensors(1).AOA = false;
sensors(2).AOA = false;
sensors(3).AOA = false;
sensors(4).AOA = false;
Meas_type = 'linear'; % Use linear measurements
for j = 1:nMC
    
    x_traj = trajectories{j};
 
    [x_kk,P_kk]  = calcEstimate(x_traj, models, sensors, surf, T, initialProbabilities, transitionMatrix, Meas_type, measInt);
    
    estimates_1{j}.x = x_kk;
    estimates_1{j}.P = P_kk;
end
% FigureSimulation(trajectories, estimates_1, 1)

%  ---------------------------------------------------------------- %
%  -----------Scenario 2 - KF with linear meas (TDOA + AOA)-------  %

estimates_2 = cell(1, nMC);

% Config scenario:
% Turn ON all TOA
sensors(1).TOA = true;
sensors(2).TOA = true;
sensors(3).TOA = true;
sensors(4).TOA = true;
% Turn on all AOA
sensors(1).AOA = true;
sensors(2).AOA = true;
sensors(3).AOA = true;
sensors(4).AOA = true;
Meas_type = 'linear'; % Use linear measurements
for j = 1:nMC
    x_traj = trajectories{j};
 
    [x_kk,P_kk]  = calcEstimate(x_traj, models, sensors, surf, T, initialProbabilities, transitionMatrix, Meas_type, measInt);
    
    estimates_2{j}.x = x_kk;
    estimates_2{j}.P = P_kk;
end
%FigureSimulation(trajectories, estimates_2, 1)

%  ---------------------------------------------------------------- %
%  -----------Scenario 3 - KF with 2 TDOAs and 3 AoAs ------------  %
estimates_3 = cell(1, nMC);

% Config scenario:
% Turn off all AOA
sensors(1).AOA = true;
sensors(2).AOA = true;
sensors(3).AOA = true;
sensors(4).AOA = false;
% Turn off one TDOA
sensors(1).TOA = true;
sensors(2).TOA = true;
sensors(3).TOA = true;
sensors(4).TOA = false;
Meas_type = 'linear'; % Use linear measurements
for j = 1:nMC
    
    x_traj = trajectories{j};
 
    [x_kk,P_kk]  = calcEstimate(x_traj, models, sensors, surf, T, initialProbabilities, transitionMatrix, Meas_type, measInt);
    
    estimates_3{j}.x = x_kk;
    estimates_3{j}.P = P_kk;
end
% FigureSimulation(trajectories, estimates_3, 1)

%  ---------------------------------------------------------------- %
%  -----------Scenario 4 - EKF with 3 TDOAs and 0 AOAs -----------  %

estimates_4 = cell(1, nMC);

% Config:
% Turn ON all AOA
sensors(1).AOA = false;
sensors(2).AOA = false;
sensors(3).AOA = false;
sensors(4).AOA = false;
% Turn ON all TOA
sensors(1).TOA = true;
sensors(2).TOA = true;
sensors(3).TOA = true;
sensors(4).TOA = true;
Meas_type = 'nonlinear'; % Use linear measurements
for j = 1:nMC
    
    x_traj = trajectories{j};
 
    [x_kk,P_kk]  = calcEstimate(x_traj, models, sensors, surf, T, initialProbabilities, transitionMatrix, Meas_type, measInt);
    
    estimates_4{j}.x = x_kk;
    estimates_4{j}.P = P_kk;
end
FigureSimulation(trajectories, estimates_4, 1)

%  ---------------------------------------------------------------- %
%  -----------Scenario 5 - EKF with 3 TDOAs and 4 AOAs -----------  %

estimates_5 = cell(1, nMC);

% Config:
% Turn on all AOA
sensors(1).AOA = true;
sensors(2).AOA = true;
sensors(3).AOA = true;
sensors(4).AOA = true;
% Turn ON all TOA
sensors(1).TOA = true;
sensors(2).TOA = true;
sensors(3).TOA = true;
sensors(4).TOA = true;
Meas_type = 'nonlinear'; % Use linear measurements
for j = 1:nMC
    
    x_traj = trajectories{j};
 
    [x_kk,P_kk]  = calcEstimate(x_traj, models, sensors, surf, T, initialProbabilities, transitionMatrix, Meas_type, measInt);
    
    estimates_5{j}.x = x_kk;
    estimates_5{j}.P = P_kk;
end
% FigureSimulation(trajectories, estimates_5, 1)

%  ---------------------------------------------------------------- %
%  -----------Scenario 6 - EKF with 2 TDOAs and 3 AOAs -----------  %

estimates_6 = cell(1, nMC);
% Config:
% Turn on all AOA
sensors(1).AOA = true;
sensors(2).AOA = true;
sensors(3).AOA = true;
sensors(4).AOA = false;
% Turn off one TDOA
sensors(1).TOA = true;
sensors(2).TOA = true;
sensors(3).TOA = true;
sensors(4).TOA = false;
Meas_type = 'nonlinear'; % Use linear measurements
for j = 1:nMC
    
    x_traj = trajectories{j};
 
    [x_kk,P_kk]  = calcEstimate(x_traj, models, sensors, surf, T, initialProbabilities, transitionMatrix, Meas_type, measInt);
    
    estimates_6{j}.x = x_kk;
    estimates_6{j}.P = P_kk;
end
% FigureSimulation(trajectories, estimates_6, 1)


save('Test_final_layour_era_NCV_DWNA_all_meas.mat');
%  ---------------------------------------------------------------- %
%  ---------------------------------------------------------------- %

%%
%  ---------------------------------------------------------------- %
%  --------------- 5.) Visualization part ------------------------  %
%  ---------------------------------------------------------------  %

% Create Eval matrix for trajectories
eval_matrix = CreateEvalMatrix(surf, quantization, trajectories);

% Treshold value for RMSE in [m]
tres = 20;

visualizeEstimateQuality(eval_matrix, estimates_1, trajectories, surf, quantization, 'rmse', tres);
visualizeEstimateQuality(eval_matrix, estimates_2, trajectories, surf, quantization, 'rmse', tres);
visualizeEstimateQuality(eval_matrix, estimates_3, trajectories, surf, quantization, 'rmse', tres);
visualizeEstimateQuality(eval_matrix, estimates_4, trajectories, surf, quantization, 'rmse', tres);
visualizeEstimateQuality(eval_matrix, estimates_5, trajectories, surf, quantization, 'rmse', tres);
visualizeEstimateQuality(eval_matrix, estimates_6, trajectories, surf, quantization, 'rmse', tres);

% 2D
visualizeEstimateQuality2D_mean(eval_matrix, estimates_1, trajectories, surf, sensors, quantization, tres);
visualizeEstimateQuality2D_mean(eval_matrix, estimates_2, trajectories, surf, sensors, quantization, tres);
visualizeEstimateQuality2D_mean(eval_matrix, estimates_3, trajectories, surf, sensors, quantization, tres);
visualizeEstimateQuality2D_mean(eval_matrix, estimates_4, trajectories, surf, sensors, quantization, tres);
visualizeEstimateQuality2D_mean(eval_matrix, estimates_5, trajectories, surf, sensors, quantization, tres);
visualizeEstimateQuality2D_mean(eval_matrix, estimates_6, trajectories, surf, sensors, quantization, tres);

%}

% Testing 
%{
% For testing

estimates_4 = cell(1, nMC);

% Config:
% Turn off all AOA
sensors(1).AOA = true;
sensors(2).AOA = true;
sensors(3).AOA = true;
sensors(4).AOA = true;
Meas_type = 'nonlinear'; % Use linear measurements
for j = 1:nMC
    
    x_traj = trajectories{j};
    
    [x_kk,P_kk]  = calcEstimate(x_traj, models, sensors, surf, T, initialProbabilities, transitionMatrix, Meas_type, measInt);
    
    estimates_4{j}.x = x_kk;
    estimates_4{j}.P = P_kk;
end
FigureSimulation(trajectories, estimates_4, 1)

%%
%  ---------------------------------------------------------------- %
%  --------------- 5.) Visualization part ------------------------  %
%  ---------------------------------------------------------------  %

% Create Eval matrix for trajectories
eval_matrix = CreateEvalMatrix(surf, quantization, trajectories);
eval_matrix2D = CreateEvalMatrix2D(surf, quantization, trajectories);
tres = 20;

%%
close all
visualizeEstimateQuality2(eval_matrix, estimates_4, trajectories, surf, quantization, 'rmse', tres);

% 2D visualization of quality
% Last parameter is for ,,choosing the layer of z axis" to figure
% visualizeEstimateQuality2D(eval_matrix, estimates_4, trajectories, surf, sensors, quantization,0)
visualizeEstimateQuality2D_mean(eval_matrix, estimates_4, trajectories, surf, sensors, quantization, tres);
% visualizeEstimateQuality2D_mean2(eval_matrix2D, estimates_4, trajectories, surf, sensors, quantization, tres)

%}


%% Visualization
%{
%  ---------------------------------------------------------------- %
%  --------------- 5.) Visualization part ------------------------  %
%  ---------------------------------------------------------------  %

% a) Visualization in 3D - comparison of each scenario
%load('testing_more_data.mat')

% Quant modification
quantization = [1000 1000 1000];

% Create Eval matrix for trajectories
eval_matrix = CreateEvalMatrix(surf, quantization, trajectories);

% Treshold value for RMSE in [m]
tres = 20;

% visualizeEstimateQuality(eval_matrix, estimates_1, trajectories, surf, quantization, 'rmse', tres);
% visualizeEstimateQuality(eval_matrix, estimates_2, trajectories, surf, quantization, 'rmse', tres);
% visualizeEstimateQuality(eval_matrix, estimates_3, trajectories, surf, quantization, 'rmse', tres);
% visualizeEstimateQuality(eval_matrix, estimates_4, trajectories, surf, quantization, 'rmse', tres);
% visualizeEstimateQuality(eval_matrix, estimates_5, trajectories, surf, quantization, 'rmse', tres);
% visualizeEstimateQuality(eval_matrix, estimates_6, trajectories, surf, quantization, 'rmse', tres);

% 2D
% visualizeEstimateQuality2D_mean(eval_matrix, estimates_1, trajectories, surf, sensors, quantization, tres);
% visualizeEstimateQuality2D_mean(eval_matrix, estimates_2, trajectories, surf, sensors, quantization, tres);
% visualizeEstimateQuality2D_mean(eval_matrix, estimates_3, trajectories, surf, sensors, quantization, tres);
% visualizeEstimateQuality2D_mean(eval_matrix, estimates_4, trajectories, surf, sensors, quantization, tres);
% visualizeEstimateQuality2D_mean(eval_matrix, estimates_5, trajectories, surf, sensors, quantization, tres);
% visualizeEstimateQuality2D_mean(eval_matrix, estimates_6, trajectories, surf, sensors, quantization, tres);


% visualizeEstimateQualityLayers(eval_matrix, estimates_6, trajectories, surf, quantization, tres, 0, 0)

%}
