function [TDOA, AOA, Linear_Meas, R_lin] = calcMeasurements(sensors, surf, true_pos)
    % Compute measurements

    % Inputs:
    %          - sensors: location of each sensor and its characteristics
    %          - true_pos: true positions (trajectory)
    %          - surf: surface information

    % Outputs:
    %          - TDOA - computed TDOAs - non-linear
    %          - AOA - computed AOAs - non-linear
    %          - Linear_Meas - linear measurements (positions)
    %          - R_lin - linear measurements covariance matrix

    % TODO - crossing the boundaries of surface is not solved yet

    % ------------------------------------------------------------------------------- %
    % ----------------------------- Configuring ------------------------------------- %
    % ------------------------------------------------------------------------------- %

    nReceivers = numel(sensors); % number of receivers
    
    % Initialize vectors for variances
    variances_TDOA = [];
    variances_AOA = [];
    % Number of TDOA and AOA available
    nTDOA = 0;
    nAOA  = 0;
    
    % Collect variances
    for i = 2:nReceivers
        if sensors(i).TOA % Need to clarify if the AOA is available
            nTDOA = nTDOA + 1;
            variances_TDOA = [variances_TDOA; sensors(i).TOAVariance];
        end
    end
    for i = 1:nReceivers
        if sensors(i).AOA % Need to clarify if the AOA is available
            nAOA = nAOA + 2;
            variances_AOA = [variances_AOA; sensors(i).AOAVariance; sensors(i).AOAVariance];
        end
    end

    % Antenna noise parameters (TDOA + AOA noise)
    Qv = diag([variances_TDOA; variances_AOA].^2); % Covariance matrix of noise 

    % Load surface info
    quantization = [surf.xQuantization, surf.yQuantization, surf.zQuantization];
    xMin = surf.xMin;
    xMax = surf.xMax;
    yMin = surf.yMin;
    yMax = surf.yMax;
    zMin = surf.zMin;
    zMax = surf.zMax;
    

    % Anonymous functions 
    h = @(emitter) calcMeas(sensors, emitter, nTDOA, nAOA); % measurement model
    J = @(emitter) assembleJacobian(sensors, emitter, nTDOA, nAOA); % Jacobian

    % Declaring of TDOAs and AOAs
    TDOA = zeros(length(true_pos),nTDOA);
    AOA = zeros(length(true_pos),(nAOA)); % it contains azimut and elevation

    % ------------------------------------------------------------------------------- %
    % ------------------------------------------------------------------------------- %

    % ------------------------------------------------------------------------------- %
    % ----------------------------- Initial condition ------------------------------- %
    % ------------------------------------------------------------------------------- %

    % Emitter true initial position
    truth = true_pos(1,1:3)'; % location of emitter

    %DeltaT = TDOAModel(sensors, truth, nTDOA);
    
    % Simulate the TDOA measurement
    TDOA(1,:) = TDOAModel(sensors, truth, nTDOA);
    
    % Simulate AOA
    AOA(1,:) = AOAModel(sensors, truth, nAOA)';
    DeltaT = [TDOA(1,:)'; AOA(1,:)'];

    % Function to be optimized
    sqfunc = @(e) (DeltaT - h(e))./sqrt(diag(Qv)); % function of which a norm is minimized
    sqfunc2 = @(e) func2outputs(e, sqfunc, J);
    
    % adjust quantization for finding init condition
    quant = 3;
    % Surface - needed for brute force method
    [X,Y,Z] = meshgrid(xMin:quantization(1) * quant:xMax, yMin:quantization(2) * quant:yMax, zMin:quantization(3) * quant:zMax);
    
    % Brute force method for finding initial condition
    funVal = zeros(size(X));
    for i = 1:numel(X)
        pos = [X(i), Y(i), Z(i)]';
        funVal(i) = sum((DeltaT - h(pos)).^2); 
    end
    % Find the minimum value and its index
    [val, ind] = min(funVal(:));
    e0 = [X(ind), Y(ind), Z(ind)]'; % This is your initial estimate based on brute-force

    % Nonlinear squares - using init condition:

    % Optimization with lsqnonlin
    lb = [xMin; yMin; zMin]; % borders
    ub = [xMax; yMax; zMax];
    options = optimoptions('lsqnonlin', 'Display', 'off', 'CheckGradients', true); % Gradient options
    % Estimation of position:
    Linear_Meas(1,:)= lsqnonlin(sqfunc2, e0, lb, ub, options);

    %Output for debug:
    %fprintf('\tsolution: [%.0f %.0f %.0f]\n\ttruth:    [%.0f %.0f %.0f]\n', emEst(1), emEst(2), emEst(3), truth(1), truth(2), truth(3))

    % Covariance matrix of Measurement
    R_lin(:,:,1) = pinv(J(Linear_Meas(1,:)')) * Qv * pinv(J(Linear_Meas(1,:)'))'; % Covariance matrix calculation
    % ------------------------------------------------------------------------------- %
    % ------------------------------------------------------------------------------- %

    % ------------------------------------------------------------------------------- %
    % ----------------------- TDOA and estimate calculation ------------------------- %
    % ------------------------------------------------------------------------------- %

    searchRadius = [quantization(1),quantization(2), quantization(3)]; % Search radius of next ,,estimate"

    for i= 2:length(true_pos)

        truth = true_pos(i,1:3)'; % location of emitter

        % Todo - min/max borders given according to latest measurement
        % discuss
        lb = [Linear_Meas(i-1,1) - searchRadius(1); Linear_Meas(i-1,2) - searchRadius(2); Linear_Meas(i-1,3) - searchRadius(3)];
        ub = [Linear_Meas(i-1,1) + searchRadius(1); Linear_Meas(i-1,2) + searchRadius(2); Linear_Meas(i-1,3) + searchRadius(3)];
        
        %DeltaT = TDOAModel(sensors, truth, nTDOA);

        % Calculate TDOA
        TDOA(i,:) = TDOAModel(sensors, truth, nTDOA);
        
        % Store AOA
        AOA(i,:) = AOAModel(sensors, truth, nAOA)';
        
        DeltaT = [TDOA(i,:)'; AOA(i,:)'];

        % Function for lsqnonlin:
        sqfunc = @(e) (DeltaT - h(e))./sqrt(diag(Qv));
        sqfunc2 = @(e) func2outputs(e, sqfunc, J);

        % Estimate using lsqnonlin
        Linear_Meas(i,:)= lsqnonlin(sqfunc2, Linear_Meas(i-1,:)', lb, ub, options);

        % Covariance matrix of measurement
        R_lin(:,:,i) = pinv(J(Linear_Meas(i,:)')) * Qv * pinv(J(Linear_Meas(i,:)'))';

        % Output for debugging:
        %fprintf('\tsolution: [%.0f %.0f %.0f]\n\ttruth:    [%.0f %.0f %.0f]\n', emEst(i,1), emEst(i,2), emEst(i,3), truth(1), truth(2), truth(3))
    end

end

% ------------------------------------------------------------------------------- %
% ------------------------------ Functions -------------------------------------- %
% ------------------------------------------------------------------------------- %
function [val1, val2] = func2outputs(inp,func1,func2)
val1 = func1(inp);
val2 = func2(inp);
end


function Meas = calcMeas(sensors, emitter, nTDOA, nAOA)
% Calculare Measurements (without noise)
    Meas = nan(nTDOA + nAOA,1);
    i = 1;
    % TDOAs
    for n = 2:(numel(sensors))
        if sensors(n).TOA
            % TOA of sensor needs to be available
            Meas(i) = sensors(n).calcTDOAMeas(emitter, sensors(1));
            i = i + 1;
        end
    end
    
    %AOAs
    for n = 1:(numel(sensors))
        if sensors(n).AOA
            % AOA of sensor needs to be available
            [az, el] = sensors(n).calcAOAMeas(emitter);
            Meas(i) = az;
            Meas(i+1) = el;
            i = i + 2;
        end
    end

end

function tdoaVal = TDOAModel(sensors, emitter, nTDOA)
    % TDOA with noise
    tdoaVal = zeros(nTDOA,1);
    for n = 1:(numel(sensors)-1)
        if sensors(n+1).TOA
            % TOA of sensors needs to be available
            tdoa1 = sensors(n+1).calcTDOA(emitter, sensors(1));
            tdoaVal(n) = tdoa1;
        end
    end
end


function AOA_val = AOAModel(sensors, emitter, nAOA)
    % Simulate AOA for all receivers
    AOA_val = zeros(nAOA, 1);
    i = 1;
    for n = 1:2:(numel(sensors)*2)
        if sensors(i).AOA
            [az, el] = sensors(i).calcAOA(emitter);
            AOA_val(n:n+1) = [az, el];
        end
        i = i+1;
    end
end

function J = assembleJacobian(sensors, emitter, nTDOA, nAOA)
% Calculate and asseble Jacobian
%   - calculation is done inside objects
    J = zeros(nTDOA + nAOA, 3);
    i = 1;
    % TDOAs
    for n = 1:(numel(sensors)-1) 
        if sensors(n+1).TOA
            J(i,:) = sensors(n+1).calcTDOAJacobian(emitter, sensors(1));
            i = i + 1;
        end
    end
    % AOAs
    for n = 1:(numel(sensors)) 
        if sensors(n).AOA
            J(i:i+1,:) = sensors(n).calcAOAJacobian(emitter);
            i = i + 2;
        end
    end
end
% ------------------------------------------------------------------------------- %
% ------------------------------------------------------------------------------- %
