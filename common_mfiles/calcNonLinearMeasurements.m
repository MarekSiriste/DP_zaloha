function [DeltaT, R] = calcNonLinearMeasurements(sensors, true_pos)
    % Compute measurements

    % Inputs:
    %          - sensors: location of each sensor and its characteristics
    %          - true_pos: true positions (trajectory)
    %          - surf: surface information

    % Outputs:
    %          - DeltaT - combination of TDOA and AOA based on availability
    %          - R - covariance matrix of measurements

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
    R = diag([variances_TDOA; variances_AOA].^2); % Covariance matrix of noise 
    
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
    
    % Simulate the TDOA measurement
    TDOA(1,:) = TDOAModel(sensors, truth, nTDOA);
    
    % Simulate AOA
    AOA(1,:) = AOAModel(sensors, truth, nAOA)';
    
    DeltaT(1,:) = [TDOA(1,:)'; AOA(1,:)'];

    % ------------------------------------------------------------------------------- %
    % ------------------------------------------------------------------------------- %

    % ------------------------------------------------------------------------------- %
    % ----------------------- TDOA and AOA calculation ------------------------------ %
    % ------------------------------------------------------------------------------- %


    for i= 2:length(true_pos)
        truth = true_pos(i,1:3)'; % location of emitter

        % Calculate TDOA
        TDOA(i,:) = TDOAModel(sensors, truth, nTDOA)';
        
        % Store AOA
        AOA(i,:) = AOAModel(sensors, truth, nAOA)';
        
        DeltaT(i,:) = [TDOA(i,:)'; AOA(i,:)'];
    end
    

end

% ------------------------------------------------------------------------------- %
% ------------------------------ Functions -------------------------------------- %
% ------------------------------------------------------------------------------- %

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
% ------------------------------------------------------------------------------- %
% ------------------------------------------------------------------------------- %
