function [x_kk,P_kk,modeProb,x_kk_j,P_kk_j] = IMM(model_prob,pij, x_prev, P_prev, models, sensors, Z, R, Meas_type, numstep,MeasInterval)
% IMM algorthm
% Inputs:
%           - model_prob - likelihood
%           - pij - transition probability
%           - x_prev - previous state
%           - P_prev - previous cov
%           - models - array of struct models
%           - sensors - receivers 
%           - Z - linear/nonlinear measurement (based on MeasType)
%           - R - linear/nonlinear measurement covariance matrix (based on
%           Meas_type)
%           - Meas_type - Linear/Nonlinear - choose between approaches
%           - numstep - actual step of filter 
%           - MeasInterval - interval between measurements that the filters consider

% Outputs:
%           - x_kk - \hat{x}(k|k) - estimate of state
%           - P_kk - P(k|k) - estimate of covariance
%           - modeProb - \mu_j(k) - Mode probability
%           - x_kk_j - state estimates of r models
%           - P_kk_j - covariance estimate of r models

    [r, ~] = size(pij); % Number of models, r should be 2 in your case

    %------------------------------------------------%
    % Step 1) Calculation of mixing probabilities
    %------------------------------------------------%
    c_j = pij * model_prob';

    % Mixing probability for each model to each model
    mu_ij =  (pij .* (model_prob' ./ c_j'))';
    
    %------------------------------------------------%
    % Step 2) Mixing 
    %------------------------------------------------%
    x_mixed = zeros(length(x_prev), r); % Initialize mixed states for each model
    P_mixed = zeros(size(P_prev, 1), size(P_prev, 2), r); % Initialize mixed covariances for each model
    
    for j = 1:r
        mu_ij_j = mu_ij(j, :)';
        x_mixed(:, j) = sum(x_prev .* mu_ij_j', 2);
        for i = 1:r
            P_mixed(:, :, j) = P_mixed(:, :, j) + mu_ij_j(i) * (P_prev(:, :, i) + (x_prev(:, i) - x_mixed(:, j)) * (x_prev(:, i) - x_mixed(:, j))');
        end
    end

    %------------------------------------------------%
    % Step 3) Mode matched filtering (Kalman Filters for each model)
    %------------------------------------------------%
    x_kk_j = zeros(length(x_prev), r);
    P_kk_j = zeros(size(P_prev, 1), size(P_prev, 2), r);
    Lambda = zeros(1, r);
    
    for i = 1:r
        % Kalman prediction step for each model
        [x_pred, P_pred] = Kalman_prediction(x_mixed(:, i), P_mixed(:, :, i), models{i}.F, models{i}.Q);
        
        if strcmp(Meas_type, 'Linear') || strcmp(Meas_type, 'linear') % Linear measuremnts available
            % Linear measurements
            % Kalman filtering step for each model
            if mod(numstep, MeasInterval) == 0 % Only MeasInterval-th measurement is performed
                [x_filt, P_filt, innovation, S] = Kalman_filter(x_pred, P_pred, Z, models{i}.H, R);
                in = 1;
            else
                % Skip the update - just prediction
                x_filt = x_pred;
                P_filt = P_pred;
                in = nan; % Innovation is not available
            end
        else
            if mod(numstep, MeasInterval) == 0 % Only MeasInterval-th measurement is performed
                [x_filt, P_filt, innovation, S] = Extended_Kalman_filter(x_pred, P_pred, models{i}.H, sensors, Z, R);
                in = 1;
            else
                % Skip the update - just prediction
                x_filt = x_pred;
                P_filt = P_pred;
                in = nan; % Innovation is not available
            end
        end
        
        % Update the state estimates, covariances, and calculate likelihood
        x_kk_j(:, i) = x_filt;
        P_kk_j(:, :, i) = P_filt;
        % Calculate likelihood if innov and S are valid
        if ~isnan(in)
            Lambda(i) = LAMBDA(innovation', S); % Calculate likelihood of each model based on the innovation
        end
    end

    %------------------------------------------------%
    % Step 4) Mode Probability Update
    %------------------------------------------------%
    if sum(Lambda .* c_j') > 0
        modeProb = (Lambda .* c_j') / sum(Lambda .* c_j');
    else
        modeProb = c_j';
    end
%     modeProb = (Lambda .* c_j') / sum(Lambda .* c_j');

    %------------------------------------------------%
    % Step 5) Estimate and Covariance Combination
    %------------------------------------------------%
    x_kk = zeros(length(x_prev), 1);
    P_kk = zeros(size(P_prev, 1), size(P_prev, 2));
    
    for i = 1:r
        x_kk = x_kk + x_kk_j(:, i) * modeProb(i);
        P_kk = P_kk + (modeProb(i) * (P_kk_j(:, :, i) + (x_kk_j(:, i) - x_kk) * (x_kk_j(:, i) - x_kk)'));
    end
end

%----------------------------------------------------------------------%
%-----------------------Kalman filter----------------------------------%
%----------------------------------------------------------------------%
function[xk,Pk,inov,S] = Kalman_filter(xkk_1, Pkk_1, Zk, H, R)
% Filtering step of Kalman Filter
% Outputs:
%           - xkk_1 - prediction state
%           - Pkk_1 - prediction covariance matrix
%           - xk - x(k|k) actual state
%           - Pk - P(k|k) actual covariance
%           - inov - innovation measurement 
%           - S -innovation cov

    % Filtering step:
    inov = (Zk - (H*xkk_1)')'; 
    S = H*Pkk_1*H' + R; 
    K = Pkk_1*H'*inv(S); % Kalman gain - TODO - pinv

    xk  = xkk_1 + K*inov; 

    I = eye(size(Pkk_1));
    Pk = (I - K*H)*Pkk_1*(I - K*H)' + K*R*K'; % Joseph form update
end

function[xkk_1,Pkk_1] = Kalman_prediction(xk_1,Pk_1,F,Q)
% Prediction step of Kalman Filter
% Outputs:
%           - xk_1 - x(k-1|k-1) state
%           - Pk_1 - P(k-1|k-1) covariance


    % Predition step:
    Pkk_1 = F*Pk_1*F' + Q; 
    xkk_1 = F*xk_1;

end

%----------------------------------------------------------------------%
%----------------------- Extended Kalman filter------------------------%
%----------------------------------------------------------------------%

function[xk, Pk, innov, S] = Extended_Kalman_filter(xkk_1, Pkk_1, H_m, sensors, DeltaT, R)
    % Inputs:
    %  - xkk_1: The prior state estimate at time k-1
    %  - Pkk_1: The prior covariance matrix at time k-1
    %  - sensors: Array of sensor objects
    %  - TDOA: Time Difference of Arrival measurements
    %  - AOA: Angle of Arrival measurements

    % To extract positions from vector of each model
    if length(xkk_1) > 6 
        pos = [1, 4, 7];
    else
        pos = [1, 3, 5];
    end
    
    xkk_1 = real(xkk_1); % TODO - problem with AOA !

    % Initialization
    H = zeros(length(DeltaT), numel(xkk_1(pos)));
    Z_pred = zeros(length(DeltaT), 1);

    % Calculate TDOA Jacobian for all available TDOA measurements
    index = 0;
    for i = 2:numel(sensors)
        if sensors(i).TOA
            H(index + 1, :) = sensors(i).calcTDOAJacobian(xkk_1(pos), sensors(1));
            Z_pred(index + 1) = sensors(i).calcTDOA(xkk_1(pos), sensors(1));
            index = index + 1;
        end
    end
    nTdoa = index; % number of available TDOA

    % Calculate Jacobians for all availables AOA measurements
    for i = 1:numel(sensors)
        if sensors(i).AOA
            AOA_Jacobian = sensors(i).calcAOAJacobian(xkk_1(pos));
            [az, el] = sensors(i).calcAOAMeas(xkk_1(pos));
            H(index + (1:2), :) = AOA_Jacobian;
            Z_pred(index + 1) = az;
            Z_pred(index + 2) = el;
            index = index + 2;
        end
    end
    
    Z = DeltaT'; % Combine TDOAs and AOAs measurements
    % innovation
    innov = Z - Z_pred;
    
    % If AOA are available
    if length(innov) ~= nTdoa
        % Wrap
        innov(nTdoa+1:end) = mod(innov(nTdoa+1:end)+pi,2*pi)-pi;
    end
    
    % Modify H matrix only for positions
    H = H * H_m;
    
    % Calculate the innovation covariance
    S = H * Pkk_1 * H' + R;

    % Calculate the Kalman gain
    K = Pkk_1 * H' * inv(S);
    
    % Update the state estimate and covariance matrix
    xk = xkk_1 + K * innov;
    Pk = (eye(size(Pkk_1)) - K * H) * Pkk_1 * (eye(size(Pkk_1)) - K * H)' + K*R*K'; % Joseph form
    
end


%----------------------------------------------------------------------%
%-----------------------Likelihood-------------------------------------%
%----------------------------------------------------------------------%
function y = LAMBDA(meas_in,S) 
% Computation of likelihood - gaussian pdf
    E = -1/2.*(meas_in*inv(S)*meas_in'); % 
    y = exp(E) / sqrt(2*pi*det(S));
end
%----------------------------------------------------------------------%
