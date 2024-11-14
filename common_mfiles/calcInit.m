function [x_init, P_init] = calcInit(model, LinMeas, R, T)
% Function calculate init of state and covariance matrix for given model
%   - the init is based on first 3 measurements
%   Inputs:
%           - model = current model to initialize
%           - LinMeas = calculated linear measurements (positions)
%           - R = covariance matrices of linear measurements
%           - T = sampling period

%   Outputs:
%           - x_init = init of state
%           - P_init = init of covariance matrix

    % Number of state - for example: Singer has 9 states, NCV only 6
    numStates = length(model.F);
    
    if numStates == 9 % Singer/DWPA
        V = [1 0 0 0 0 0 0 0 0;
             0 0 0 1 0 0 0 0 0;
             0 0 0 0 0 0 1 0 0;
             0 1 0 0 0 0 0 0 0;
             0 0 0 0 1 0 0 0 0;
             0 0 0 0 0 0 0 1 0;
             0 0 1 0 0 0 0 0 0;
             0 0 0 0 0 1 0 0 0;
             0 0 0 0 0 0 0 0 1];
        A = [   0      0          1
             -1/(2*T)  0       1/(2*T)
             1/T^2   -2/T^2    1/T^2];  
        % First three measurements are used
        Z = [LinMeas(1,:) LinMeas(2,:) LinMeas(3,:)]'; 
        tmp_R = blkdiag(R(:,:,1), R(:,:,2),R(:,:,3));
    else % NCV/DWNA
        V = [1 0 0 0 0 0;
             0 0 0 1 0 0;
             0 1 0 0 0 0;
             0 0 0 0 1 0;
             0 0 1 0 0 0;
             0 0 0 0 0 1];
        A = [ 0      1
             -1/T  1/T];  
        % First two measurements are used
        Z = [LinMeas(1,:) LinMeas(2,:)]'; 
        
        tmp_R =  blkdiag(R(:,:,1), R(:,:,2));
    end

    A_3D = blkdiag(A, A, A);
    
    B = A_3D * V;
    
    x_init = B*Z;
    
    P_init = B * tmp_R * B';

end