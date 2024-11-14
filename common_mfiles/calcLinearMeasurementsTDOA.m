function [Linear_Meas, R_lin] = calcLinearMeasurementsTDOA(sensors, true_pos, T)
    % Compute linear measurements from TDOA only using matlab function
    % helperTDOA2pos - NOT RECOMMENDED 
    % the function does not work correctly - TODO - delete later

    % Inputs:
    %          - sensors: location of each sensor and its characteristics
    %          - true_pos: true positions (trajectory)

    % Outputs:
    %          - Linear_Meas - linear measurements (positions)
    %          - R_lin - linear measurements covariance matrix

    nReceivers = numel(sensors); % Number of receivers

    Linear_Meas = zeros(length(true_pos), 3); 
    R_lin = zeros(3, 3, length(true_pos)); 

    speedOfLight = 299792458; % Speed of light
    
    % Extract the reference sensor's position
    referenceSensor = sensors(1);
    referenceLoc = referenceSensor.position;


    for i = 1:length(true_pos)
        
        truth = true_pos(i,1:3)'; 
        d = zeros(length(sensors)-1, 1);
        delta = zeros(length(sensors)-1, 1);
        S = zeros(length(sensors)-1, 3);
        for j = 2:length(sensors)
            tdoaMeasurement = sensors(j).calcTDOA(truth, referenceSensor); % Simulated TDOA measurement
            receiverLoc = sensors(j).position;

            d(j-1) = tdoaMeasurement * speedOfLight; % TDOA to distance
            delta(j-1) = norm(receiverLoc - referenceLoc)^2 - d(j-1)^2;
            S(j-1, :) = (receiverLoc - referenceLoc)';
        end
        
       % Pseudo-inverse of S
        Swstar = pinv(S);

        % Assemble the quadratic range equation
        STS = (Swstar'*Swstar);
        a = 4 - 4*d'*STS*d;
        b = 4*d'*STS*delta;
        c = -delta'*STS*delta;

        Rs = zeros(2,1);
        % Imaginary solution, return a location outside coverage
        if b^2 < 4*a*c 
            pos = 1e10*ones(3,1);
            R_lin(:,:,i) = 1e10*eye(3);
            return;
        end

        % Two range values
        Rs(1) = (-b + sqrt(b^2 - 4*a*c))/(2*a);
        Rs(2) = (-b - sqrt(b^2 - 4*a*c))/(2*a);

        % If one is negative, use the positive solution
        if prod(Rs) < 0
            Rs = Rs(Rs > 0);
            pos = 1/2*Swstar*(delta - 2*Rs(1)*d) + referenceLoc;
        else % Use range which minimize the error
            xs1 = 1/2*Swstar*(delta - 2*Rs(1)*d);
            xs2 = 1/2*Swstar*(delta - 2*Rs(2)*d);
            e1 = norm(delta - 2*Rs(1)*d - 2*S*xs1);
            e2 = norm(delta - 2*Rs(2)*d - 2*S*xs2);
            if e1  > e2
                pos = xs2 + referenceLoc;
            else
                pos = xs1 + referenceLoc;
            end
        end
        
        % Store the results
        Linear_Meas(i, :) = pos;
        
        % Calculate covariance matrix of the position estimate
        H = zeros(nReceivers-1,3);
        S_matrix = zeros(nReceivers-1,nReceivers-1); % Measurement noise matrix
        for k = 1:nReceivers-1
            e1 = pos - sensors(k+1).position;
            e2 = pos - referenceLoc;
            H(k,:) = (e1'/norm(e1) - e2'/norm(e2)) / T;
            S_matrix(k,k) = sensors(k+1).TOAVariance;
        end

        Pinv = H' / S_matrix * H;
        if rank(Pinv) >= 3
            R_lin(:,:,i) = inv(Pinv);
        else
            R_lin(:,:,i) = eye(3) * 1e10; 
        end
    end
end
