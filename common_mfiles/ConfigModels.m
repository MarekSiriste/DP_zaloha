function [NCV, DWNA, DWPA, SINGER,initialProbabilities, transitionMatrix, T] = ConfigModels()
    % Function loads configuration of the models
    % Sampling period
    T = 0.03;
    
    
    % Nearly constant velocity model
    NCV.description = 'Nearly constant velocity model';
    % 1D:
    NCV.F_alpha = [1, T; 0, 1];
    NCV.q = 100;
    NCV.Q_alpha = NCV.q * [T^3/3, T^2/2; T^2/2, T];
    NCV.mu_alpha = [0, 0];
    % Define it for 3D:
    NCV.F = blkdiag(NCV.F_alpha, NCV.F_alpha, NCV.F_alpha);
    NCV.H = [1 0 0 0 0 0
             0 0 1 0 0 0
             0 0 0 0 1 0]; % Consider only position as a measurement
    NCV.Q =  blkdiag(NCV.Q_alpha, NCV.Q_alpha, NCV.Q_alpha);
    NCV.mu = [NCV.mu_alpha, NCV.mu_alpha, NCV.mu_alpha];
    % Initial condition distribution
    % 1D:
    NCV.mu_x0_alpha = [50 10]';
    NCV.Q_x0_alpha = [1 1
                     1 4];
    % 3D:
    NCV.mu_x0 = [NCV.mu_x0_alpha; NCV.mu_x0_alpha; NCV.mu_x0_alpha];
    NCV.Q_x0 = blkdiag(NCV.Q_x0_alpha, NCV.Q_x0_alpha, NCV.Q_x0_alpha);

    % Discrete white noise acceleration model
    DWNA.description = 'Discrete white noise acceleration model';
    % 1D:
    DWNA.F_alpha = [1, T; 0, 1];
    DWNA.q = 30;
    DWNA.Q_alpha = DWNA.q * [T^4/4, T^3/2; T^3/2, T^2];
    DWNA.mu_alpha = [0, 0];
    % 3D:
    DWNA.F = blkdiag(DWNA.F_alpha, DWNA.F_alpha, DWNA.F_alpha);
    DWNA.H = [1 0 0 0 0 0
             0 0 1 0 0 0
             0 0 0 0 1 0];
    DWNA.Q = blkdiag(DWNA.Q_alpha, DWNA.Q_alpha, DWNA.Q_alpha);
    DWNA.mu = [DWNA.mu_alpha, DWNA.mu_alpha, DWNA.mu_alpha];
    % Initial condition distribution
    % 1D:
    DWNA.mu_x0_alpha = [50 10]';
    DWNA.Q_x0_alpha = [1 1
                1 4];
    % 3D:
    DWNA.mu_x0 = [DWNA.mu_x0_alpha; DWNA.mu_x0_alpha; DWNA.mu_x0_alpha];
    DWNA.Q_x0 = blkdiag(DWNA.Q_x0_alpha, DWNA.Q_x0_alpha, DWNA.Q_x0_alpha);

    % Discrete Wiener process acceleration model
    DWPA.description = 'Discrete Wiener process acceleration model';
    % 1D:
    DWPA.F_alpha = [1, T, T^2/2; 0, 1, T; 0, 0, 1];
    DWPA.q = 30;
    DWPA.Q_alpha =  DWPA.q * [T^4/4, T^3/2, T^2/2;
                              T^3/2, T^2/2, T;
                              T^2/2, T,     1]; 
    DWPA.mu_alpha = [0, 0, 0];
    
    % 3D:
    DWPA.F = blkdiag(DWPA.F_alpha, DWPA.F_alpha, DWPA.F_alpha);
    DWPA.Q = blkdiag(DWPA.Q_alpha, DWPA.Q_alpha, DWPA.Q_alpha);
    DWPA.mu = [DWPA.mu_alpha, DWPA.mu_alpha, DWPA.mu_alpha];
    DWPA.H = [1 0 0 0 0 0 0 0 0
              0 0 0 1 0 0 0 0 0
              0 0 0 0 0 0 1 0 0];
    
    % Initial condition distribution
    DWPA.mu_x0_alpha = [10 0 0]';
    DWPA.Q_x0_alpha = [1 1 0
                 1 4 0
                 0 0 4];
    DWPA.mu_x0 = [DWPA.mu_x0_alpha; DWPA.mu_x0_alpha; DWPA.mu_x0_alpha];
    DWPA.Q_x0 = blkdiag(DWPA.Q_x0_alpha, DWPA.Q_x0_alpha, DWPA.Q_x0_alpha);
    
    % Singer acceleration model
    SINGER.description = 'Singer acceleration model';
    tau_m = 10;%0.5 ; 
    sigma_m = 0.1;
    beta = 1/tau_m;
    rho_m = exp(-beta * T);
    % 1D:
    SINGER.F_alpha = [1, T, (1/beta^2) * (-1 + beta * T + rho_m); ...
                      0, 1, (1/beta) * (1 - rho_m); ...
                      0, 0, rho_m];
    q11 = (1/(2*beta^5)) * (1 - exp(-2*beta*T) + 2*beta*T + (2*beta^3*T^3)/3 - 2*beta^2*T^2 - 4*beta*T*exp(-beta*T));
    q12 = (1/(2*beta^4)) * (exp(-2*beta*T) + 1 - 2*exp(-beta*T) + 2*beta*T*exp(-beta*T) - 2*beta*T + beta^2*T^2);
    q13 = (1/(2*beta^3)) * (1 - exp(-2*beta*T) - 2*beta*T*exp(-beta*T));
    q22 = (1/(2*beta^3)) * (4*exp(-beta*T) - exp(-2*beta*T) - 3 + 2*beta*T);
    q23 = (1/(2*beta^2)) * (exp(-2*beta*T) + 1 - 2*exp(-beta*T));
    q33 = (1/(2*beta)) * (1 - exp(-2*beta*T));
    SINGER.Q_alpha = (2*sigma_m^2/tau_m) * [q11, q12, q13; ...
                                            q12, q22, q23; ...
                                            q13, q23, q33];
    SINGER.mu_alpha = [0, 0, 0];

    % 3D:
    SINGER.F = blkdiag(SINGER.F_alpha, SINGER.F_alpha, SINGER.F_alpha);
    SINGER.Q = blkdiag(SINGER.Q_alpha, SINGER.Q_alpha, SINGER.Q_alpha);
    SINGER.mu = [SINGER.mu_alpha, SINGER.mu_alpha, SINGER.mu_alpha];
    SINGER.H = [1 0 0 0 0 0 0 0 0
                0 0 0 1 0 0 0 0 0
                0 0 0 0 0 0 1 0 0];
            
    % Initial condition distribution
    % 1D:
    SINGER.mu_x0_alpha = [10 0 0]';
    SINGER.Q_x0_alpha = [1 1 0
                         1 4 0
                         0 0 4];
    % 3D:
    SINGER.mu_x0 = [SINGER.mu_x0_alpha; SINGER.mu_x0_alpha; SINGER.mu_x0_alpha];
    SINGER.Q_x0 = blkdiag(SINGER.Q_x0_alpha, SINGER.Q_x0_alpha, SINGER.Q_x0_alpha);
    
    
    initialProbabilities = [0.3, 0.7];  % Equal probability to start in each state
    transitionMatrix = [0.95, 0.05;       % Transition probabilities
                        0.05, 0.95];

end
