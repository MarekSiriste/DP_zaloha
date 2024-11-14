function visualizeEstimateQuality_sim(eval_matrix, estimates, trajectories, surf, quantization, metric, tres, simulation)
   % function to visualize the quality of estimates in the map
   % The map is separated into the ,,cubes"
   % Each cube will be colored based on the quality of estimate
   % Draw just one simulation for showcase of the functionality

    % Determine the number of cubes along each dimension
    nx = length(surf.xMin:quantization(1):surf.xMax) - 1;
    ny = length(surf.yMin:quantization(2):surf.yMax) - 1;
    nz = length(surf.zMin:quantization(3):surf.zMax) - 1;
    
    if length(estimates{1}.x(:,1)) == 9 % [position, velocity, acc]
        pos_index = [1 4 7];    % Define to extract only position from vectors
    else
        pos_index = [1 3 5];
    end
    
    threshold = tres;
    
    
    x_r = trajectories{simulation}; % True positions
    est_r = estimates{simulation}.x(pos_index,:)'; % Estimated positions
    
    withinBoundsIdx = x_r(:,1) >= surf.xMin & x_r(:,1) <= surf.xMax & ...
                  x_r(:,2) >= surf.yMin & x_r(:,2) <= surf.xMax & ...
                  x_r(:,3) >= surf.zMin & x_r(:,3) <= surf.zMax;
              
    x_r = x_r(withinBoundsIdx, 1:3);
    est_r = est_r(withinBoundsIdx, :);

    % Start plotting
    figure;
    hold on; grid on; axis equal;
    xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
    title('Map of Estimates Quality for a chosen sample trajectory');
    
    %True pos
    plot3(x_r(:,1), x_r(:,2), x_r(:,3), 'Color', [1, 0, 0], 'LineWidth', 3, 'LineStyle', '-', 'DisplayName', 'True Position');
%     plot3(x_r(:,1), x_r(:,2), x_r(:,3), 'Color', 'black', 'LineWidth', 3, 'LineStyle', '-', 'DisplayName', 'True Position'); % For bad view

    % Estimation
    plot3(est_r(:,1), est_r(:,2), est_r(:,3), 'Color', [0, 0, 1], 'LineWidth', 3, 'LineStyle', '--', 'DisplayName', 'Estimated Position');
    
    legend('show');
    
    
    % Each ,,cube" of eval matrix will have average trace/diff value of Cov matrix
    avgMetric_rmse = nan(nx, ny, nz);
    avgMetric_anees = nan(nx, ny, nz);
    
    err = nan;
    for ix = 1:nx
        for iy = 1:ny
            for iz = 1:nz
                if ~isempty(eval_matrix{ix, iy, iz})
                    rmse = [];
                    anees = 0;
                    for idx = 1:length(eval_matrix{ix, iy, iz}.trajectory)
                       % Going throught all trajectories saved in ,,cube"
                        nSim = eval_matrix{ix, iy, iz}.trajectory(idx);
                        if nSim == simulation
                            pointIdx = eval_matrix{ix, iy, iz}.pointIndex(idx);
                            trajPoint = trajectories{nSim}(pointIdx, 1:3);
                            estPoint = estimates{nSim}.x(pos_index, pointIdx);
                            % Calculate squared difference for RMSE
                            rmse = norm(estPoint - trajPoint')^2;
                            rmse = [rmse, rmse];

                            % Calculate of ANEES:
                            err = estPoint - trajPoint';
                            anees = anees + (err ./ (estimates{nSim}.P(pos_index, pos_index, pointIdx) * err));
                        end
                    end
                       if ~isnan(err)
                        avgMetric_rmse(ix, iy, iz) = sqrt(sum(rmse)/length(rmse));
                        avgMetric_anees(ix, iy, iz) = norm(anees) / (length(anees) * length(err));
                       end
                end
            end
        end
    end

    if strcmp(metric, 'anees')
        avgMetric = log(avgMetric_anees);
    else 
        avgMetric = avgMetric_rmse;
    end
    
    % Visualization

    % Color each cube based on the chosen metric and thresholds
    for ix = 1:nx
        for iy = 1:ny
            for iz = 1:nz
                if ~isnan(avgMetric(ix, iy, iz))
                    if strcmp(metric, 'anees')
                        %ANEES: % TODO - not done yet - just testing
                        if abs(avgMetric(ix, iy, iz)) < 1.2
                            cubeColor = [0, 0, 1]; % Blue = good quality (near 1)
                        elseif avgMetric(ix, iy, iz) > 1.2
                            cubeColor = [1, 0, 0]; % Red = overconfident estimates (>1)
                        else
                            cubeColor = [0, 1, 0]; % Green = underconfident estimates (<1)
                        end
                    else % RMSE:
                            if avgMetric(ix, iy, iz) < 0.1 * threshold
                                cubeColor = [0, 1, 0]; % Green for good
                            elseif avgMetric(ix, iy, iz) < threshold
                                % Interpolate between green and yellow
                                interp = (avgMetric(ix, iy, iz) - 0.1 * threshold) / (0.9 * threshold);
                                cubeColor = [1, 1 - interp, 0]; % Transition from green to yellow
                           elseif avgMetric(ix, iy, iz) < 2 * threshold
                                % Interpolate between yellow and red
                                interp = (avgMetric(ix, iy, iz) - threshold) / threshold;
                                cubeColor = [1, 1 - interp, 0]; % Transition from yellow to red
                            else
                                cubeColor = [1, 0, 0]; % Red for very poor
                            end
                    end
                    drawCubes(surf, quantization, ix, iy, iz, cubeColor);
                end
            end
        end
    end
  
    n = 256;
    greenToYellow = [linspace(0, 1, n/2)', ones(n/2, 1), zeros(n/2, 1)];
    yellowToRed = [ones(n/2, 1), linspace(1, 0, n/2)', zeros(n/2, 1)];
    customCMap = [greenToYellow; yellowToRed];

    colorbar;
    colormap(customCMap);


    hColorbar = colorbar;
    set(hColorbar.Label, 'Interpreter', 'latex');
    hColorbar.Limits = [0 1];
    hColorbar.Ticks = [0, 0.5, 1]; 
    lab1 = [num2str(threshold) ' [m]'];
    lab2 = [num2str(2 * threshold) ' [m]'];
    hColorbar.TickLabels = {'0', lab1, lab2};

    view(3);
    %view(30, 15);
%     saveas(gcf,'fig8_vis_good','svg')
%     
%     view(60, 15);
%     saveas(gcf,'fig9_vis_bad','svg')
end

function drawCubes(surf, quantization, ix, iy, iz, cubeColor)
    % Calculate the base coordinates for the current cube
    baseX = surf.xMin + (ix - 1) * quantization(1);
    baseY = surf.yMin + (iy - 1) * quantization(2);
    baseZ = surf.zMin + (iz - 1) * quantization(3);
    
    % Define cube vertices based on base point
    cubeVertices = [baseX, baseY, baseZ;                                        % 1
                    baseX + quantization(1), baseY, baseZ;                      % 2
                    baseX + quantization(1), baseY + quantization(2), baseZ;    % 3
                    baseX, baseY + quantization(2), baseZ;                      % 4
                    baseX, baseY, baseZ + quantization(3);                      % 5
                    baseX + quantization(1), baseY, baseZ + quantization(3);    % 6
                    baseX + quantization(1), baseY + quantization(2), baseZ + quantization(3); % 7
                    baseX, baseY + quantization(2), baseZ + quantization(3)];   % 8

    %         8--------7
    %        /|       /|
    %       / |      / |
    %     5--------6   |
    %     |  4-----|--3
    %     | /      | /
    %     |/       |/
    %     1--------2

    % Draw the cube with the defined color
    drawCube(cubeVertices, cubeColor);
end

function drawCube(vertices, cubeColor)
    % Draw cube with defined vertices and color
    
    % To clarify each faces of cube:
    %         8--------7
    %        /|       /|
    %       / |      / |
    %     5--------6   |
    %     |  4-----|--3
    %     | /      | /
    %     |/       |/
    %     1--------2
    faces = [1 2 3 4;   % Bottom face
             5 6 7 8;   % Top face
             1 5 8 4;   % Side face
             2 6 7 3;   % Opposite side face
             1 2 6 5;   % Front face
             3 4 8 7];  % Back face

    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', cubeColor, ...
          'EdgeColor', 'k', 'FaceAlpha', 0.1, 'HandleVisibility', 'off'); % Adjust 'FaceAlpha' for transparency
end