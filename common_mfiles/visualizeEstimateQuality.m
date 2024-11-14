function visualizeEstimateQuality(eval_matrix, estimates, trajectories, surf, quantization, metric, tres)
   % function to visualize the quality of estimates in the map
   % The map is separated into the ,,cubes"
   % Each cube will be colored based on the quality of estimate

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
    
    % Each ,,cube" of eval matrix will have average trace/diff value of Cov matrix
    avgMetric_rmse = nan(nx, ny, nz);
    avgMetric_anees = nan(nx, ny, nz);
    
    for ix = 1:nx
        for iy = 1:ny
            for iz = 1:nz
                if ~isempty(eval_matrix{ix, iy, iz})
                    rmse = [];
                    anees = 0;
                    for idx = 1:length(eval_matrix{ix, iy, iz}.trajectory)
                       % Going throught all trajectories saved in ,,cube"
                        nSim = eval_matrix{ix, iy, iz}.trajectory(idx);
                        pointIdx = eval_matrix{ix, iy, iz}.pointIndex(idx);
                        trajPoint = trajectories{nSim}(pointIdx, 1:3);
                        estPoint = estimates{nSim}.x(pos_index, pointIdx);
                        % Calculate squared difference for RMSE
                        rmse = [rmse, norm(estPoint - trajPoint')^2];
                        
                        % Calculate of ANEES:
                        err = estPoint - trajPoint';
                        anees = anees + (err ./ (estimates{nSim}.P(pos_index, pos_index, pointIdx) * err));
                    end
                        avgMetric_rmse(ix, iy, iz) = sqrt(sum(rmse)/length(rmse));
                        avgMetric_anees(ix, iy, iz) = norm(anees) / (length(anees) * length(err));
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
    figure; 
    hold on; grid on; axis equal;
    xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');
    
    if strcmp(metric, 'rmse')
        % TODO - ANEES
        n = 256;
        greenToYellow = [linspace(0, 1, n/2)', ones(n/2, 1), zeros(n/2, 1)];
        yellowToRed = [ones(n/2, 1), linspace(1, 0, n/2)', zeros(n/2, 1)];
        customCMap = [greenToYellow; yellowToRed];

        colorbar;
        colormap(customCMap);

        hColorbar = colorbar;
        set(hColorbar.Label, 'Interpreter', 'latex');
        hColorbar.Label.String = 'RMSE';
        hColorbar.Limits = [0 1];
        hColorbar.Ticks = [0, 0.5, 1]; 
        lab1 = [num2str(threshold) ' [m]'];
        lab2 = [num2str(2 * threshold) ' [m]'];
        hColorbar.TickLabels = {'0', lab1, lab2};
    end
    
    % Color each cube based on the chosen metric and thresholds
    for ix = 1:nx
        for iy = 1:ny
            for iz = 1:nz
                if ~isnan(avgMetric(ix, iy, iz))
                    if strcmp(metric, 'anees')
                        %ANEES:% TODO - not done yet - just testing
                        if abs(avgMetric(ix, iy, iz)) < 1.2
                            cubeColor = [0, 0, 1]; % Blue = good quality (near 1)
                        elseif avgMetric(ix, iy, iz) > 1.2
                            cubeColor = [1, 0, 0]; % Red = overconfident estimates (>1)
                        else
                            cubeColor = [0, 1, 0]; % Green = underconfident estimates (<1)
                        end
                    else % RMSE:
                        normalizedRMSE = avgMetric(ix, iy, iz) / (2 * threshold);
                        normalizedRMSE = min(max(normalizedRMSE, 0), 1);  % Clamp values to [0, 1]

                        % Map the normalized RMSE to the current colormap
                        colormapIndex = round(normalizedRMSE * (size(customCMap, 1) - 1)) + 1;
                        cubeColor = customCMap(colormapIndex, :);
                    end
                    drawCubes(surf, quantization, ix, iy, iz, cubeColor);
                end
            end
        end
    end
    
    view(3); % Adjust the view to 3D
%     saveas(gcf,'Singer_DWPA_all_scenario4','svg')
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
          'EdgeColor', 'k', 'FaceAlpha', 0.2); % Adjust 'FaceAlpha' for transparency
end