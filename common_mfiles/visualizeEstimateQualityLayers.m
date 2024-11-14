function visualizeEstimateQualityLayers(eval_matrix, estimates, trajectories, surf, quantization, tres, selectedY, selectedZ)
% Visualize quality in 2D of selected Y and Z layer
% The visualization of quality is based on RMSE
% Two layers concept - only for testing

nx = length(surf.xMin:quantization(1):surf.xMax) - 1;
ny = length(surf.yMin:quantization(2):surf.yMax) - 1;
nz = length(surf.zMin:quantization(3):surf.zMax) - 1;

zIndices = (surf.zMin:quantization(3):surf.zMax);
[~, selectedZIndex] = min(abs(zIndices - selectedZ));

yIndices = (surf.yMin:quantization(2):surf.yMax);
[~, selectedYIndex] = min(abs(yIndices - selectedY));


if length(estimates{1}.x(:,1)) == 9
    pos_index = [1 4 7];
else
    pos_index = [1 3 5];
end

threshold = tres;

% Initialize a matrix for average metric
avgMetric = nan(nx, ny);
% Calculation for Z layer
for ix = 1:nx
    for iy = 1:ny
        if ~isempty(eval_matrix{ix, iy, selectedZIndex}) % Use selectedZIndex
            rmse = [];
            for idx = 1:length(eval_matrix{ix, iy, selectedZIndex}.trajectory)
                % Going through all trajectories saved in "cube"
                nSim = eval_matrix{ix, iy, selectedZIndex}.trajectory(idx);
                pointIdx = eval_matrix{ix, iy, selectedZIndex}.pointIndex(idx);
                trajPoint = trajectories{nSim}(pointIdx, 1:3);
                estPoint = estimates{nSim}.x(pos_index, pointIdx);
                % Calculate  RMSE
                rmse = norm(estPoint - trajPoint')^2;
                rmse = [rmse, rmse];
            end
               avgMetric(ix, iy) = sqrt(sum(rmse)/length(rmse));
        end
    end
end

figure;
hold on; grid on;
xlabel('x [m]'); ylabel('y [m]');
title(['Layer Z = ', num2str(zIndices(selectedZIndex))]);
xlim([(surf.xMin - quantization(1)), (surf.xMax + quantization(1))])
ylim([(surf.yMin - quantization(2)), (surf.yMax + quantization(2))])

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

% Color each cell based on RMSE
for ix = 1:nx
    for iy = 1:ny
        if ~isnan(avgMetric(ix, iy))
            color = getColorForValue(avgMetric(ix, iy), threshold, customCMap);
            rectangle('Position', [surf.xMin + (ix-1)*quantization(1), surf.yMin + (iy-1)*quantization(2), quantization(1), quantization(2)], ...
                      'FaceColor', color, 'EdgeColor', 'k'); % Color interpolated
        end
    end
end



% saveas(gcf,'fig15_2D_layers_1','epsc')

avgMetric = nan(nx, nz);
% Calculation for Z layer
for ix = 1:nx
    for iz = 1:nz
        if ~isempty(eval_matrix{ix, selectedYIndex, iz}) % Use selectedZIndex
            rmse = [];
            for idx = 1:length(eval_matrix{ix, selectedYIndex, iz}.trajectory)
                % Going through all trajectories saved in "cube"
                nSim = eval_matrix{ix, selectedYIndex, iz}.trajectory(idx);
                pointIdx = eval_matrix{ix, selectedYIndex, iz}.pointIndex(idx);
                trajPoint = trajectories{nSim}(pointIdx, 1:3);
                estPoint = estimates{nSim}.x(pos_index, pointIdx);
                % Calculate  RMSE
                rmse = norm(estPoint - trajPoint')^2;
                rmse = [rmse, rmse];
            end
               avgMetric(ix, iz) = sqrt(sum(rmse)/length(rmse));
        end
    end
end
figure
hold on; grid on;
xlabel('x [m]'); ylabel('z [m]');
title(['Layer Y = ', num2str(yIndices(selectedYIndex))]);
xlim([(surf.xMin - quantization(1)), (surf.xMax + quantization(1))])
ylim([(surf.yMin - quantization(2)), (surf.yMax + quantization(2))])

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

% Color each cell based on RMSE
for ix = 1:nx
    for iz = 1:nz
        if ~isnan(avgMetric(ix, iz))
            color = getColorForValue(avgMetric(ix, iy), threshold, customCMap);
            % Correctly place and size rectangle
            rectangle('Position', [surf.xMin + (ix-1) * quantization(1), surf.zMin + (iz-1) * quantization(3), quantization(1), quantization(3)], ...
                      'FaceColor', color, 'EdgeColor', 'k'); 
        end
    end
end

xlim([(surf.xMin - quantization(1)), (surf.xMax + quantization(1))])
ylim([(surf.zMin - quantization(3)), (surf.zMax + quantization(3))])

% saveas(gcf,'fig15_2D_layers_2','epsc')
end

function color = getColorForValue(value, threshold, customCMap)
    % Normalize the RMSE value to a range of [0, 1]
    normalizedValue = value / (2 * threshold);
    normalizedValue = min(max(normalizedValue, 0), 1);  % <min, max> 

    % Calculate the index for the colormap
    colormapIndex = round(normalizedValue * (size(customCMap, 1) - 1)) + 1;
    color = customCMap(colormapIndex, :);
end
