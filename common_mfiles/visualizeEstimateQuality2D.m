function visualizeEstimateQuality2D(eval_matrix, estimates, trajectories, surf, sensors, quantization, tres, selectedZ)
% Visualize quality in 2D of selected Z level
% The visualization of quality is based on RMSE

nx = length(surf.xMin:quantization(1):surf.xMax) - 1;
ny = length(surf.yMin:quantization(2):surf.yMax) - 1;

zIndices = (surf.zMin:quantization(3):surf.zMax);
[~, selectedZIndex] = min(abs(zIndices - selectedZ));


if length(estimates{1}.x(:,1)) == 9
    pos_index = [1 4 7];
else
    pos_index = [1 3 5];
end

threshold = tres;

% Initialize a matrix for average metric
avgMetric = nan(nx, ny);
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
%title(['Map of Estimates Quality at Z = ', num2str(zIndices(selectedZIndex))]);
colors = {'#000000', '#0072BD', '#7E2F8E', '#EDB120'};

% Plot each sensor as a point in the 3D space
for i = 1:numel(sensors)
    if sensors(i).TOA || sensors(i).TOA
        plot3(sensors(i).x, sensors(i).y, sensors(i).z, "^", ...
            'MarkerSize', 12, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', colors{i}, ...
            'DisplayName', sensors(i).name);
    end
end
legend('Location', 'best');

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



% saveas(gcf,'fig14_2D','epsc')
end

function color = getColorForValue(value, threshold, customCMap)
    % Normalize the RMSE value to a range of [0, 1]
    normalizedValue = value / (2 * threshold);
    normalizedValue = min(max(normalizedValue, 0), 1);  % <min, max> 

    % Calculate the index for the colormap
    colormapIndex = round(normalizedValue * (size(customCMap, 1) - 1)) + 1;
    color = customCMap(colormapIndex, :);
end



