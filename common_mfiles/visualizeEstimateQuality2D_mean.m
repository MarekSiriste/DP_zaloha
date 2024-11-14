function visualizeEstimateQuality2D_mean(eval_matrix, estimates, trajectories, surf, sensors, quantization, tres)
% Visualize quality in 2D
% The visualization of quality is based on RMSE
% averages all Z levels

nx = length(surf.xMin:quantization(1):surf.xMax) - 1;
ny = length(surf.yMin:quantization(2):surf.yMax) - 1;
nz = length(surf.zMin:quantization(3):surf.zMax)-1;



if length(estimates{1}.x(:,1)) == 9
    pos_index = [1 4 7];
else
    pos_index = [1 3 5];
end

threshold = tres;

avgMetric = nan(nx, ny);
for ix = 1:nx
    for iy = 1:ny
        rmse = [];
        for iz = 1:nz
            if ~isempty(eval_matrix{ix, iy, iz})
                for idx = 1:length(eval_matrix{ix, iy, iz}.trajectory)
                    nSim = eval_matrix{ix, iy, iz}.trajectory(idx);
                    pointIdx = eval_matrix{ix, iy, iz}.pointIndex(idx);
                    trajPoint = trajectories{nSim}(pointIdx, 1:3);
                    estPoint = estimates{nSim}.x(pos_index, pointIdx);
                    if isnan(sum(estPoint)) % NAN issue in est
                        continue
                    else
                        rmse = [rmse, norm(estPoint - trajPoint')^2];
                    end
                end
            end
        end
        if ~isempty(rmse)
            avgMetric(ix, iy) = sqrt(sum(rmse) / length(rmse));
        else
            avgMetric(ix, iy) = NaN;  % Assign NaN if rmse is empty (no valid data)
        end
    end
end

figure;
hold on; grid on;
xlabel('x [m]'); ylabel('y [m]');
title(['Map of Estimates Quality Across All Z-Levels']);
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

% saveas(gcf,'fig_SinDWPA_all_z_4','epsc')
end

function color = getColorForValue(value, threshold, customCMap)
    % Normalize the RMSE value to a range of [0, 1]
    normalizedValue = value / (2 * threshold);
    normalizedValue = min(max(normalizedValue, 0), 1);  % <min, max> 

    % Calculate the index for the colormap
    colormapIndex = round(normalizedValue * (size(customCMap, 1) - 1)) + 1;
    color = customCMap(colormapIndex, :);
end



