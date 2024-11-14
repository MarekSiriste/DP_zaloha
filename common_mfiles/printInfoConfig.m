function printInfoConfig(config)
    % Print info about configuration
    disp('Surface Information:');
    disp(['x-axis area: [', num2str(config.surface.xMin), ', ', num2str(config.surface.xMax) ']',' quantization: ', num2str(config.surface.xQuantization)])
    disp(['y-axis area: [', num2str(config.surface.yMin), ', ', num2str(config.surface.yMax) ']',' quantization: ', num2str(config.surface.yQuantization)])
    disp(['z-axis area: [', num2str(config.surface.zMin), ', ', num2str(config.surface.zMax) ']',' quantization: ', num2str(config.surface.zQuantization)])
    disp('----------------------');
    disp('Sensor Information:');
    for i = 1:numel(config.sensors)
        disp(['Sensor Name: ', config.sensors(i).name]);
        disp(['Position: (', num2str(config.sensors(i).x), ', ', num2str(config.sensors(i).y), ', ', num2str(config.sensors(i).z), ')']);
        disp(['AOA Available: ', num2str(config.sensors(i).AOA)]);
        disp(['AOA Mean: ', num2str(config.sensors(i).AOAMean)]);
        disp(['AOA Variance: ', num2str(config.sensors(i).AOAVariance)]);
        disp(['TOA available: ', num2str(config.sensors(i).TOA)]);
        disp(['TOA Mean: ', num2str(config.sensors(i).TOAMean)]);
        disp(['TOA Variance: ', num2str(config.sensors(i).TOAVariance)]);
        disp('----------------------');
    end
    
    % Figure the sensors
    figure;
    hold on;
    grid on;
    colors = {'#000000', '#0072BD', '#7E2F8E', '#EDB120'};
    % Plot each sensor as a point in the 3D space
    for i = 1:numel(config.sensors)
         plot3(config.sensors(i).x, config.sensors(i).y, config.sensors(i).z, "^", ...
            'MarkerSize', 10, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', colors{i}, ... % Random color for sensor
            'DisplayName', config.sensors(i).name);
    end

    % Setting the axis limits based on surfaceLimits
    %   - not recommended for receivers out of area of surface
    % axis([config.surface.xMin config.surface.xMax config.surface.yMin config.surface.yMax config.surface.zMin config.surface.zMax]);

    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    legend('Location', 'bestoutside');
    title('Map of sensors');
    %saveas(gcf,'fig','epsc')

    view(0,90)  % XY
    % pause(2)
    % view(135, 30);

end

