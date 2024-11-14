% ------------------------------------------%
% ---------- Configuration file  -----------%
% ------------------------------------------%

function config = ConfigLoader(ncon)
% ncon.... number of config to load

switch(ncon)
    case 1
        config = LoadConfig1();
    case 2
        config = LoadConfig2();
    case 3
        config = LoadConfig3();
    otherwise
        error('Unsupported configuration number: %d', ncon);
end

end

function config = LoadConfig1()
% Surface
surf = Surface();
surf.xMin = -100;
surf.xMax = 100;
surf.yMin = -100;
surf.yMax = 200;
surf.zMin = 0;
surf.zMax = 200;
surf.xQuantization = 20;
surf.yQuantization = 20;
surf.zQuantization = 20;

config.surface = surf;

% Sensors
sensor1 = Sensor();
sensor1.name = 'Reference sensor';
sensor1.x = 0;
sensor1.y = 0;
sensor1.z = 100;
sensor1.AOA = true;
sensor1.AOAMean = 0;
sensor1.AOAVariance = ((pi/180) / 3 )^2;
sensor1.TOA = true;
sensor1.TOAMean = 0;
sensor1.TOAVariance = 2e-10;

sensor2 = Sensor();
sensor2.name = 'Sensor 2';
sensor2.x = 300;
sensor2.y = 300;
sensor2.z = 150;
sensor2.AOA = true;
sensor2.AOAMean = 0;
sensor2.AOAVariance = ((pi/180) / 3 )^2;
sensor2.TOA = true;
sensor2.TOAMean = 0;
sensor2.TOAVariance = 1e-10;

sensor3 = Sensor();
sensor3.name = 'Sensor 3';
sensor3.x = -300;
sensor3.y = 300;
sensor3.z = 50;
sensor3.AOA = true;
sensor3.AOAMean = 0;
sensor3.AOAVariance = ((pi/180) / 3 )^2;
sensor3.TOA = true;
sensor3.TOAMean = 0;
sensor3.TOAVariance = 2e-10;

sensor4 = Sensor();
sensor4.name = 'Sensor 4';
sensor4.x = -300;
sensor4.y = -300;
sensor4.z = 200;
sensor4.AOA = true;
sensor4.AOAMean = 0;
sensor4.AOAVariance = ((pi/180) / 3 )^2;
sensor4.TOA = true;
sensor4.TOAMean = 0;
sensor4.TOAVariance = 3e-10;

config.sensors = [sensor1, sensor2, sensor3, sensor4];

% Choose variant
config.initStates = ConfigInitStates(surf, 2);

end


function config = LoadConfig2()
% Config for testing
% Surface
surf = Surface();
surf.xMin = -100;
surf.xMax = 100;
surf.yMin = -100;
surf.yMax = 100;
surf.zMin = 0;
surf.zMax = 200;
surf.xQuantization = 20;
surf.yQuantization = 20;
surf.zQuantization = 20;

config.surface = surf;

% Sensors
sensor1 = Sensor();
sensor1.name = 'Reference sensor';
sensor1.x = 0;
sensor1.y = 0;
sensor1.z = 100;
sensor1.AOA = true;
sensor1.AOAMean = 0;
sensor1.AOAVariance = ((pi/180) / 3 )^2;
sensor1.TOA = true;
sensor1.TOAMean = 0;
sensor1.TOAVariance = 2e-10;

sensor2 = Sensor();
sensor2.name = 'Sensor 2';
sensor2.x = 200;
sensor2.y = 200;
sensor2.z = 150;
sensor2.AOA = true;
sensor2.AOAMean = 0;
sensor2.AOAVariance = ((pi/180) / 3 )^2;
sensor2.TOA = true;
sensor2.TOAMean = 0;
sensor2.TOAVariance = 1e-10;

sensor3 = Sensor();
sensor3.name = 'Sensor 3';
sensor3.x = -200;
sensor3.y = 200;
sensor3.z = 50;
sensor3.AOA = true;
sensor3.AOAMean = 0;
sensor3.AOAVariance = ((pi/180) / 3 )^2;
sensor3.TOA = true;
sensor3.TOAMean = 0;
sensor3.TOAVariance = 2e-10;

sensor4 = Sensor();
sensor4.name = 'Sensor 4';
sensor4.x = -200;
sensor4.y = -200;
sensor4.z = 200;
sensor4.AOA = true;
sensor4.AOAMean = 0;
sensor4.AOAVariance = ((pi/180) / 3 )^2;
sensor4.TOA = true;
sensor4.TOAMean = 0;
sensor4.TOAVariance = 3e-10;

config.sensors = [sensor1, sensor2, sensor3, sensor4];

% Choose variant
config.initStates = ConfigInitStates(surf, 2);

end

function config = LoadConfig3()
% Config for testing
% Surface
surf = Surface();
surf.xMin = -5000;
surf.xMax = 5000;
surf.yMin = -5000;
surf.yMax = 5000;
surf.zMin = 0;
surf.zMax = 10000;
surf.xQuantization = 500;
surf.yQuantization = 500;
surf.zQuantization = 500;

config.surface = surf;

% Sensors
sensor1 = Sensor();
sensor1.name = 'Reference sensor';
sensor1.x = 0;
sensor1.y = 0;
sensor1.z = 0;
sensor1.AOA = true;
sensor1.AOAMean = 0;
sensor1.AOAVariance = ((pi/180) / 3 )^2;
sensor1.TOA = true;
sensor1.TOAMean = 0;
sensor1.TOAVariance = 2e-10;

sensor2 = Sensor();
sensor2.name = 'Sensor 2';
sensor2.x = -10483.39;
sensor2.y = 15093.59;
sensor2.z = 3.24;
sensor2.AOA = true;
sensor2.AOAMean = 0;
sensor2.AOAVariance = ((pi/180) / 3 )^2;
sensor2.TOA = true;
sensor2.TOAMean = 0;
sensor2.TOAVariance = 1e-10;

sensor3 = Sensor();
sensor3.name = 'Sensor 3';
sensor3.x = 14472.85;
sensor3.y = 7020.92;
sensor3.z = 53.89;
sensor3.AOA = true;
sensor3.AOAMean = 0;
sensor3.AOAVariance = ((pi/180) / 3 )^2;
sensor3.TOA = true;
sensor3.TOAMean = 0;
sensor3.TOAVariance = 2e-10;

sensor4 = Sensor();
sensor4.name = 'Sensor 4';
sensor4.x = -3148.59;
sensor4.y = -18307.34;
sensor4.z = 213.20;
sensor4.AOA = true;
sensor4.AOAMean = 0;
sensor4.AOAVariance = ((pi/180) / 3 )^2;
sensor4.TOA = true;
sensor4.TOAMean = 0;
sensor4.TOAVariance = 3e-10;

config.sensors = [sensor1, sensor2, sensor3, sensor4];

% Choose variant
config.initStates = ConfigInitStates(surf, 3);

end


function init_states = ConfigInitStates(surf, variant)
% Several variants to define init states


switch (variant)
    case 1
        % Variant 1 - just 4 sides of cube
        %           - each init points into the opposite side
        
        % To not include influence of init condition in evaluation
        extend_surf = 2; % Extend surface little bit
        
        % The range of each axis:
        xRange = (surf.xMin - extend_surf):surf.xQuantization:(surf.xMax + extend_surf);
        yRange = (surf.yMin - extend_surf):surf.yQuantization:(surf.yMax + extend_surf);
        zRange = surf.zMin:surf.zQuantization:surf.zMax;
        
        % First face of the ,,cube"
        xMinFace = surf.xMin;
        [xGrid1, yGrid1, zGrid1] = meshgrid(xMinFace, yRange, zRange);
        
        % Second face of the ,,cube"
        yMinFace = surf.yMin;
        [xGrid2, yGrid2, zGrid2] = meshgrid(xRange, yMinFace, zRange);
        
        % Third face of the ,,cube"
        xMaxFace = surf.xMax;
        [xGrid3, yGrid3, zGrid3] = meshgrid(xMaxFace, yRange, zRange);
        
        % Fourth face of the ,,cube"
        yMaxFace = surf.yMax;
        [xGrid4, yGrid4, zGrid4] = meshgrid(xRange, yMaxFace, zRange);
        
        %     % fifth face, fixed at zMax
        %     zMaxFace = surf.zMax;
        %     [xGrid5, yGrid5, zGrid5] = meshgrid(xRange, yRange, zMaxFace);
        
        % Define init positions (two faces of cube)
        initPositions = [xGrid1(:), yGrid1(:), zGrid1(:);
            xGrid2(:), yGrid2(:), zGrid2(:);
            xGrid3(:), yGrid3(:), zGrid3(:);
            xGrid4(:), yGrid4(:), zGrid4(:)];
        
        
        % Center of the cube
        centerX = (surf.xMax + surf.xMin) / 2;
        centerY = (surf.yMax + surf.yMin) / 2;
        
        
        % velocities toward the opposite edges of the x and y dimensions
        oppositeX = centerX + (centerX - initPositions(:, 1)); % Opposite of current x position
        oppositeY = centerY + (centerY - initPositions(:, 2)); % Opposite of current y position
        
        % Vector points to the opposite sides (x and y) with zero velocity in z axis
        vectors = [oppositeX - initPositions(:,1), oppositeY - initPositions(:,2), zeros(size(initPositions, 1), 1)];
        
        %     % Vectors points to the center of cube (but with 0 velocity in z axis)
        %     vectors = [centerX - initPositions(:,1), centerY - initPositions(:,2), zeros(size(initPositions, 1), 1)];
        
        % normalizing
        norms = sqrt(sum(vectors.^2, 2));
        
        % Replace zeros in norms with 1 to avoid division by zero
        norms(norms == 0) = 1;
        normalizedVectors = vectors ./ norms;
        
        initVelocities = normalizedVectors;
        
        % Assume acc as 0 for all cases
        initAccelerations = zeros(size(initPositions));
        
    case 2
        % Variant 2 - perimeter
        %           - each init points in the center of cube
        % The range of each axis:
        xRange = surf.xMin:surf.xQuantization:surf.xMax;
        yRange = surf.yMin:surf.yQuantization:surf.yMax;
        zRange = surf.zMin:surf.zQuantization:surf.zMax;
        [xGrid, yGrid, zGrid] = meshgrid(xRange, yRange, zRange);
        all_pos = [xGrid(:), yGrid(:), zGrid(:)];
        isPerimeter = (xGrid == surf.xMin | xGrid == surf.xMax | ...
            yGrid == surf.yMin | yGrid == surf.yMax | ...
            zGrid == surf.zMin | zGrid == surf.zMax);
        perimeterPositions = all_pos(isPerimeter, :);
        %{
        scatter3(perimeterPositions(:,1), perimeterPositions(:,2), perimeterPositions(:,3), 'filled');
        axis equal;
        xlabel('X'); ylabel('Y'); zlabel('Z');
        title('Initial Positions within Cube');
        %}
        initPositions = perimeterPositions;
        
        % Center of the cube
        centerX = (surf.xMax + surf.xMin) / 2;
        centerY = (surf.yMax + surf.yMin) / 2;
        centerZ = (surf.zMax + surf.zMin) / 2;
        
        % Vectors points to the center of cube
        vectors = [centerX - initPositions(:,1), centerY - initPositions(:,2), centerZ - initPositions(:,3)];
        
        % normalizing
        norms = sqrt(sum(vectors.^2, 2));
        
        % Replace zeros in norms with 1 to avoid division by zero
        norms(norms == 0) = 1;
        normalizedVectors = vectors ./ norms;
        
        initVelocities = normalizedVectors;
        
        % Assume acc as 0 for all cases
        initAccelerations = zeros(size(initPositions));
        
    case 3
        % Variant 3
        %       - full grid
        % Configuration of all init states
        [xGrid, yGrid, zGrid] = meshgrid(surf.xMin:(surf.xQuantization):surf.xMax, ...
            surf.yMin:(surf.yQuantization):surf.yMax, ...
            surf.zMin:(surf.zQuantization):surf.zMax);
        initPositions = [xGrid(:), yGrid(:), zGrid(:)];
        
        % Center of the cube
        centerX = (surf.xMax + surf.xMin) / 2;
        centerY = (surf.yMax + surf.yMin) / 2;
        centerZ = (surf.zMax + surf.zMin) / 2;
        
        % Vectors points to the center of cube
        vectors = [centerX - initPositions(:,1), centerY - initPositions(:,2), centerZ - initPositions(:,3)];
        
        % normalizing
        norms = sqrt(sum(vectors.^2, 2));
        
        % Replace zeros in norms with 1 to avoid division by zero
        norms(norms == 0) = 1;
        normalizedVectors = vectors ./ norms;
        
        initVelocities = normalizedVectors;
        
        % Assume acc as 0 for all cases
        initAccelerations = zeros(size(initPositions));
    otherwise
        disp('No variant was chosen - full grid will be calculated')
        %       - full grid - even more detailed
        % Configuration of all init states
        [xGrid, yGrid, zGrid] = meshgrid(surf.xMin:(surf.xQuantization/2):surf.xMax, ...
            surf.yMin:(surf.yQuantization/2):surf.yMax, ...
            surf.zMin:(surf.zQuantization/2):surf.zMax);
        initPositions = [xGrid(:), yGrid(:), zGrid(:)];
        
        % Center of the cube
        centerX = (surf.xMax + surf.xMin) / 2;
        centerY = (surf.yMax + surf.yMin) / 2;
        centerZ = (surf.zMax + surf.zMin) / 2;
        
        % Vectors points to the center of cube
        vectors = [centerX - initPositions(:,1), centerY - initPositions(:,2), centerZ - initPositions(:,3)];
        
        % normalizing
        norms = sqrt(sum(vectors.^2, 2));
        
        % Replace zeros in norms with 1 to avoid division by zero
        norms(norms == 0) = 1;
        normalizedVectors = vectors ./ norms;
        
        initVelocities = normalizedVectors;
        
        % Assume acc as 0 for all cases
        initAccelerations = zeros(size(initPositions));
end


init_states = [initPositions(:, 1), initVelocities(:, 1), initAccelerations(:, 1), ...
    initPositions(:, 2), initVelocities(:, 2), initAccelerations(:, 2), ...
    initPositions(:, 3), initVelocities(:, 3), initAccelerations(:, 3)];

end
