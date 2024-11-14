function eval_matrix = CreateEvalMatrix(surf, quantization, trajectories)
% Creates evaluation matrix
% Split surface into the ,,cubes" of dimensions given by the quantization
% The function goes throught all trajectories and map the points inside
% the cubes

% For example:
%       - surf = [xMin = -10, xMax = 10, yMin = -10,yMax = 10, zMin = -10, zMax = 10]
%       - quant = [1, 1, 1]
%       - It creates eval_Matrix of dimension 21x21x21 (cubes are 1x1x1)
%       - lets consider point x = [0.5, 3, -9.5]
%       - it will be saved as a trajectory to the {11, 14, 1}

% Define matrix of evaluation dimension:
nx = length(surf.xMin:quantization(1):surf.xMax);
ny = length(surf.yMin:quantization(2):surf.yMax);
nz = length(surf.zMin:quantization(3):surf.zMax);

eval_matrix = cell(nx, ny, nz);

% Define boundaries of the map
xBounds = [surf.xMin, surf.xMax];
yBounds = [surf.yMin, surf.yMax];
zBounds = [surf.zMin, surf.zMax];


% Define anonymous function for calculating indices
getIndicesForPoint = @(point) getEvalMatrixIndices(point, xBounds, yBounds, zBounds, quantization);

%indeces1 = getEvalMatrixIndices(x_r, [surf.xMin surf.xMax], [surf.yMin surf.yMax], [surf.zMin surf.zMax], quantization);

for nSim = 1:length(trajectories)
    x_r = trajectories{nSim}; % Current trajectory
    
    % For each point in the trajectory
    for i = 1:size(x_r, 1)
        point = x_r(i, 1:3); % Extract point - only position
        indices = getIndicesForPoint(point); % Calculate indexes of matrix
        
        % Is point in the range of defined map/surface ?
        if all(indices >= 1 & [indices(1) <= nx, indices(2) <= ny, indices(3) <= nz])
            % Create structure of trajectory and index of actual point
            % which lies in the one of the quantized ,,cubes"
            if isempty(eval_matrix{indices(1), indices(2), indices(3)})
                eval_matrix{indices(1), indices(2), indices(3)} = struct('trajectory', [], 'pointIndex', []);
            end
            eval_matrix{indices(1), indices(2), indices(3)}.trajectory = [eval_matrix{indices(1), indices(2), indices(3)}.trajectory; nSim];
            eval_matrix{indices(1), indices(2), indices(3)}.pointIndex = [eval_matrix{indices(1), indices(2), indices(3)}.pointIndex; i];
        end
    end
end
end

function indices = getEvalMatrixIndices(x_r, x_b, y_b, z_b, quant)
% Function for calculation of position inside evaluation matrix    

    % For each dimension
    x_pos = x_r(:, 1); % For x
    y_pos = x_r(:, 2); % For y
    z_pos = x_r(:, 3); % For z
    
    % Calculate indices for x, y, z
    indices_x = arrayfun(@(x) floor((x - x_b(1)) ./ quant(1)) + 1, x_pos);
    indices_y = arrayfun(@(y) floor((y - y_b(1)) ./ quant(2)) + 1, y_pos);
    indices_z = arrayfun(@(z) floor((z - z_b(1)) ./ quant(3)) + 1, z_pos);
    
    % Combine indices
    indices = [indices_x, indices_y, indices_z];
end

