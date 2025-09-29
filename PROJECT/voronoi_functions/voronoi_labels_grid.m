function L = voronoi_labels_grid(X, Y, robot_pos, free_mask)
    % VORONOI_LABELS_GRID Function to compute Voronoi labels on a grid
    % DESCRIPTION: takes every cell of the grid and assigns it a number
    % that is either:
    %     - robot(i)'s id : if robot(i) is the closest robot to it
    %     - 0 : if the cell is occupied by something (free_mask determined)
    % Input Arguments:
    %     X - X-coordinates of the grid
    %     Y - Y-coordinates of the grid
    %     robot_pos - positions of the robots
    %     free_mask - binary mask indicating free cells
    % Output Arguments:
    %     L - Voronoi labels for each grid cell

    % Get the dimensions of the grid
    [m,n] = size(X);
    % Number of robots
    n_rob = size(robot_pos,1);
    % Initialize distance array (between each grid point and each robot)
    % with infinity (3D matrix containing distance of rob_i from all
    % points)
    D2 = inf(m,n,n_rob);
    % Calculate squared distances from each robot to each grid point
    for i = 1:n_rob
        dx = X - robot_pos(i,1); % x distance of rob_i from all grid points
        dy = Y - robot_pos(i,2); % y distance of rob_i from all grid points
        D2(:,:,i) = dx.^2 + dy.^2; % overall distance of rob_i from all grid points
    end
    % Assign labels based on the closest robot: index of closest rob_i
    % (which  is also the id) becomes the cell id
    [~, L] = min(D2, [], 3);
    % Set labels to 0 for non-free cells
    L(~free_mask) = 0;
end
