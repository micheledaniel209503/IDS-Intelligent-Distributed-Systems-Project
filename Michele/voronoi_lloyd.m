function [L, areas, masses, centroids] = voronoi_lloyd(Robots, Robots_voronoi, X, Y, free_mask, sigma_lineup, sigma_transport)
% DESCRIPTION
%   - Uses "voronoi_labels_grid" to assign each point of a cell to a robot
%   - Builds a pdf (either gaussian or uniform) after reading each robot's working state and target
%   - Computes area, mass and centroid per each cell (each robot)
%
% INPUT
%   Robots(i): class
%       .state  = [x y theta]
%       .working_state  = 'f' | 'l' | 't'
%       .target = [xt yt] 
%   Robots_voronoi: vector of indeces to select robots that will be given a Voronoi cell (maybe optional)
%   X,Y        : meshgrid
%   free_mask  : logic, true = cell is available for tasselation | false = cell is unavailable (obstacle/forbidden area...)
% OUTPUT
%   L          : Voronoi labels
%   areas      : area for each robot
%   masses     : weighted mass of the cell
%   centroids  : centroid position
    
    % Robots_voronoi can be either the indeces or the id's of the robots (since they coincide)
    N = length(Robots_voronoi); % n° robots
    idx_use = Robots_voronoi(:); % indeces of the robots in the list of robots that are going to be given a Voronoi cell
    
    % grid spacing
    dx = mean(diff(unique(X(1,:)))); % horizontal dim of a cell
    dy = mean(diff(unique(Y(:,1)))); % vertical dim of a cell
    cell_area = dx * dy; % area of each grid's cell (elmentary cell)
    
    % positions of the selected robots
    robot_pos = zeros(length(idx_use), 2);
    for k = 1:length(idx_use)
        robot_pos(k,:) = Robots(idx_use(k)).state(1:2);
    end

    % Voronoi tasselation for the selected robots
    L = voronoi_labels_grid(X, Y, robot_pos, free_mask);

    % pdf for each robot(Phi(:,:,i))
    [m,n] = size(X); % size of map
    Phi = zeros(m,n,N); % pdf inizialization (for each robot)
    for k = 1:N % per each robot
        r = Robots(idx_use(k)); % select k-th robot
        working_state = lower(r.working_state);
        switch working_state
            case 'l' % lineup
                mu = pick_target(r, robot_pos(k,:)); % mean value of Phi must be the robot's target --> lineup point
                sg = sigma_lineup;
                Phi(:,:,k) = gauss2d(X, Y, mu, sg^2*eye(2)); % build Phi (bivariate Gaussian)

            case 't' % transport
                mu = pick_target(r, robot_pos(k,:)); % mean value of Phi must be the robot's target --> trajectory defined
                sg = opts.sigma_transport;
                Phi(:,:,k) = gauss2d(X, Y, mu, sg^2*eye(2)); % build Phi (bivariate Gaussian)

            case 'f' % free
                    Phi(:,:,k) = ones(m,n);  % uniform pdf --> optimal coverage (classic Voronoi)
            otherwise % error management
                error("Unrecognized working state '%working_state' for Robots(%d).", working_state, idx_use(k));
        end
    end

    % area, mass, centroid for each Voronoi cell (N_cells = N_robots)
    areas     = zeros(N,1);
    masses    = zeros(N,1);
    centroids = nan(N,2);

    for k = 1:N % for each robot (each cell)
        mask_k = (L == k) & free_mask; % extract grid's points assigned to k-th robot (mask of k-th robot) which are also not occupied by obstacles (free_mask)
        areas(k) = nnz(mask_k) * cell_area; % compute area associated to mask of k-th robot

            if areas(k) == 0 % no area assigned to k-th robot (unlikely)
                masses(k)    = 0;
                centroids(k,:)= [NaN NaN];
                continue;
            end

        w = Phi(:,:,k); % k-th robot's pdf
        w(~mask_k) = 0; % must be zero outside mask --> only integrate it on the mask : this is coherent with theory (only integrate on the voronoi cell's domain)

        M = sum(w(:)) * cell_area; % integral of Phi(x,y)*dxdy  =  sum(Phi)*dA
        masses(k) = M;

            if M <= 1e-30 % !!! this has to be fine tuned !!! if mass is null --> no centroid
                centroids(k,:) = [NaN NaN];
            else % define centroid
                Cx = sum(sum(X .* w)) * cell_area / M; % look at theory
                Cy = sum(sum(Y .* w)) * cell_area / M;
                centroids(k,:) = [Cx, Cy];
            end
    end
end

% PICK TARGET: picks the target stored in rob.target. If not present,
% target becomes 'fallback_pos' (current position in this script).
function mu = pick_target(robot, fallback_pos)
    if all(isfinite(robot.target(:))) % if target is set
        mu = robot.target(:).';
    else                              % if NO target is set
        mu = fallback_pos(:).';
    end
end
