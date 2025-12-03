function [L, areas, masses, centroids, Phi, Wrs_set] = voronoi_lloyd_ring_dyna_unc_reactive(Robots, Robots_voronoi, X, Y, free_mask, sigma_point, sigma_ring, R, use_reactive)
% DESCRIPTION
%   - Uses "voronoi_labels_grid" to assign each point of a cell to a robot
%   - Builds a pdf (either gaussian or uniform) after reading each robot's working state and target
%   - Computes area, mass and centroid per each cell (each robot)
%
% INPUT
%   Robots(i): class
%       .state  = [x y theta]
%       .working_state  = 'p' | 'f' | 'l' | 'r'
%       .target = [xt yt] 
%   Robots_voronoi: vector of indeces to select robots that will be given a Voronoi cell (maybe optional)
%   X,Y        : meshgrid
%   free_mask  : logic, true = cell is available for tasselation | false = cell is unavailable (obstacle/forbidden area...)
%   sigma_ring : value for the ring pdf
%   R          : ring radius
%   use_reactive: flag --> if false --> use standard voronoi
% OUTPUT
%   L          : Voronoi labels
%   areas      : area for each robot
%   masses     : weighted mass of the cell
%   centroids  : centroid position
    
    % Robots_voronoi can be either the indeces or the id's of the robots (since they coincide)
    N = length(Robots_voronoi); % n° robots
    idx_use = Robots_voronoi(:); % indeces of the robots in the list of robots that are going to be given a Voronoi cell
    obs = ~free_mask; % obstacles mask
    Wrs_set = cell(N,1); % initialization
    
    % grid spacing
    dx = mean(diff(unique(X(1,:)))); % horizontal dim of a cell
    dy = mean(diff(unique(Y(:,1)))); % vertical dim of a cell
    cell_area = dx * dy; % area of each grid's cell (elmentary cell)
    
    % positions of the selected robots
    robot_pos_est = zeros(length(idx_use), 2);
    for k = 1:length(idx_use)
        robot_pos_est(k,:) = Robots(idx_use(k)).state_est(1:2);
    end

    % Voronoi tasselation for the selected robots
    L = voronoi_labels_grid(X, Y, robot_pos_est, free_mask);

    % Reactive control for static obstacles (optional)
    if use_reactive
        ds = 0.5*min( mean(diff(unique(X(1,:)))), mean(diff(unique(Y(:,1)))) ); % sampling step for the line of sight check
        for k = 1:N % for each robot
            pr = robot_pos_est(k,:);
            Vmask = (L == k) & free_mask; % voronoi cell
            Wrs_set{k} = build_Wrs_fast_ray2(pr, Vmask, obs, X, Y, Robots(k).sr, ds); % build the reduced cell
        end
    end


    % pdf for each robot(Phi(:,:,i))
    [m,n] = size(X); % size of map
    Phi = zeros(m,n,N); % pdf inizialization (for each robot)
    for k = 1:N % per each robot
        r = Robots(idx_use(k)); % select k-th robot
        working_state = lower(r.working_state);
        switch working_state
            case 'p' % point (adaptive sigma)
                sigma0 = sigma_point;
                sigmaMax = 5*sigma0;
                mu = pick_target(r, robot_pos_est(k,:)); % mean value of Phi must be the robot's target --> lineup point
                p  = robot_pos_est(k,:);  % current robot position
                d  = norm(p - mu);  % distance to the target
                sigma_k = ((d + R)/R) * sigma0;  % adaptive sigma = sigma0*(1 + d/R)
                if sigma_k > sigmaMax
                   sigma_k = sigmaMax;                   % saturation
                end
                Phi(:,:,k) = gauss2d(X, Y, mu, sigma_k^2, sigma_k^2); % build Phi (bivariate Gaussian)

            case 'r' % ring (adaptive sigma)
                assert(R>0,'R must be > 0 for adaptive sigma.');
                sigma0   = sigma_ring;    % sigma near ring
                sigmaMax = 3*sigma0;      % limit on sigma
                mu_pkg = pick_target(r, robot_pos_est(k,:));
                p  = robot_pos_est(k,:); 
                d  = abs(norm(p - mu_pkg) - R); % distance from ring
                sigma_k = ((d + R)/R) * sigma0; % adaptive sigma
                if sigma_k > sigmaMax
                    sigma_k = sigmaMax;  % saturation
                end
                Phi(:,:,k) = ring2d(X, Y, mu_pkg, R, sigma_k, 1e-15);
                
            case 't' % transportation TO BE TESTED
                    Phi(:,:,k) = 0; % no centroid
            case 'f' % free
                    Phi(:,:,k) = ones(m,n);  % uniform pdf --> optimal coverage (classic Voronoi)
            case 'i' % inbound attraction
                    mu = [max(X(1,:))/2 5];
                    sigma0 = sigma_point*4;
                    sigmaMax = 4*sigma0;
                    p  = robot_pos_est(k,:);  % current robot position
                    d  = norm(p - mu);  % distance to the target
                    sigma_k = ((d + R)/R) * sigma0;  % adaptive sigma = sigma0*(1 + d/R)
                    if sigma_k > sigmaMax
                        sigma_k = sigmaMax;                   % saturation
                    end
                    Phi(:,:,k) = gauss2d(X, Y, mu, sigma_k^2, sigma_k^2);
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
        % Reactive control --> use Wrs instead of whole voronoi cell
        if use_reactive
            mask_k = Wrs_set{k};
        end

        w = Phi(:,:,k); % k-th robot's pdf
        w(~mask_k) = 0; % must be zero outside mask --> only integrate it on the mask : this is coherent with theory (only integrate on the voronoi cell's domain)

        M = sum(w(:)) * cell_area; % integral of Phi(x,y)*dxdy  =  sum(Phi)*dA
        masses(k) = M;

            if M <= 1e-50 || lower(Robots(idx_use(k)).working_state)=='t'% !!! this has to be fine tuned !!! if mass is null --> no centroid
                centroids(k,:) = [NaN NaN];
            else % define centroid
                Cx = sum(sum(X .* w)) * cell_area / M; % look at theory
                Cy = sum(sum(Y .* w)) * cell_area / M;
                centroids(k,:) = [Cx, Cy];
            end
            % projection of centroid on Wrs
            if use_reactive && all(isfinite(centroids(k,:)))
                centroids(k,:) = project_centroid(centroids(k,:), Wrs_set{k}, X, Y);
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
