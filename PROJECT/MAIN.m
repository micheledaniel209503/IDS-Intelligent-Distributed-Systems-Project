%% IDS PROJECT - WAREHOUSE
% GIANLUCA LONA / MICHELE DANIEL

% NOTE: Running the whole file will close the figures related to the single
% steps and will only display the simulations. If one wants to analyze the
% intermediate steps, execute the subsections one at a time.


clc
close all
clear all

% add path to the class, functions folders
addpath(genpath(pwd));

toRad = pi / 180; % Conversion factor from degrees to radians
toDeg = 180 / pi; % Conversion factor from radians to degrees

%% Parameters
tolerance_on_lineup_form_radius = 0.3; % radial tolerance for the initial position of robots in lineup phase

%% Building the environment
% --------DEFINITIONS---------
width = 100;
height = 80;
% ---------MAP-----------
map.W = width;
map.H = height;
map.polygon = [0 0; map.W 0; map.W map.H; 0 map.H];
% inbound zone
zones.inbound.W = 30;
zones.inbound.H = 10;
zones.inbound.polygon = [
    map.W/2-zones.inbound.W/2 0,
    map.W/2+zones.inbound.W/2 0, 
    map.W/2+zones.inbound.W/2 zones.inbound.H, 
    map.W/2-zones.inbound.W/2 zones.inbound.H
    ];
%% package spawn
Npack = 1; % number of packages
Packages(Npack,1) = pkg();
for k = 1:Npack
    Packages(k).id = k;
    [Packages(k).state(1), Packages(k).state(2)] = spawn_pkg(zones.inbound);
end
%% Robots spawn
% create the robots, placing them in random positions, avoiding the inbound
% zone
Nrobots = 15;
Robots(Nrobots,1) = rob(); % list of rob objects

% dimension factor to set the width of the robots position spawn in x
dimensionFactor = 0.8; % Example factor for robot width adjustment
init_pos = zeros(Nrobots, 2);

for k = 1:Nrobots
    Robots(k).id = k;
    valid_pos = false;
    while ~valid_pos
        % random spawn in the map (close to the inbound zone)
        x = map.W/2 * (1-dimensionFactor) + dimensionFactor * map.W * rand();
        y = zones.inbound.H + randn * 0.03 * map.H;
        init_pos(k, :) = [x y];
        % check if inside inbound zone
        in_x = (x >= 0);
        in_y = (y >= zones.inbound.H);

        if  (in_x && in_y)
            Robots(k).state(1) = x;
            Robots(k).state(2) = y;
            valid_pos = true;
        end
    end
end

% Robot sensing radius (for the consensus)
for i=1:Nrobots
    Robots(i).sr = 15;
end
% Number of messages exchanged
n_msg = 12 ;

% Measurement noise
sigma = 0.5; % std of distance noise

%% Initial measure of package position (robots only if they see it)
x0 = nan(Nrobots,2); % [x_est, y_est]
distances = zeros(Nrobots, 1);

for i = 1:Nrobots
    dx = Packages(1).state(1) - Robots(i).state(1);
    dy = Packages(1).state(2) - Robots(i).state(2);
    dist = sqrt(dx^2 + dy^2);
    distances(i) = dist;

    if dist <= Robots(i).sr
        % Noisy distance measurement -> backproject as rough estimate
        noisy_dist = dist + sigma*randn();
        % Assume robot points directly toward the package (simplification)
        est_x = Robots(i).state(1) + noisy_dist * dx/dist;
        est_y = Robots(i).state(2) + noisy_dist * dy/dist;
        x0(i,:) = [est_x, est_y];
    end
end

% Find robots that can see the package
idxNonNaN = find(~isnan(x0(:,1)));
nNonNaN = numel(idxNonNaN);
disp('Robots that see the package (IDs): '); disp(idxNonNaN);
disp('Number of robots that see the package: '); disp(nNonNaN);

%% --------------------------------
close all; % close previous figures





%  --------------------------------

%% CONSENSUS ALGORITHM (broadcast one at a time) - only among robots that see the package
% Preallocate storage for all iterations (for all robots; NaN for non-seers)
xStore = nan(Nrobots,2,n_msg+1);
xStore(:,:,1) = x0;

% If no robot sees the package or only one sees it, no active consensus
if nNonNaN == 0
    warning('No robot sees the package.');
    return;
end

% Extract the initial subset (only those that see)
% xStore_sub holds the estimates of the participating robots only
xStore_sub = nan(nNonNaN,2,n_msg+1);
x0_sub = x0(idxNonNaN,:);            % initial estimates [nNonNaN x 2]
xStore_sub(:,:,1) = x0_sub;

% alpha relative to the participating subset (ensures double-stochastic matrices)
% alpha_sub = 1/(nNonNaN-1) is the standard choice to guarantee nonnegative diagonal
alpha_sub = 1/(nNonNaN-1);
Q_tot = eye(nNonNaN);

for t = 1:n_msg
    
    % choose the transmitter in round-robin among the robots that see the package
    % i_local is a LOCAL index (1..nNonNaN), not the global robot id
    i_local = mod(t-1, nNonNaN) + 1;

    % Build Qi_sub (size nNonNaN x nNonNaN)
    % Qi_sub implements the cyclic broadcast gossip: the selected node broadcasts,
    % others update as x_j <- (1-alpha)*x_j + alpha*x_i, and the broadcaster's
    % own row is set so that Qi_sub is doubly stochastic (so average is preserved).
    Qi_sub = zeros(nNonNaN);
   
    % Distinguish the case when we have only 2 values, if we keep the
    % standard algorithm, the solution oscillates
    if nNonNaN == 2
        Qi_sub = 1/2 * ones(2)
    else
        for j = 1:nNonNaN
            if j ~= i_local
                Qi_sub(j,j)      = 1 - alpha_sub;   % keeps part of own value
                Qi_sub(j,i_local)= alpha_sub;       % receives fraction from i_local
                Qi_sub(i_local,j)= alpha_sub;       % symmetric term to maintain double-stochasticity
            end
        end
    % close the balance on the broadcaster diagonal entry
    Qi_sub(i_local,i_local) = 1 - alpha_sub*(nNonNaN-1); % ensures rows sum to 1
    end

    % Debug: show which local index transmits and the Qi_sub matrix
    %disp(['Transmitter = ', num2str(i_local)]);
    %disp(Qi_sub);
    Q_tot = Qi_sub*Q_tot;

    x_next_sub = Qi_sub * xStore_sub(:,:,t);   % [nNonNaN x 2]
    xStore_sub(:,:,t+1) = x_next_sub;

    % Also update the global storage so we can plot global positions later
    % Non-seers remain unchanged (NaN), participating robots are overwritten
    xStore(:,:,t+1) = xStore(:,:,t);
    xStore(idxNonNaN,:,t+1) = x_next_sub;
end
%disp(Q_tot);

% assign to the package its own state estimate
Packages(1).state_est = mean(xStore_sub(:,:,end),1);

disp(['Real package position (x, y): ', num2str(Packages(1).state(1)), '  ', num2str(Packages(1).state(2))]);
disp(['Est. package position (x, y): ', num2str(xStore_sub(1,1,end)), '  ',num2str(xStore_sub(1,2,end))]);
disp(['error (distance) [cm]: ', num2str(100*sqrt((Packages(1).state(1)-xStore_sub(1,1,end))^2+(Packages(1).state(2)-xStore_sub(1,2,end))^2))])



%% Plot results on the map 
figure(1); clf; hold on; axis equal;
% Map boundary
plot([0 map.W map.W 0 0],[0 0 map.H map.H 0],'k-'); % outline

% Inbound zone
fill(zones.inbound.polygon(:,1),zones.inbound.polygon(:,2), ...
     [0.9 0.9 1],'EdgeColor','b');

% Robots
S = vertcat(Robots.state);   % N×2: each row = [x y] state of robot
plot(S(:,1), S(:,2), 'ko', 'MarkerFaceColor','w');
%text([Robots.state(1)]+1,[Robots.state(2)]+1,string(1:Nrobots))

% Package true position
hPkg = plot(Packages(1).state(1),Packages(1).state(2),'rs','MarkerFaceColor','r','MarkerSize',10);

% Final consensus estimates (only those who participated)
hEst = plot(xStore(:,1,end),xStore(:,2,end),'bo','MarkerFaceColor','b');

% --- Highlight robots that sense the package ---
theta = linspace(0,2*pi,100);
for i = idxNonNaN(:,1)'
        cx = Robots(i).state(1) + Robots(i).sr*cos(theta);
        cy = Robots(i).state(2) + Robots(i).sr*sin(theta);
        plot(cx,cy,'k:','LineWidth',1)  % dotted circle showing sensing radius
end

% --- Legend (only package and estimates) ---
legend([hPkg, hEst], {'Package (true)','Consensus estimates'}, 'Location','best')
title('Map: robot positions, package and final estimates')

%% Plot: evolution of x and y estimates for each iteration
if nNonNaN >= 1
    iterations = 0:n_msg; % t=0 (initial) ... n_msg
    figure(2); clf;

    % Subplot X
    subplot(2,1,1); hold on; grid on;
    for k = 1:nNonNaN
        plot(iterations, squeeze(xStore_sub(k,1,:)), 'LineWidth', 1.5);
    end
    % Add true X value as reference line
    yline(Packages(1).state(1), 'k--', 'LineWidth', 2, 'DisplayName','True X');

    xlabel('Iteration (message step)');
    ylabel('Estimate x (coordinate)');
    title('Evolution of X estimates for robots that see the package');
    legendStrings = arrayfun(@(r) sprintf('Robot %d', r), idxNonNaN, 'UniformOutput', false);
    legend([legendStrings], 'Location','best');

    % Subplot Y
    subplot(2,1,2); hold on; grid on;
    for k = 1:nNonNaN
        plot(iterations, squeeze(xStore_sub(k,2,:)), 'LineWidth', 1.5);
    end
    % Add true Y value as reference line
    yline(Packages(1).state(2), 'k--', 'LineWidth', 2, 'DisplayName','True Y');

    xlabel('Iteration (message step)');
    ylabel('Estimate y (coordinate)');
    title('Evolution of Y estimates for robots that see the package');
    legend([legendStrings], 'Location','best');
end

%% --------------------------------
close all; % close previous figures





%  --------------------------------
%% ORGANIZE TRANSPORTATION
% Based on surface area of the package --> n° robots needed
Packages(1).s = 20; % [m^2] surface area of the package
robot_sc = 4; % [m^2] surface capacity of a robot

disp('--- Organizing package transportation ---');
RN = ceil(Packages(1).s / robot_sc);
disp(['Robots needed to carry package 1: ', num2str(RN)])

% Which robots will go? The ones that are closer
% Considering that, after consensous on package position, all robots can share
% the information regardless of which ones see the package
% Find the indices of the first RN smallest distances
% this only works because robots id = robots indices
[~, sortedIndices] = sort(distances);
Robots_selected_id = sortedIndices(1:RN);

% Associate the package to the robots, using the field item_id
for i=1:length(Robots_selected_id)
    Robots(Robots_selected_id(i)).item_id = 1; % change with more than one package
end

disp('Indices of robots selected to carry the package: '); 
disp(Robots_selected_id);

for j = 1:length(Robots)
    Robots(j).state(3) = 0; % [rad] set robot's orientation to 0
end

radius = ceil(sqrt(Packages(1).s / pi)*1.5); % compute radius from Packages(1).s circular surface

disp('Radius of the package ');
disp([radius]);
Packages(1).r = radius;

%% LLOYD-VORONOI simulation setup parameters
% Simulation parameters
dt = 0.1; % [s] time step
T_sim = 100; % [s] sim time
iter_sim = T_sim / dt; % Calculate the number of iterations
disp(['Simulation Setup...']);

% map grid
dx = 0.5;
dy = 0.5;
xv = 0:dx:100;   % X axis
yv = 0:dy:80;    % Y axis
[X,Y] = meshgrid(xv, yv); % grid
free_mask = true(size(X)); % free grid (no obstacles) --> matrix of ones


Robots_voronoi = [Robots(:).id]'; % all robots in the tasselation
sigma_ring = 1.25;
sigma_lineup = 1.5;
sigma_transport = 1.0;
u_sat = 2.5/3.6; % [m/s]
omega_sat = 30*toRad; % [rad/s]
Kp_u     = 0.8;          % gain
Kp_theta = 3.0;          % gain
r_goal   = 0.0;            % [m] distance from target at which to stop
theta_goal = deg2rad(0); % [rad] final robot orientation

% optional keep theta bounded
wrap = @(a) atan2(sin(a),cos(a));
traj = zeros(iter_sim,3,numel(Robots_selected_id));  % log

%% Debugging the pdfs
ring_pdf = ring2d(X, Y, Packages(1).state, Packages(1).r, 2, 1e-3);
figure
surf(xv,yv,ring_pdf); shading interp; colorbar;

% opts = struct('logscale', true, 'mode','ring', 'mu_pkg', mu_pkg, 'R', Packages(1).r);
% figure; plot_phi_debug(Phi, X, Y, L, Robots, idx_use, centroids, opts);
%% INTRODUCING UNCERTAINTY
% let's introduce uncertainty and use the EKF

% Define the system matrix
% Unicycle dynamics
fun2 = @(x, y, theta, vel, omega, dT) [x + vel * cos(theta) * dT; y + vel * sin(theta) * dT; theta + omega * dT]; 
A = @(x, y, theta, vel, omega, dT) [1, 0, -vel * sin(theta) * dT; 0, 1, vel * cos(theta) * dT; 0, 0, 1];
G = @(x, y, theta, vel, omega, dT)[ dT*cos(theta), 0; dT*sin(theta), 0; 0,  dT ];
% Uncertainties
noise_std = 0.2; % [m] std on the distance measures with the anchors
sigma_theta = 0.02; % [RAD] std on the orientation measure 
Q = 0.1*eye(2); % process noise covariance (set to zero if no uncertainty on the model)

% Choose N of anchors
Nanchors = 5;
anchors = zeros(Nanchors,2);
anchors(:,1) = map.W * rand(Nanchors,1);
anchors(:,2) = map.H * rand(Nanchors,1);
for a = 1:Nanchors
    ax = anchors(a,1); 
    ay = anchors(a,2);
end

%% Generate random obstacles in the grid
Nobs = 3;
obstacle_dim = 13.0; % [m] characteristic dimension of the obstacles (global)
min_dist = 12; % [m] distance between centers of the obstacles

Obstacles = spawn_obstacles(Nobs, map, obstacle_dim, min_dist);

% Adapt the grid
free_mask = true(size(X)); % initialization : free grid (no obstacles) --> matrix of ones
[Xg, Yg] = deal(X, Y); % alias
for k = 1:numel(Obstacles)
     o = Obstacles(k);
     ox = o.state(1); oy = o.state(2);
     switch o.type
          case 'c'
                r = (o.l)/2;
                mask = (Xg - ox).^2 + (Yg - oy).^2 <= r^2;
          case 's'
                half = o.l/2;
                mask = (abs(Xg - ox) <= half) & (abs(Yg - oy) <= half);
     end
     free_mask(mask) = false;
end

%% initial position of robots
for j = 1:length(Robots_voronoi)
    idx = Robots_voronoi(j); % select robot id
    Robots(idx).state(1:2) = init_pos(j, :); % initialize robot positions true
    Robots(idx).state_est(1:3) = [NaN NaN Robots(idx).state(3)];  % initialize robot positions estimate
end

for i = 1:length(Robots_selected_id)
    idx = Robots_selected_id(i);
    Robots(idx).working_state = 'r'; % ring working state test
    Robots(idx).target = Packages(1).state(1:2);
end

% Select two robots that are NOT in Robots_selected_id and send them on two
% points
Robots_unselected_id = setdiff(1:Nrobots, Robots_selected_id);
idx1 = Robots_unselected_id(1);
idx2 = Robots_unselected_id(end);
Robots(idx1).working_state = 'p';
Robots(idx2).working_state = 'p';
Robots(idx1).target = [70, 70];
Robots(idx2).target = [80, 10];


%% Initialisation for the estimation algorithm (recursive least square)
for i = 1:length(Robots_voronoi)
    distances = sqrt(sum((anchors - [Robots(i).state(1), Robots(i).state(2)]).^2, 2));
    distances_noisy = distances + noise_std * randn(Nanchors, 1);

    [H,z,C] = trilateration(anchors, distances_noisy, 0.1);
    P = (H'*C^-1*H)^-1;
    state_ls = P*H'*C^-1*z;
    Robots(i).state_est(1) = state_ls(1);
    Robots(i).state_est(2) = state_ls(2);
    Robots(i).state_est(3) = Robots(i).state(3); % initial estimate on initial theta = exact initial theta
  
    P_values = [P, zeros(2, 1); zeros(1, 3)];
    P_values(3,3) = 0.1;
    Robots(i).P = P_values;
    
end

trace_P = zeros(iter_sim); % initialization;

%% LINEUP SIMULATION with SENSING RADIUS

disp('--- STARTING LINEUP...');

USEREACTIVE = true;   % flag: if false --> use standard voronoi (no cutting for the cells, no obstacle avoidance)
DRAWCONTOUR = true;   % flag: draw contours of the voronoi cells
rsense      = 10.0;  % [m] sensing radius for the robots (can be reduced in order to manage the collision avoidance in a different way)

for i=1:numel(Robots)
    Robots(i).sr = rsense;
end

[L, areas, masses, centroids, Phi, Wrs_set] = voronoi_lloyd_ring_dyna_unc_reactive(Robots, Robots_voronoi, X, Y, free_mask, sigma_lineup, sigma_ring, Packages(1).r, USEREACTIVE);

% plot
RN   = numel(Robots_voronoi);

figure(101); clf
imagesc(xv, yv, free_mask); set(gca,'YDir','normal'); colormap(gray); axis equal tight; hold on
hTitle = title(sprintf('RING LINEUP (WITH UNCERTAINTIES): 0.0%%'));
cmap = lines(numel(Wrs_set));
if DRAWCONTOUR == true
    for k = 1:numel(Wrs_set)
        if any(Wrs_set{k}(:))
            contour(xv, yv, double(Wrs_set{k}), [0.5 0.5], ...
                    'Color', cmap(k,:), 'LineWidth', 1.4);
        end
    end
end
hold on

    % plot stuff
    hRob  = gobjects(RN,1);
    hCent = gobjects(RN,1);
    center = Packages(1).state(1:2); % plot ring
    tt = linspace(0, 2*pi, 361);
    xr = center(1) + Packages(1).r*cos(tt);
    yr = center(2) + Packages(1).r*sin(tt);
    for j = 1:RN % initial step
        id = Robots_voronoi(j);
        p  = Robots(id).state(1:2);
        hRob(j)  = plot(p(1), p(2), 'o', 'MarkerSize',6, ...
                        'MarkerFaceColor',cmap(j,:), 'Color','k');
        hCent(j) = plot(NaN, NaN, 'r+', 'MarkerSize',10, 'LineWidth',1.2);
        hRing = plot(xr, yr, 'k--', 'LineWidth', 1);
        hPkg = plot(center(1), center(2), 'bs', 'MarkerSize', 10);
    end

    % refresh plot every iteration
    refresh_every = 1;          % robot/centrois
    refresh_contours_every = 20; % contours of Wrs


for k = 1:iter_sim

    % state estimation
    for j = 1:numel(Robots_voronoi)
        % Measure of distances wrt anchors and orientation
        distances = sqrt(sum((anchors - [Robots(j).state(1), Robots(j).state(2)]).^2, 2));
        distances_noisy = distances + noise_std * randn(Nanchors, 1);
        % trilateration
        [H_tril,z_tril,R_tril] = trilateration(anchors, distances_noisy, noise_std);
        z_theta = Robots(j).state(3) + sigma_theta*randn();
        % Extend H to include theta state and add orientation measurement row
        m = size(H_tril,1);
        H_extended = [H_tril, zeros(m,1)];  
        H_theta = [0, 0, 1];                   
    
        % Stack H and z
        H = [H_extended; H_theta];         
        z = [z_tril; z_theta];

        R = blkdiag(R_tril, sigma_theta^2);

        % Estimation of position and orientation using EKF (based on the last control input)
        [Robots(j).state_est, Robots(j).P] = EKF_function(Robots(j).state_est, Robots(j).u, Robots(j).omega, fun2, A, G, z, H, Robots(j).P, Q, R, dt);
        % Update the trace of the covariance matrix
        trace_P(k) = trace_P(k) + trace(Robots(j).P);
    end

    % compute Voronoi partitioning + centroids for each robot
    [L, areas, masses, centroids, Phi, Wrs_set] = voronoi_lloyd_ring_dyna_unc_reactive( ...
        Robots, Robots_voronoi, X, Y, free_mask, ...
        sigma_lineup, sigma_ring, Packages(1).r, USEREACTIVE);

    % check the robots lineup for every package
    for i = 1:Npack
        % Get all robots carrying package i
        robots_id_itemId_i = [];  % initialize empty list

        for r = 1:length(Robots)
            % Skip robots with empty or NaN item_id
            if isempty(Robots(r).item_id) || isnan(Robots(r).item_id)
                continue;
            end
        
            % Check if this robot carries the i-th package
            if Robots(r).item_id == i
                robots_id_itemId_i = [robots_id_itemId_i, Robots(r).id];
            end
        end

        % Skip if no robots are assigned to this package
        if isempty(robots_id_itemId_i)
            continue;
        end
        
        % Check alignment for each robot
        allAligned = true; 
        
        for id = robots_id_itemId_i
            % Skip if robot is already in transport phase
            if Robots(id).working_state == 't'
                continue;
            end
            
            % Check alignment for this robot
            isAligned = isRobotAligned(Robots(id), Packages(i).state_est, Packages(i).r, tolerance_on_lineup_form_radius);
            % If even one robot is not aligned, stop checking further
            if ~isAligned
                allAligned = false;
                break;
            end
        end
        
        % If all robots are aligned, switch their state to t
        if allAligned
            for id = robots_id_itemId_i
                Robots(id).working_state = 't';
            end
        end
    end

    % robot control + dynamics
    for j = 1:numel(Robots_voronoi)
        id = Robots_voronoi(j);
        ci = centroids(j,:);

        % control
        if Robots(j).working_state ~= 't'
            [u, omega] = ROB_control(Robots(id).state_est, ci, ...
                                     u_sat, omega_sat, Kp_u, Kp_theta, ...
                                     r_goal, theta_goal);
        end
        if Robots(j).working_state == 't'
            u=0;
            omega=0;
        end
        Robots(id).u = u;
        Robots(id).omega = omega;
        
        % dyna
        %Robots(id).state = fun(Robots(id).state, u, omega, dt);
        Robots(id).state = fun2(Robots(id).state(1), Robots(id).state(2), Robots(id).state(3), u, omega, dt);

        % wrap theta
        Robots(id).state(3) = wrap(Robots(id).state(3));
    end

        % plot
        if mod(k, refresh_every) == 0
            perc = 100*k/iter_sim;
            set(hTitle, 'String', sprintf('SIMULATION (WITH UNCERTAINTIES): %.1f%%', perc));
            % robots + centroids
            for j = 1:RN
                id = Robots_voronoi(j);
                px = Robots(id).state(1); py = Robots(id).state(2);
                set(hRob(j), 'XData', px, 'YData', py); % robot marker
                cx = centroids(j,1); cy = centroids(j,2);
                set(hCent(j), 'XData', cx, 'YData', cy);
            end
    
            if (mod(k, refresh_contours_every) == 0 && DRAWCONTOUR == true)
                ax = gca;
                hOld = findall(ax,'Type','Contour');
                delete(hOld);

                for kk = 1:numel(Wrs_set)
                    if any(Wrs_set{kk}(:))
                        contour(xv, yv, double(Wrs_set{kk}), [0.5 0.5], ...
                        'Color', cmap(kk,:), 'LineWidth', 1.4);
                    end
                end
            end
            drawnow limitrate nocallbacks
        end
end

disp('--- LINEUP ENDED !');