%% IDS PROJECT - WAREHOUSE
% GIANLUCA LONA / MICHELE DANIEL

% NOTE: Running the whole file will close the figures related to the single
% steps and will only display the simulations. If one wants to analyze the
% intermediate steps, execute the subsections one at a time.


clc
close all
clear all

%% Auxiliary plot flags
PLOT_CONSENSUS = 1;   % plot consensus estimate
PLOT_STATEERROR = 1;  % plot state error in time
PLOT_TRACEP = 1;      % plot state covariance matrix
PLOT_FORM_ERR = 1;    % flag: plot formation error metrics
USEREACTIVE = true;   % flag: if false --> use standard voronoi (no cutting for the cells, no obstacle avoidance)
DRAWCONTOUR = true;   % flag: draw contours of the voronoi cells
%% Plot variables
N_iter = 1000;
estState = zeros(3,N_iter);
realState = zeros(3,N_iter);
trP = zeros(1, N_iter);
trP_delivered = [];
form_err = zeros(1, N_iter);
form_err_idx = 1;
% add path to the class, functions folders
addpath(genpath(pwd));

toRad = pi / 180; % Conversion factor from degrees to radians
toDeg = 180 / pi; % Conversion factor from radians to degrees

%% Parameters
%% Number of messages exchanged in the consensus algorithm
n_msg = 12 ;
%% Measure Uncertainties
sigma_UWB = 0.2; % [m] std on the distance measures with the anchors
sigma_theta = 3*toRad; % [RAD] std on the orientation measure 
sigma_LiDAR = 0.1; % [m] std of distance noise
%% Model Uncertainties
% Define the system matrix
% Unicycle dynamics
fun2 = @(x, y, theta, vel, omega, dT) [x + vel * cos(theta) * dT; y + vel * sin(theta) * dT; theta + omega * dT]; 
A = @(x, y, theta, vel, omega, dT) [1, 0, -vel * sin(theta) * dT; 0, 1, vel * cos(theta) * dT; 0, 0, 1];
G = @(x, y, theta, vel, omega, dT)[ dT*cos(theta), 0; dT*sin(theta), 0; 0,  dT ];

% set process noise variances for v and omega
sigma_v = 0.2;        % [m/s] 
sigma_omega = 0.1;    % [rad/s] 
Q = diag([sigma_v^2, sigma_omega^2]);

%% Formation Control Parameters
u_max = 1;         % [m/s]
omega_max = 0.2;     % [rad/s]

w_form  = 0.2; % FORMATION TASK
w_att = 5.0; % TARGET TASK
w_obs = 10.0;   % OBSTACLE AVOIDANCE

%% Simulation parameters
dt = 0.1; % [s] time step
T_sim = 1000; % [s] sim time
iter_sim = T_sim / dt; % Calculate the number of iterations

%% Building the environment
% --------DEFINITIONS---------
width = 100; % [m]
height = 80; % [m]
% ---------MAP-----------
map.W = width;
map.H = height;
map.polygon = [0 0; map.W 0; map.W map.H; 0 map.H];
% inbound zone
zones.inbound.W = 30; % [m]
zones.inbound.H = 10; % [m]
zones.inbound.polygon = [
    map.W/2-zones.inbound.W/2 0,
    map.W/2+zones.inbound.W/2 0, 
    map.W/2+zones.inbound.W/2 zones.inbound.H, 
    map.W/2-zones.inbound.W/2 zones.inbound.H
    ];
% outbound zone
zones.outbound.W = 20; % [m]
zones.outbound.H = 10; % [m]
zones.outbound1.polygon = [
    0+5, map.H - zones.outbound.H;
    zones.outbound.W + 5, map.H - zones.outbound.H;
    zones.outbound.W + 5, map.H;
    0+5, map.H
];

zones.outbound2.polygon = [
    map.W/2 - zones.outbound.W/2, map.H - zones.outbound.H;
    map.W/2 + zones.outbound.W/2, map.H - zones.outbound.H;
    map.W/2 + zones.outbound.W/2, map.H;
    map.W/2 - zones.outbound.W/2, map.H
];

zones.outbound3.polygon = [
    map.W - zones.outbound.W - 5, map.H - zones.outbound.H;
    map.W - 5, map.H - zones.outbound.H;
    map.W - 5, map.H;
    map.W - zones.outbound.W - 5, map.H
];
%% Anchors
% Choose N of anchors
Nanchors = 4;
anchors = zeros(Nanchors,2);
anchors(:,1) = map.W * rand(Nanchors,1);
anchors(:,2) = map.H * rand(Nanchors,1);
for a = 1:Nanchors
    ax = anchors(a,1); 
    ay = anchors(a,2);
end
%% Generate random obstacles in the grid
Nobs = 4;
obstacle_dim = 8.0; % [m] characteristic dimension of the obstacles
min_dist = 10 + obstacle_dim; % [m] distance between centers of the obstacles
% spawn obstacles
Obstacles = spawn_obstacles(Nobs, map, obstacle_dim, min_dist);
%% package spawn
Npack = 3; % number of packages
current_pkg = 1;
Packages(Npack,1) = pkg();
for k = 1:Npack
    Packages(k).id = k;
    [Packages(k).state(1), Packages(k).state(2)] = spawn_pkg(zones.inbound);
end

% error on localization of each pkg (euclidean norm between eestimate and
% real)
pkg_loc_error = NaN(Npack,1); % [m]
%% Robots spawn
% create the robots, placing them in random positions, avoiding the inbound
% zone
Nrobots = 9; % number of robots in the map
Robots(Nrobots,1) = rob(); % list of rob objects

%% Initialialisation of robots position in the map
% dimension factor to set the width of the robots position spawn in x
dimensionFactor = 0.8;
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
            Robots(k).state(3) = -pi + 2 * pi * rand();
            valid_pos = true;
        end
    end
end

% ROBOT TRAJECTORIES
traj_x = nan(Nrobots, iter_sim+1);
traj_y = nan(Nrobots, iter_sim+1);


% log working_state (i, f, l, r, t, etc.)
state_log = repmat('u', Nrobots, iter_sim+1); % 'u' = unknown

% initial state
for j = 1:Nrobots
    traj_x(j,1)    = Robots(j).state(1);
    traj_y(j,1)    = Robots(j).state(2);
    if isempty(Robots(j).working_state)
        state_log(j,1) = 'u';
    else
        state_log(j,1) = Robots(j).working_state;
    end
end

%% WAREHOUSE TOPOLOGY PLOT
Packages(1).s = 20;
figure('Color','w'); hold on; axis equal; box on;
xlabel('x [m]'); ylabel('y [m]');
title('Warehouse topology');
% background
hMap = fill(map.polygon(:,1), map.polygon(:,2), [0.95 0.95 0.95], ...
'EdgeColor','k', 'LineWidth',1.5);
% inbound zone
hIn = fill(zones.inbound.polygon(:,1), zones.inbound.polygon(:,2), ...
'b', 'FaceAlpha',0.15, 'EdgeColor','b', 'LineWidth',1.5);
% outbound zones
hOut1 = fill(zones.outbound1.polygon(:,1), zones.outbound1.polygon(:,2), ...
'r', 'FaceAlpha',0.15, 'EdgeColor','r', 'LineWidth',1.5);
hOut2 = fill(zones.outbound2.polygon(:,1), zones.outbound2.polygon(:,2), ...
'r', 'FaceAlpha',0.15, 'EdgeColor','r', 'LineWidth',1.5);
hOut3 = fill(zones.outbound3.polygon(:,1), zones.outbound3.polygon(:,2), ...
'r', 'FaceAlpha',0.15, 'EdgeColor','r', 'LineWidth',1.5);
% obstacles
hObs_proto = [];
for k = 1:numel(Obstacles)
o  = Obstacles(k);
ox = o.state(1);
oy = o.state(2);
switch o.type
    case 'c' % circle
r = o.l/2;
th = linspace(0, 2*pi, 60);
x_o = ox + r*cos(th);
y_o = oy + r*sin(th);
hTmp = fill(x_o, y_o, 'k', 'EdgeColor','k');
case 's' % square
half = o.l/2;
x_o = [ox-half, ox+half, ox+half, ox-half];
y_o = [oy-half, oy-half, oy+half, oy+half];
hTmp = fill(x_o, y_o, 'k', 'EdgeColor','k');
end
if isempty(hObs_proto)
hObs_proto = hTmp;
end
end
% robots
Nrobots = numel(Robots);
rob_pos = zeros(Nrobots,2);
for i = 1:Nrobots
rob_pos(i,:) = Robots(i).state(1:2);
end
hRob = plot(rob_pos(:,1), rob_pos(:,2), 'o', ...
'MarkerSize',8, 'MarkerFaceColor',[0 0.6 1], ...
'MarkerEdgeColor','k');
% package
pkg_pos  = Packages(1).state_est(1:2);
pkg_side = sqrt(Packages(1).s);
half_s   = pkg_side/2;
x_pkg = [pkg_pos(1)-half_s, pkg_pos(1)+half_s, ...
pkg_pos(1)+half_s, pkg_pos(1)-half_s];
y_pkg = [pkg_pos(2)-half_s, pkg_pos(2)-half_s, ...
pkg_pos(2)+half_s, pkg_pos(2)+half_s];
hPkg = fill(x_pkg, y_pkg, [0.6 0.6 0.6], ...
'EdgeColor','k', 'LineWidth',1.5);
% anchors
hAnc = plot(anchors(:,1), anchors(:,2), '^', ...
    'MarkerSize',8, 'MarkerFaceColor','y', ...
    'MarkerEdgeColor','k', 'LineWidth',1.2);

xlim([0 map.W]); ylim([0 map.H]);
% legend
legend([hIn, hOut1, hObs_proto, hRob, hPkg, hAnc], ...
    {'Inbound zone', 'Outbound zones', 'Obstacles', 'Agents', 'Package', 'Anchors'}, ...
    'Location','bestoutside');


%% VORONOI control parameters
Robots_voronoi = [Robots(:).id]'; % all robots in the tasselation
sigma_ring = 1.25;
sigma_point = 1.5;
sigma_transport = 1.0;
u_sat = 3/3.6; % [m/s]
omega_sat = 30*toRad; % [rad/s]
Kp_u     = 0.8;          % gain
Kp_theta = 3.0;          % gain
r_goal   = 0.0;            % [m] distance from target at which to stop
theta_goal = deg2rad(0); % [rad] final robot orientation

% optional keep theta bounded
wrap = @(a) atan2(sin(a),cos(a));
%% Initialisation of state estimation for robots
for j = 1:length(Robots_voronoi)
    idx = Robots_voronoi(j);

    %% Initialisation for the estimation algorithm (weighted least square)
    % true distances to anchors
    distances = sqrt(sum((anchors - [Robots(idx).state(1), Robots(idx).state(2)]).^2, 2));
    % noisy distance measurements
    distances_noisy = distances + sigma_UWB * randn(Nanchors, 1);

    % trilateration
    [H_xy, z_xy, C_xy] = trilateration(anchors, distances_noisy, sigma_UWB);

    % simulated noisy measurement
    theta_meas = Robots(idx).state(3) + sigma_theta * randn(); 
    theta_meas = atan2(sin(theta_meas), cos(theta_meas));

    % build augmented H, z, C for [x; y; theta]
    H_theta = [0 0 1];           
    z_theta = theta_meas;       
    C_theta = sigma_theta^2; 

    H_temp = [H_xy, zeros(size(H_xy,1),1)]; 

    % combined matrices
    H = [H_temp; H_theta];      
    z = [z_xy; z_theta];         
    C = blkdiag(C_xy, C_theta); 

    % compute least-squares estimate and covariance for [x;y;theta]
    P = (H'*C^-1*H)^-1;
    state_ls = P*H'*C^-1*z; 

    % assign estimates
    Robots(idx).state_est(1) = state_ls(1);
    Robots(idx).state_est(2) = state_ls(2);
    Robots(idx).state_est(3) = atan2(sin(state_ls(3)), cos(state_ls(3)));

    Robots(idx).P = P;

end

%% LLOYD-VORONOI simulation setup parameters
disp(['Simulation Setup...']);

% map grid
dx = 0.5;
dy = 0.5;
xv = 0:dx:100;   % X axis
yv = 0:dy:80;    % Y axis
[X,Y] = meshgrid(xv, yv); % grid
free_mask = true(size(X)); % free grid (no obstacles) --> matrix of ones


%% PLACE OBSTACLES in the MAP
% Adapt the grid
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

% initialise simulation plot
figure(); clf


% SIMULATING Npack DELIVERIES
for pkg_i = 1:Npack

    for j = 1:numel(Robots)
        Robots(j).sr = 20;
    end

    current_pkg = pkg_i;
    
    % Initial measure of package position (robots only if they see it)
    x0 = nan(Nrobots,2); % [x_est, y_est]
    distances = zeros(Nrobots, 1);

    for i = 1:Nrobots
        dx = Packages(current_pkg).state(1) - Robots(i).state(1);
        dy = Packages(current_pkg).state(2) - Robots(i).state(2);
        dist = sqrt(dx^2 + dy^2);
        distances(i) = dist;

        if dist <= Robots(i).sr
            % Noisy distance measurement
            noisy_dist = dist + sigma_LiDAR * randn();
            % Assume robot points directly toward the package
            est_x = Robots(i).state_est(1) + noisy_dist * dx/dist;
            est_y = Robots(i).state_est(2) + noisy_dist * dy/dist;
            x0(i,:) = [est_x, est_y];
        end
    end

    % Find robots that can see the package
    idxNonNaN = find(~isnan(x0(:,1)));
    nNonNaN = numel(idxNonNaN);
    disp('Robots that see the package (IDs): '); disp(idxNonNaN);
    disp('Number of robots that see the package: '); disp(nNonNaN);


    %% CONSENSUS ALGORITHM (broadcast one at a time) - only among robots that see the package
    % Preallocate storage for all iterations
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
    alpha_sub = 1/(nNonNaN);
    Q_tot = eye(nNonNaN);

    for t = 1:n_msg

        % choose the transmitter in round-robin among the robots that see the package
        i_local = mod(t-1, nNonNaN) + 1;

        % Qi_sub implements the cyclic broadcast gossip: the selected node broadcasts,
        % others update as x_j <- (1-alpha)*x_j + alpha*x_i, and the broadcaster's
        Qi_sub = zeros(nNonNaN);

        if nNonNaN == 2
            Qi_sub = 1/2 * ones(2);
        else
            for j = 1:nNonNaN
                if j ~= i_local
                    Qi_sub(j,j)      = 1 - alpha_sub;   % keeps part of own value
                    Qi_sub(j,i_local)= alpha_sub;       % receives fraction from i_local
                    Qi_sub(i_local,j)= alpha_sub;       % symmetric term to maintain double-stochasticity
                end
            end

            Qi_sub(i_local,i_local) = 1 - alpha_sub*(nNonNaN-1); % ensures rows sum to 1
        end

        Q_tot = Qi_sub*Q_tot;

        x_next_sub = Qi_sub * xStore_sub(:,:,t);
        xStore_sub(:,:,t+1) = x_next_sub;
        
        % store the evolution of the estimated position for plot
        if nNonNaN >= 1 && pkg_i == 2
            xStore(:,:,t+1) = xStore(:,:,t);
            xStore(idxNonNaN,:,t+1) = x_next_sub;
        end
    end

    % assign to the package its own state estimate
    Packages(current_pkg).state_est = mean(xStore_sub(:,:,end),1);
    % pkg localization error (euclidean norm between real and estimate)
    err_pkg = sqrt( (Packages(current_pkg).state(1) - Packages(current_pkg).state_est(1))^2 + (Packages(current_pkg).state(2) - Packages(current_pkg).state_est(2))^2 );
    pkg_loc_error(current_pkg) = err_pkg; % [m] error for this pkg

    disp(['Real package position (x, y): ', num2str(Packages(current_pkg).state(1)), '  ', num2str(Packages(current_pkg).state(2))]);
    disp(['Est. package position (x, y): ', num2str(xStore_sub(1,1,end)), '  ',num2str(xStore_sub(1,2,end))]);
    disp(['error (distance) [cm]: ', num2str(100*sqrt((Packages(current_pkg).state(1)-xStore_sub(1,1,end))^2+(Packages(current_pkg).state(2)-xStore_sub(1,2,end))^2))])


    tolerance_on_lineup_form_radius = 0.20 + 0.8*err_pkg;

    for j = 1:numel(Robots)
        Robots(j).sr = 12;
    end

    %% ORGANIZE TRANSPORTATION
    % Based on surface area of the package --> n° robots needed
    Packages(current_pkg).s = 20; % [m^2] surface area of the package
    robot_sc = 4; % [m^2] surface capacity of a robot

    disp('--- Organizing package transportation ---');
    RN = ceil(Packages(current_pkg).s / robot_sc);
    disp(['Robots needed to carry package 1: ', num2str(RN)])

    % Which robots will go? The ones that are closer
    % Considering that, after consensous on package position, all robots can share
    % the information regardless of which ones see the package
    [~, sortedIndices] = sort(distances);
    Robots_selected_id = sortedIndices(1:RN);

    for i=1:length(Robots_selected_id)
        Robots(Robots_selected_id(i)).item_id = current_pkg;
    end

    disp('Indices of robots selected to carry the package: ');
    disp(Robots_selected_id);

    % Give the selected robots the proper working state
    for i = 1:length(Robots_selected_id)
        idx = Robots_selected_id(i);
        Robots(idx).working_state = 'r'; % ring working state test
        Robots(idx).target = Packages(current_pkg).state_est(1:2);
    end

    radius = ceil(sqrt(Packages(current_pkg).s / pi)*1.5); % compute radius from Packages(current_pkg).s circular surface
    Packages(current_pkg).r = radius;

    %% TARGET of the package, in the OUTBOUND zone
    outbound_id = randi([1, 3]); % randomly select an outbound zone
    if outbound_id == 1
        zones.outbound = zones.outbound1;
    elseif outbound_id == 2
        zones.outbound = zones.outbound2;
    elseif outbound_id == 3
        zones.outbound = zones.outbound3;
    end
    [x,y] = spawn_target(zones.outbound);
    Packages(current_pkg).target = [x, y];

    %% SIMULATION START


    center_pkg = NaN(Npack, 2);

    [L, areas, masses, centroids, Phi, Wrs_set] = voronoi_lloyd_ring_dyna_unc_reactive(Robots, Robots_voronoi, X, Y, free_mask, sigma_point, sigma_ring, Packages(current_pkg).r, USEREACTIVE);

    % plot
    RN   = numel(Robots_voronoi);

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

    % inbound/outbound polygons
    fill(zones.inbound.polygon(:,1),  zones.inbound.polygon(:,2),  'b', 'FaceAlpha', 0.1, 'EdgeColor','none');
    fill(zones.outbound1.polygon(:,1), zones.outbound1.polygon(:,2), 'r', 'FaceAlpha', 0.1, 'EdgeColor','none');
    fill(zones.outbound2.polygon(:,1), zones.outbound2.polygon(:,2), 'r', 'FaceAlpha', 0.1, 'EdgeColor','none');
    fill(zones.outbound3.polygon(:,1), zones.outbound3.polygon(:,2), 'r', 'FaceAlpha', 0.1, 'EdgeColor','none');

    % plot handles
    hRob  = gobjects(RN,1);
    hCent = gobjects(RN,1);
    hPkg  = gobjects(Npack,1);
    hRing = gobjects(Npack,1);
    hTrgt = gobjects(Npack,1);

    % initialize
    for j = 1:RN
        id = Robots_voronoi(j);
        p  = Robots(id).state(1:2);
        hRob(j)  = plot(p(1), p(2), 'o', 'MarkerSize',6, 'MarkerFaceColor',cmap(j,:), 'Color','k');
        hCent(j) = plot(NaN, NaN, 'r+', 'MarkerSize',10, 'LineWidth',1.2);
    end

    tt   = linspace(0, 2*pi, 361); % ring plot

    xr = Packages(current_pkg).state_est(1) + Packages(current_pkg).r*cos(tt);
    yr = Packages(current_pkg).state_est(2) + Packages(current_pkg).r*sin(tt);
    hRing(current_pkg) = plot(xr, yr, 'k--', 'LineWidth', 1);
    hPkg(current_pkg)  = plot(xr, yr, 'bs', 'MarkerSize', 10, 'MarkerFaceColor','b');
    hTrgt(current_pkg) = plot(Packages(current_pkg).target(1), Packages(current_pkg).target(2), 'bo', 'MarkerSize', 10);


    % refresh plot every iteration
    refresh_every = 1;          % robot/centrois
    refresh_contours_every = 20; % contours of Wrs

clear center
disp('--- STARTING A NEW DELIVERY...');

for k = 1:iter_sim
        % For the auxiliary plot
        if k <= N_iter && pkg_i == 1
            realPos = Robots(1).state;
            realState(:,k) = realPos(:);
            estPos = Robots(1).state_est;
            estState(:,k) = estPos(:);
            trP(k) = trace(Robots(1).P);
        end

        % state estimation
        for j = 1:numel(Robots_voronoi)
            % Measure of distances wrt anchors and orientation
            distances = sqrt(sum((anchors - [Robots(j).state(1), Robots(j).state(2)]).^2, 2));
            distances_noisy = distances + sigma_UWB * randn(Nanchors, 1);
            % trilateration
            [H_tril,z_tril,R_tril] = trilateration(anchors, distances_noisy, sigma_UWB);
            H_tril = [H_tril, zeros(size(H_tril,1),1)];
            % Orientation measure
            z_theta = Robots(j).state(3) + sigma_theta*randn();
            H_theta = [0, 0, 1];
            R_theta = sigma_theta^2;

            % Estimation of position and orientation using EKF
            [Robots(j).state_est, Robots(j).P] = EKF_function(Robots(j).state_est, Robots(j).u, Robots(j).omega, fun2, A, G, H_tril, z_tril, R_tril, H_theta, z_theta, R_theta, Robots(j).P, Q, dt, k);

        end

        % compute Voronoi partitioning + centroids for each robot
        [L, areas, masses, centroids, Phi, Wrs_set] = voronoi_lloyd_ring_dyna_unc_reactive( ...
            Robots, Robots_voronoi, X, Y, free_mask, ...
            sigma_point, sigma_ring, Packages(current_pkg).r, USEREACTIVE);

        % check the robots lineup
        % Get all robots carrying package i
        robots_id_itemId_i = [];  % initialize empty list

        for r = 1:length(Robots)
            % Skip robots with empty or NaN item_id
            if isempty(Robots(r).item_id) || isnan(Robots(r).item_id)
                continue;
            end

            % Check if this robot carries the i-th package
            if Robots(r).item_id == current_pkg
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
            isAligned = isRobotAligned(Robots(id), Packages(current_pkg).state, Packages(current_pkg).r, tolerance_on_lineup_form_radius);
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

        %% FORMATION CONTROL
        for id = robots_id_itemId_i
            if Robots(id).working_state ~= 't'
                break;
            end

            [p_des, center] =  TRANSPORT_control(Robots, robots_id_itemId_i, id, Packages(current_pkg).r, Packages(current_pkg).target, w_form, w_att, w_obs, u_max, dt, Obstacles);
            [u, omega] = ROB_control(Robots(id).state_est, p_des(1,:), ...
                u_sat, omega_sat, Kp_u, Kp_theta, ...
                r_goal, theta_goal);
            Robots(id).u = u;
            Robots(id).omega = omega;
        end
        % compute the formation error
        if Robots(robots_id_itemId_i(1)).working_state == 't' && current_pkg == 1 && form_err_idx <= N_iter
            form_err(form_err_idx) = formation_error(Robots, robots_id_itemId_i, Packages(current_pkg).r);
            form_err_idx = form_err_idx + 1;

        end

        if (exist('center','var') && all(isfinite(center)))
            center_pkg(current_pkg,:) = center(:).'; % package position
            Packages(current_pkg).state_est = center_pkg(current_pkg,:);
        end
        % check delivery
        tol_target = 0.2; % tolerance on distance [m]
        delivered = check_delivered(Packages(current_pkg).state_est, Packages(current_pkg).target, zones.outbound.polygon, tol_target);

        if delivered % CHANGE STATE OF THE ACTIVE ROBOTS
            % change state of Packages(current_pkg)
            Packages(current_pkg).delivered = true;
            % change working state of active robots in the list
            for j = 1:numel(Robots_selected_id)
                Robots(Robots_selected_id(j)).working_state = 'i'; % inbound attraction will guide them to the inbound zone
            end
            % green marker on target
            plot(Packages(current_pkg).target(1), Packages(current_pkg).target(2), 'gp', 'MarkerSize',12, 'MarkerFaceColor','g');
            drawnow;
            break;  % exit current simulation

        end


        % robot control + dynamics
        for j = 1:numel(Robots_voronoi)
            id = Robots_voronoi(j);
            ci = centroids(j,:);

            % control
            if Robots(id).working_state ~= 't'
                [u, omega] = ROB_control(Robots(id).state_est, ci, ...
                    u_sat, omega_sat, Kp_u, Kp_theta, ...
                    r_goal, theta_goal);

                Robots(id).u = u;
                Robots(id).omega = omega;
            end

            % compute dynamics for each robot (update of the real state)
            Robots(id).state = fun2(Robots(id).state(1), Robots(id).state(2), Robots(id).state(3), Robots(id).u, Robots(id).omega, dt);

            % wrap theta
            Robots(id).state(3) = wrap(Robots(id).state(3));
        end

        % plot
        if mod(k, refresh_every) == 0
            perc = 100*k/iter_sim;
            set(hTitle, 'String', sprintf('SIMULATION PROGRESS (wrt max. allowable time): %.1f%%', perc));

            % robots + centroids
            for j = 1:RN
                id = Robots_voronoi(j);
                px = Robots(id).state(1); py = Robots(id).state(2);
                set(hRob(j), 'XData', px, 'YData', py); % robot marker
                cx = centroids(j,1); cy = centroids(j,2);
                set(hCent(j), 'XData', cx, 'YData', cy);
            end

            % package + ring
            cx = Packages(current_pkg).state_est(1); cy = Packages(current_pkg).state_est(2);
            if all(isfinite([cx,cy]))
                set(hPkg(current_pkg), 'XData', cx, 'YData', cy);
                xr = cx + Packages(current_pkg).r*cos(tt);
                yr = cy + Packages(current_pkg).r*sin(tt);
                set(hRing(current_pkg), 'XData', xr, 'YData', yr);
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


if pkg_i == Npack
    disp('--- ALL DELIVERIES CONCLUDED');
end

end



%% Auxiliary plots

%% Evolution of package position during consensus
if PLOT_CONSENSUS
    iterations = 0:n_msg; 
    figure(); clf;
    
    % Subplot X
    subplot(2,1,1); hold on; grid on;
    for k = 1:nNonNaN
        plot(iterations, squeeze(xStore_sub(k,1,:)), 'LineWidth', 1.5);
    end
    yline(Packages(current_pkg).state(1), 'k--', 'LineWidth', 2, 'DisplayName','True X');
    
    xlabel('Iteration (message step)');
    ylabel('Estimate x (coordinate)');
    title('Evolution of X estimates for robots that see the package');
    legendStrings = arrayfun(@(r) sprintf('Robot %d', r), idxNonNaN, 'UniformOutput', false);
    legend(legendStrings, 'Location','best');
    
    % Subplot Y
    subplot(2,1,2); hold on; grid on;
    for k = 1:nNonNaN
        plot(iterations, squeeze(xStore_sub(k,2,:)), 'LineWidth', 1.5);
    end
    % true Y value as reference line
    yline(Packages(current_pkg).state(2), 'k--', 'LineWidth', 2, 'DisplayName','True Y');
    
    xlabel('Iteration (message step)');
    ylabel('Estimate y (coordinate)');
    title('Evolution of Y estimates for robots that see the package');
    legend(legendStrings, 'Location','best');
end

%% State estimation error over time
if PLOT_STATEERROR
    % Compute estimation error
    errorState = realState - estState;
    % Eliminate broken columns
    valid_cols = any(realState ~= 0, 1);
    % Position error (euclidean norm in x,y)
    posErr = sqrt(errorState(1,valid_cols).^2 + errorState(2,valid_cols).^2);

    % Error metrics
    mean_state_pos_error  = mean(posErr); % [m] mean error on position
    mean_state_theta_error = mean(abs(errorState(3,valid_cols)));  % [rad] mean error on orientation

    % Number of iterations
    iters = 1:size(realState,2);
    
    % Plot errors
    figure;
    subplot(2,1,1);
    plot(iters, abs(errorState(1,:)), 'LineWidth', 1); hold on;
    plot(iters, abs(errorState(2,:)), 'LineWidth', 1);
    ylabel('Position Error [m]');
    title('State Estimation Error');
    legend('Error in x','Error in y');
    grid on;
    
    subplot(2,1,2);
    plot(iters, abs(atan2(sin(errorState(3,:)), cos(errorState(3,:)))), 'LineWidth', 1);
    ylabel('Error in \theta [rad]');
    xlabel('Iterations');
    grid on;
end

%% Covariance matrix P trace
if PLOT_TRACEP
    figure;
    plot(iters, trP,'LineWidth', 1);
    ylabel('Covariance');
    xlabel('Iterations')
    title('Trace of covariance matrix P');
    grid on;
end

%% Formation error (average on robots)
if PLOT_FORM_ERR
    figure;
    plot(1:length(form_err), form_err,'LineWidth', 1);
    ylabel('Formation error [m]');
    xlabel('Iterations')
    grid on;
end
%% TRAJECTORIES

figure('Color','w'); hold on; axis equal; box on;
xlabel('x [m]'); ylabel('y [m]');
title('Robot trajectories with working states');

% warehouse/background
fill(map.polygon(:,1), map.polygon(:,2), [0.95 0.95 0.95], ...
    'EdgeColor','k', 'LineWidth',1.5);

% inbound/outbound
fill(zones.inbound.polygon(:,1), zones.inbound.polygon(:,2), ...
    'b', 'FaceAlpha',0.15, 'EdgeColor','b', 'LineWidth',1.0);

fill(zones.outbound1.polygon(:,1), zones.outbound1.polygon(:,2), ...
    'r', 'FaceAlpha',0.10, 'EdgeColor','r', 'LineWidth',1.0);
fill(zones.outbound2.polygon(:,1), zones.outbound2.polygon(:,2), ...
    'r', 'FaceAlpha',0.10, 'EdgeColor','r', 'LineWidth',1.0);
fill(zones.outbound3.polygon(:,1), zones.outbound3.polygon(:,2), ...
    'r', 'FaceAlpha',0.10, 'EdgeColor','r', 'LineWidth',1.0);

% obstacles
for k = 1:numel(Obstacles)
    o  = Obstacles(k);
    ox = o.state(1);  oy = o.state(2);
    switch o.type
        case 'c'
            r = o.l/2;
            th = linspace(0, 2*pi, 60);
            x_o = ox + r*cos(th);
            y_o = oy + r*sin(th);
            fill(x_o, y_o, 'k', 'EdgeColor','k');
        case 's'
            half = o.l/2;
            x_o = [ox-half, ox+half, ox+half, ox-half];
            y_o = [oy-half, oy-half, oy+half, oy+half];
            fill(x_o, y_o, 'k', 'EdgeColor','k');
    end
end

% anchors
plot(anchors(:,1), anchors(:,2), '^', ...
    'MarkerSize',8, 'MarkerFaceColor','y', ...
    'MarkerEdgeColor','k', 'LineWidth',1.2);
% pkg
pkg_pos  = Packages(1).state_est(1:2);
pkg_pos = [cx cy];
pkg_side = sqrt(Packages(1).s);
half_s   = pkg_side/4;
x_pkg = [pkg_pos(1)-half_s, pkg_pos(1)+half_s, ...
pkg_pos(1)+half_s, pkg_pos(1)-half_s];
y_pkg = [pkg_pos(2)-half_s, pkg_pos(2)-half_s, ...
pkg_pos(2)+half_s, pkg_pos(2)+half_s];
hPkg = fill(x_pkg, y_pkg, [0.6 0.6 0.6], ...
'EdgeColor','k', 'LineWidth',1.5);

% colors per working state
col_inbound   = [0.6 0.6 0.6];   % 'i' inbound
col_ring      = [0.85 0.33 0.10];% 'r' lineup-ring
col_transport = [0.47 0.67 0.19];% 't' transport
col_other     = [0 0.4470 0.7410]; % others

for i = 1:Nrobots
    isSelected = exist('Robots_selected_id','var') && ismember(i, Robots_selected_id);

    x_i = traj_x(i,:);
    y_i = traj_y(i,:);

    % inbound
    mask = (state_log(i,:) == 'i');
    if any(mask)
        x_plot = x_i; x_plot(~mask) = NaN;
        y_plot = y_i; y_plot(~mask) = NaN;

        plot(x_plot, y_plot, '-', 'Color', col_inbound, 'LineWidth', 0.8);

        % markers
        idx_valid = find(mask);
        plot(x_i(idx_valid(1)),  y_i(idx_valid(1)),  'o', ...
            'MarkerSize', isSelected*2+4, 'MarkerFaceColor','w','MarkerEdgeColor','k');
        plot(x_i(idx_valid(end)), y_i(idx_valid(end)), 's', ...
            'MarkerSize', isSelected*2+5, 'MarkerFaceColor',col_inbound,'MarkerEdgeColor','k');
    end

    % lineup-ring
    mask = (state_log(i,:) == 'r');
    if any(mask)
        x_plot = x_i; x_plot(~mask) = NaN;
        y_plot = y_i; y_plot(~mask) = NaN;

        lw = isSelected*1.2 + 1.2;
        plot(x_plot, y_plot, '-', 'Color', col_ring, 'LineWidth', lw);

        idx_valid = find(mask);
        plot(x_i(idx_valid(1)),  y_i(idx_valid(1)),  'o', ...
            'MarkerSize', isSelected*2+4, 'MarkerFaceColor','w','MarkerEdgeColor','k');
        plot(x_i(idx_valid(end)), y_i(idx_valid(end)), 's', ...
            'MarkerSize', isSelected*2+5, 'MarkerFaceColor',col_ring,'MarkerEdgeColor','k');
    end

    % transportation
    mask = (state_log(i,:) == 't');
    if any(mask)
        x_plot = x_i; x_plot(~mask) = NaN;
        y_plot = y_i; y_plot(~mask) = NaN;

        lw = isSelected*1.2 + 1.2;
        plot(x_plot, y_plot, '-', 'Color', col_transport, 'LineWidth', lw);

        idx_valid = find(mask);
        plot(x_i(idx_valid(1)),  y_i(idx_valid(1)),  'o', ...
            'MarkerSize', isSelected*2+4, 'MarkerFaceColor','w','MarkerEdgeColor','k');
        plot(x_i(idx_valid(end)), y_i(idx_valid(end)), 's', ...
            'MarkerSize', isSelected*2+5, 'MarkerFaceColor',col_transport,'MarkerEdgeColor','k');
    end

    % other states
    mask = ~(state_log(i,:) == 'i' | state_log(i,:) == 'r' | state_log(i,:) == 't');
    mask = mask & ~(all(isnan(x_i) & isnan(y_i)));
    if any(mask)
        x_plot = x_i; x_plot(~mask) = NaN;
        y_plot = y_i; y_plot(~mask) = NaN;

        plot(x_plot, y_plot, '-', 'Color', col_other, 'LineWidth', 0.8);

        idx_valid = find(mask);
        plot(x_i(idx_valid(1)),  y_i(idx_valid(1)),  'o', ...
            'MarkerSize', isSelected*2+4, 'MarkerFaceColor','w','MarkerEdgeColor','k');
        plot(x_i(idx_valid(end)), y_i(idx_valid(end)), 's', ...
            'MarkerSize', isSelected*2+5, 'MarkerFaceColor',col_other,'MarkerEdgeColor','k');
    end
end


xlim([0 map.W]); ylim([0 map.H]);

h_inb = plot(NaN,NaN,'-','Color',col_inbound,   'LineWidth',1.0);
h_ring_sel = plot(NaN,NaN,'-','Color',col_ring, 'LineWidth',2.2);
h_trans_sel = plot(NaN,NaN,'-','Color',col_transport,'LineWidth',2.2);
h_other = plot(NaN,NaN,'-','Color',col_other,   'LineWidth',1.0);

legend([h_inb, h_ring_sel, h_trans_sel, h_other], ...
    {'Inbound / idle', ...
     'Lineup (robots carrying the package)', ...
     'Transport (robots carrying the package)', ...
     'Other states'}, ...
    'Location','bestoutside');