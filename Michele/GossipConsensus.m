% Simple example for linear consensus: Once a packet is in the inbound
% zone, the robots that "see" it, exchange message with the other one to
% estimate the packet position. Main hypotheses:
% - All robots can comunicate with each other
% - Only one robot is able to share information at time
% NB If we consider es. robot 1, at time instant 1, it receives data from
% all other robots, but the other are not able to comunicate except with 1.
% Robot 1 send its meaure to all the others.
% References: Cyclic Broadcast Gossip Algorithm
% 
clc
close all
clear all

toRad = pi / 180; % Conversion factor from degrees to radians
toDeg = 180 / pi; % Conversion factor from radians to degrees
 %% Building the environment (only inbound zone)
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
Npack = 1; % For now we consider just 1 package
Packages(Npack,1) = pkg();
for k = 1:Npack
    Packages(k).id = k;
    [Packages(k).state(1), Packages(k).state(2)] = spawn_pkg(zones.inbound);
    Packages(k).v = 0.5 + 0.2*randn(1); % volumetric dimension of the package
end

%% Robots spawn
% create the robots, placing them in random positions, avoiding the inbound
% zone
Nrobots = 15;
Robots(Nrobots,1) = rob();

%% Dimension factor to set the width of the robots position spawn in x

dimensionFactor = 0.8; % Example factor for robot width adjustment

for k = 1:Nrobots
    Robots(k).id = k;
    valid_pos = false;
    while ~valid_pos
        % random spawn in the map (close to the inbound zone)
        x = map.W/2 * (1-dimensionFactor) + dimensionFactor * map.W * rand();
        y = zones.inbound.H + randn * 0.05 * map.H;

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

% Robot sensing radius
robot_sr = 15; 
% Number of messages exchanged
n_msg = 12 ;

% Measurement noise
sigma = 0.5; % std of distance noise

% Initial measure of package position (robots only if they see it)
x0 = nan(Nrobots,2); % [x_est, y_est]
distances = zeros(Nrobots, 1);

for i = 1:Nrobots
    dx = Packages(1).state(1) - Robots(i).state(1);
    dy = Packages(1).state(2) - Robots(i).state(2);
    dist = sqrt(dx^2 + dy^2);
    distances(i) = dist;

    if dist <= robot_sr
        % Noisy distance measurement -> backproject as rough estimate
        noisy_dist = dist + sigma*randn();
        % Assume robot points directly toward the package (simplification)
        est_x = Robots(i).state(1) + noisy_dist * dx/dist;
        est_y = Robots(i).state(2) + noisy_dist * dy/dist;
        x0(i,:) = [est_x, est_y];
    end
end
%% Find robots that can see the package
idxNonNaN = find(~isnan(x0(:,1)));
nNonNaN = numel(idxNonNaN);
disp('Robots that see the package (IDs): '); disp(idxNonNaN);
disp('Number of robots that see the package: '); disp(nNonNaN);


%% Consensus algorithm (broadcast one at a time) - only among robots that see the package
% Preallocate storage for all iterations (for all robots; NaN for non-seers)
xStore = nan(Nrobots,2,n_msg+1);
xStore(:,:,1) = x0;

% If no robot sees the package or only one sees it, no active consensus
if nNonNaN == 0
    warning('No robot sees the package.');
    return;
% elseif nNonNaN == 1
%     warning('only one robot sees the package.');
%     return;
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
    disp(['Transmitter = ', num2str(i_local)]);
    disp(Qi_sub);
    Q_tot = Qi_sub*Q_tot;

    x_next_sub = Qi_sub * xStore_sub(:,:,t);   % [nNonNaN x 2]
    xStore_sub(:,:,t+1) = x_next_sub;

    % Also update the global storage so we can plot global positions later
    % Non-seers remain unchanged (NaN), participating robots are overwritten
    xStore(:,:,t+1) = xStore(:,:,t);
    xStore(idxNonNaN,:,t+1) = x_next_sub;
end
disp(Q_tot);

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
        cx = Robots(i).state(1) + robot_sr*cos(theta);
        cy = Robots(i).state(2) + robot_sr*sin(theta);
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

%% Organize transportation
% Based on surface area of the package --> n° robots needed
clc;
Packages(1).s = 4; % [m^2] surface area of the package
robot_sc = 1/2; % [m^2] surface capacity of a robot

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

disp('Indices of robots selected to carry the package: '); 
disp(Robots_selected_id);

% Plot: highlight selected robots
figure(1)
selX = S(Robots_selected_id, 1);
selY = S(Robots_selected_id, 2);

hSel = plot(selX, selY, 'ro', 'MarkerFaceColor','g', 'MarkerSize',10); 
legend([hPkg, hEst, hSel], ...
       {'Package (true)', 'Consensus estimates', 'Selected robots'}, ...
       'Location','best');


%% TEST: ROBOT CONTROL
% test: one robot must go from its initial position to pkg 1
% robot params
rob_id = Robots_selected_id(5); % select first robot from the list
active_rob = Robots(rob_id);
active_rob.state(3) = 0; % [rad] robot's initial orientation
state_rob = [active_rob.state];
% disp(['Moving robot id: ', num2str(rob_id), '; robot properties at start: ']);
% disp(active_rob);
u_sat = 3/3.6; % [m/s]
omega_sat = 15*toRad; % [rad/s]
Kp_u     = 0.5;          % gain
Kp_theta = 3.0;          % gain
r_goal   = 0.2;            % [m] distance from target at which to stop
theta_goal = deg2rad(90); % [rad] final robot orientation

% package params
active_pkg = Packages(1);
state_pkg = [Packages(1).state(1); Packages(1).state(2)];
disp('Package position: ');
disp(state_pkg)

% Simulation parameters
dt = 0.05; % [s] time step
T_sim = 70; % [s] sim time
iter_sim = T_sim / dt; % Calculate the number of iterations
disp(['--- Simulation ---']);
disp(['Number of iterations: ', num2str(iter_sim)]);

% robot dynamics
fun = @(state, u, omega) [state(1) + u * cos(state(3)) * dt; state(2) + u * sin(state(3)) * dt; state(3) + omega * dt]; 

% % start sim
% traj = zeros(iter_sim, 3);  % [x, y, theta]
% for k = 1:iter_sim
%     % control
%     [u, omega] = ROB_control(active_rob.state, active_pkg.state, ...
%                              u_sat, omega_sat, ...
%                              Kp_u, Kp_theta, r_goal, theta_goal);
%     % update robot state
%     active_rob.state = fun(active_rob.state, u, omega);
%     traj(k, :) = active_rob.state'; % write on trajectory
% 
%     % end of work condition
%     % d = norm(state_pkg - state_rob(1:2));
%     % e_goal = atan2(sin(theta_goal - state_rob(3)), cos(theta_goal - state_rob(3)));
%     % if (d < r_goal && (theta_goal - state_rob(3) < deg2rad(3)))
%     %     traj = traj(1:k,:);
%     %     disp(['Arrived & aligned in ', num2str(k*dt), ' s']);
%     %     break;
%     % end
% end

% figure(1); hold on;
% 
% % trajectory
% hTraj = plot(traj(:,1), traj(:,2), 'k-', 'LineWidth', 1.5);
% % end point
% hEnd = plot(active_rob.state(1), active_rob.state(2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
% 
% 
% legend([hPkg, hEst, hSel, hTraj, hEnd], ...
%        {'Package (true)', 'Consensus estimates', 'Selected robots', 'Robot trajectory', 'End position'}, ...
%        'Location','best');
% 
% disp('---Simulation results---');
% disp(['Robot state: X = ', num2str(active_rob.state(1)), ' Y = ', num2str(active_rob.state(2)), ' ang = ', num2str(active_rob.state(3)*toDeg)]);

%%
% figure(4)
% plot(1:length(traj), traj(:, 3).*toDeg)
% xlabel('iteration [-]')
% ylabel('ang [deg]')

%% TEST: all selected robots go to package
% let's make all robots travel towards target
% Set all robots' orientation to 0
for j = 1:length(Robots)
    Robots(j).state(3) = 0; % [rad] set robot's orientation to 0
end



% start sim
% trajectories = zeros(iter_sim, 3, length(Robots));  % [x, y, theta] per each  robot
% for k = 1:iter_sim
%     for j = 1:length(Robots_selected_id)
%         idx = Robots_selected_id(j); % select id
%         % control
%         [u, omega] = ROB_control(Robots(idx).state, Packages(1).state, ...
%                                  u_sat, omega_sat, ...
%                                  Kp_u, Kp_theta, r_goal, theta_goal);
%         % update robot state
%         Robots(idx).state = fun(Robots(idx).state, u, omega);
%         trajectories(k, :, j) = Robots(idx).state.'; % write on trajectory
% 
%     end
% end

% figure(1); hold on;
% for j = 1:length(Robots_selected_id)
%     idx = Robots_selected_id(j);
%     % trajectory
%     plot(trajectories(:, 1, j), trajectories(:, 2, j), 'k-', 'LineWidth', 1.5);
%     % end point
%     plot(Robots(idx).state(1), Robots(idx).state(2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'y');
%     legend off;
%     disp(['Robot id: ', num2str(Robots_selected_id(j)),' --- arrival angle: ' , num2str(Robots(Robots_selected_id(j)).state(3)*toDeg), ' [°]'])
% end

% % plot the angle history of the most distant robot
% figure
% plot(1:length(trajectories), trajectories(:, 3, 6).*toDeg)
% xlabel('iteration [-]')
% ylabel('ang [deg]')
% title('History of orientation of the most distance robot')

%% INITIAL FORMATION
% Based on RN number of robots needed, divide the circle centered on the 
% package position in RN slices and select equidistant points in this circle for the initial robot line-up

r_goal = 0.0; % reach the exact point

radius = Packages(1).s*0.75;

lineup_points = ROB_lineup(RN, Packages(1).state, radius); % target points for the initial line-up procedure

% start sim: ROBOT LINE-UP
trajectories = zeros(iter_sim, 3, length(Robots));  % [x, y, theta] per each  robot
for k = 1:iter_sim
    for j = 1:length(Robots_selected_id) % from closest to most distant
        idx = Robots_selected_id(j); % select id
        % control
        [u, omega] = ROB_control(Robots(idx).state, lineup_points(j, :), ...
                                 u_sat, omega_sat, ...
                                 Kp_u, Kp_theta, r_goal, theta_goal);
        % update robot state
        Robots(idx).state = fun(Robots(idx).state, u, omega);
        trajectories(k, :, j) = Robots(idx).state.'; % write on trajectory

    end
end

figure(1); hold on;
for j = 1:length(Robots_selected_id)
    idx = Robots_selected_id(j);
    % trajectory
    plot(trajectories(:, 1, j), trajectories(:, 2, j), 'k-', 'LineWidth', 0.75);
    % end point
    plot(Robots(idx).state(1), Robots(idx).state(2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'y');    
    % target points
    plot(lineup_points(:, 1), lineup_points(:, 2), 'r*', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    legend off;
    disp(['Robot id: ', num2str(Robots_selected_id(j)),' --- arrival angle: ' , num2str(Robots(Robots_selected_id(j)).state(3)*toDeg), ' [°]'])
end

% plot the angle history of the most distant robot (to check)
figure
plot(1:length(trajectories), trajectories(:, 3, length(Robots_selected_id)).*toDeg)
xlabel('iteration [-]')
ylabel('ang [deg]')
title('History of orientation of the most distance robot')