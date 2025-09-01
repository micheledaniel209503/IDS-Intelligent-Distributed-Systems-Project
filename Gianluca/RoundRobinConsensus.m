% Simple example for linear consensus: Once a packet is in the inbound
% zone, the robots that "see" it, exchange message with the other one to
% estimate the packet position. Main hypotheses:
% - All robots can comunicate with each other
% - Only ONE robot is able to share information at time with another SINGLE
% robot --> ROUND ROBIN ALGORITHM
% NB If we consider es. robot 1, at time instant 1, it sends data to
% all other robots, and the other are not able to comunicate except with 1.
% Robot 1 send its measure to all the others, while it doesn't receive anything.
% References: Round Robin Algorithm
% 
clc
close all
clear all

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
    [Packages(k).x, Packages(k).y] = spawn_pkg(zones.inbound);
    Packages(k).v = 0.5 + 0.2*randn(1); % volumetric dimension of the package
end

%% Robots spawn
% create the robots, placing them in random positions, avoiding the inbound
% zone
Nrobots = 15;
Robots(Nrobots,1) = rob();

%% Dimension factor to set the width of the robots position spawn in x

dimensionFactor = 0.5; % Example factor for robot width adjustment

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
            Robots(k).x = x;
            Robots(k).y = y;
            valid_pos = true;
        end
    end
end

% Robot sensing radius
robot_sr = 15; 

% Measurement noise
sigma = 0.5; % std of distance noise

% Initial measure of package position (robots only if they see it)
x0 = nan(Nrobots,2); % [x_est, y_est]
for i = 1:Nrobots
    dx = Packages(1).x - Robots(i).x;
    dy = Packages(1).y - Robots(i).y;
    dist = sqrt(dx^2 + dy^2);

    if dist <= robot_sr
        % Noisy distance measurement -> backproject as rough estimate
        noisy_dist = dist + sigma*randn();
        % Assume robot points directly toward the package (simplification)
        est_x = Robots(i).x + noisy_dist * dx/dist;
        est_y = Robots(i).y + noisy_dist * dy/dist;
        x0(i,:) = [est_x, est_y];
    end
end
%% Find robots that can see the package
idxNonNaN = find(~isnan(x0(:,1)));
nNonNaN = numel(idxNonNaN);
disp('Robots that see the package (IDs): '); disp(idxNonNaN);
disp('Number of robots that see the package: '); disp(nNonNaN);
%% !!!!IMPORTANT!!!
% The number of messages exchanges is NOT a free choice. In fact the
% total transition matrix is doubly stochastic only if the sequence Q1*...Q
% is repeated an INTEGER number of times. es. Qn*...*Q1*Qn IS NOT DOUBLY
% STOCHASTIC!!!

% In this case a n_msg of 3 times the number of 'active' robots is chosen
n_msg = 3*nNonNaN;


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

%% Use the function to build the Qn,Qn-1,....Q1 transition matrices
Q_tot = eye(nNonNaN); % Initialize total transition matrix as identity
[~,Qlist,~] = RoundRobinQmatrix(nNonNaN);

for t = 1:n_msg
    
    % choose the transmitter in round-robin among the robots that see the package
    % i_local is a LOCAL index (1..nNonNaN), not the global robot id
    i_local = mod(t-1, nNonNaN) + 1;

    % Debug: show which local index transmits and the Qi_sub matrix
    disp(['Transmitter = ', num2str(i_local)]);
    disp(Qlist{i_local});
    Q_tot = Qlist{i_local}*Q_tot;

    x_next_sub = Qlist{i_local} * xStore_sub(:,:,t);   % [nNonNaN x 2]
    xStore_sub(:,:,t+1) = x_next_sub;

    % Also update the global storage so we can plot global positions later
    % Non-seers remain unchanged (NaN), participating robots are overwritten
    xStore(:,:,t+1) = xStore(:,:,t);
    xStore(idxNonNaN,:,t+1) = x_next_sub;
end
disp(Q_tot);

disp(['Real package position (x, y): ', num2str(Packages(1).x), '  ', num2str(Packages(1).y)]);
disp(['Est. package position (x, y): ', num2str(xStore_sub(1,1,end)), '  ',num2str(xStore_sub(1,2,end))]);
disp(['error (distance) [cm]: ', num2str(100*sqrt((Packages(1).x-xStore_sub(1,1,end))^2+(Packages(1).y-xStore_sub(1,2,end))^2))])



%% Plot results on the map 
figure(1); clf; hold on; axis equal;
% Map boundary
plot([0 map.W map.W 0 0],[0 0 map.H map.H 0],'k-'); % outline

% Inbound zone
fill(zones.inbound.polygon(:,1),zones.inbound.polygon(:,2), ...
     [0.9 0.9 1],'EdgeColor','b');

% Robots
plot([Robots.x],[Robots.y],'ko','MarkerFaceColor','w')
%text([Robots.x]+1,[Robots.y]+1,string(1:Nrobots))

% Package true position
hPkg = plot(Packages(1).x,Packages(1).y,'rs','MarkerFaceColor','r','MarkerSize',10);

% Final consensus estimates (only those who participated)
hEst = plot(xStore(:,1,end),xStore(:,2,end),'bo','MarkerFaceColor','b');

% --- Highlight robots that sense the package ---
theta = linspace(0,2*pi,100);
for i = idxNonNaN(:,1)'
        cx = Robots(i).x + robot_sr*cos(theta);
        cy = Robots(i).y + robot_sr*sin(theta);
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
    yline(Packages(1).x, 'k--', 'LineWidth', 2, 'DisplayName','True X');

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
    yline(Packages(1).y, 'k--', 'LineWidth', 2, 'DisplayName','True Y');

    xlabel('Iteration (message step)');
    ylabel('Estimate y (coordinate)');
    title('Evolution of Y estimates for robots that see the package');
    legend([legendStrings], 'Location','best');
end
