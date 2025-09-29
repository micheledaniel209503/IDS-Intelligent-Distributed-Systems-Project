%% DEMO SCRIPT: trajectory planning with obstacle avoidance (simplified)
clear; clc; close all;
%% Map setup
map.W = 100;
map.H = 100;
map.polygon = [0 0; map.W 0; map.W map.H; 0 map.H];

%% Generate random obstacles
Nobs = 10;
obstacle_dim = 5.0;   % [m] size of obstacles
minDist = 10;          % minimum clearance between obstacles
obstacles = spawn_obstacles(Nobs, map, obstacle_dim, minDist);

%% Trajectory parameters
startPoint   = [5 5];
goalPoint    = [map.W/1.05 map.H/1.05];
numPoints    = 100;   % number of samples along trajectory
minClearance = 4.0;   % [m] safety distance from obstacles

%% Plan trajectory
trajectory = planTrajectory([map.W map.H], obstacles, startPoint, goalPoint, numPoints, minClearance);

%% Plot map
figure(1); clf; hold on; axis equal;
xlim([0 map.W]); ylim([0 map.H]);
plot([0 map.W map.W 0 0],[0 0 map.H map.H 0],'k-'); % outline

% Obstacles
for k = 1:numel(obstacles)
    obs_k = obstacles(k);
    switch obs_k.type
        case 'c' % circle
            r = obs_k.l/2;
            ang = linspace(0,2*pi,50);
            xc = obs_k.state(1) + r*cos(ang);
            yc = obs_k.state(2) + r*sin(ang);
            fill(xc, yc, obs_k.facecolor, 'EdgeColor','k');
        case 's' % square
            halfSide = obs_k.l/2;
            x = obs_k.state(1) + halfSide*[-1 -1 1 1];
            y = obs_k.state(2) + halfSide*[-1 1 1 -1];
            fill(x, y, obs_k.facecolor, 'EdgeColor','k');
    end
end

% Trajectory
plot(trajectory(:,1), trajectory(:,2), 'b-', 'LineWidth', 2);
plot(startPoint(1), startPoint(2), 'go', 'MarkerSize',10, 'MarkerFaceColor','g'); % start
plot(goalPoint(1), goalPoint(2), 'ro', 'MarkerSize',10, 'MarkerFaceColor','r'); % goal

title('Trajectory Planning with Obstacles (Penalty-based Optimization)');
xlabel('X [m]'); ylabel('Y [m]');
grid on;


