clc
close all
clear all

% Simulation parameters
dt = 0.01;                                      % Time step
T_sim = 15;                                     % Simulation time

% --------DEFINITIONS---------
width = 100;
height = 80;


% ---------ROBOTS---------
% state: [x,y,theta]
% inputs: [u,omega]
rob1.state = [4,2,3*pi/2];
rob2.state = [-4,2,3*pi/2];
rob3.state = [0,4,pi];

% ---------MAP-----------
map.W = width;
map.H = height;
map.polygon = [0 0; map.W 0; map.W map.H; 0 map.H];

% inbound/outbound zones
zones.inbound.W = 30;
zones.inbound.H = 10;
zones.inbound.polygon = [
    map.W/2-zones.inbound.W/2 0,
    map.W/2+zones.inbound.W/2 0, 
    map.W/2+zones.inbound.W/2 zones.inbound.H, 
    map.W/2-zones.inbound.W/2 zones.inbound.H
    ]
zones.outbound1.H = 10;
zones.outbound1.W = 5;
zones.outbound1.polygon = [
    0 map.H/2-zones.outbound1.H/2,
    0 map.H/2+zones.outbound1.H/2,
    zones.outbound1.W map.H/2+zones.outbound1.H/2,
    zones.outbound1.W map.H/2-zones.outbound1.H/2
    ]
zones.outbound2.W = 10;
zones.outbound2.H = 5;
zones.outbound2.polygon = [
    map.W/2-zones.outbound2.W/2 map.H,
    map.W/2+zones.outbound2.W/2 map.H,
    map.W/2+zones.outbound2.W/2 map.H-zones.outbound2.H,
    map.W/2-zones.outbound2.W/2 map.H-zones.outbound2.H,
    ]
zones.outbound3.H = 10;
zones.outbound3.W = 5;
zones.outbound3.polygon = [
    map.W map.H/2-zones.outbound3.H/2,
    map.W map.H/2+zones.outbound3.H/2,
    map.W-zones.outbound3.W map.H/2+zones.outbound3.H/2,
    map.W-zones.outbound3.W map.H/2-zones.outbound3.H/2
    ]

% package list
Npack = 10;
Packages(Npack,1) = pkg();
for k = 1:Npack
    Packages(k).id = k;
    Packages(k).PL = Npack+1-k;
    [Packages(k).state(1), Packages(k).state(2)] = spawn_pkg(zones.inbound);
end

% robot list
Nrob = 3;
Robots(Nrob,1) = rob();
Robots(1).state = [50, 30, 0]; Robots(1).id = 1; % the rest is all set
Robots(2).state = [60, 50, 0]; Robots(2).id = 2; % the rest is all set
Robots(3).state = [20, 20, 0]; Robots(3).id = 3; % the rest is all set

% ----------PLOT-----------
figure
hold on; axis equal; grid on
% map
fill(map.polygon(:,1), map.polygon(:,2), [0.9 0.9 0.9], 'EdgeColor', 'k')
% inbound1
fill(zones.inbound.polygon(:,1), zones.inbound.polygon(:,2), [0.8 0.95 1], 'EdgeColor', 'b')
% outbound1
fill(zones.outbound1.polygon(:,1), zones.outbound1.polygon(:,2), [0.8 0.95 1], 'EdgeColor', 'r')
% outbound2
fill(zones.outbound2.polygon(:,1), zones.outbound2.polygon(:,2), [0.8 0.95 1], 'EdgeColor', 'r')
% outbound3
fill(zones.outbound3.polygon(:,1), zones.outbound3.polygon(:,2), [0.8 0.95 1], 'EdgeColor', 'r')
% packages
for k = 1:Npack
    plot(Packages(k).state(1), Packages(k).state(2), 'ks', 'MarkerSize',6, 'MarkerFaceColor',[1 1 1]);
end
% robots
for k = 1:Nrob
    plot(Robots(k).state(1), Robots(k).state(2), 'ko', 'MarkerSize',6, 'MarkerFaceColor',[1 1 1]);
end
xlabel('X [m]')
ylabel('Y [m]')



% Define the system dynamics
fun = @(state, u, omega) [state(1) + u*cos(theta)*dt; state(2) + u*sin(theta)*dt; state(3) + omega*dt];

%% SIMULATION

r_goal = 0.5; % [m] radius in which we consider the robot arrived at pkg
Kp_u = 1.0; % [-] gain for speed P-controller
Kp_theta = 1.0; % [-] gain for angular speed P-controller

% Send robot 1 to take pkg with maximum priority
[max_PL, max_index] = max([Packages.PL]);
target_pkg = Packages(max_index);
active_rob = Robots(1);

% --- Handle per traiettoria e marker robot ---
pathX = active_rob.state(1);
pathY = active_rob.state(2);
hPath = plot(pathX, pathY, 'k-', 'LineWidth', 1.5);        % traiettoria
hRob  = plot(active_rob.state(1), active_rob.state(2), ...
             'ko', 'MarkerSize', 6, 'MarkerFaceColor', [1 1 1]); % marker

% (opz.) freccia di heading
heading_len = 2.0;
hHead = quiver(active_rob.state(1), active_rob.state(2), ...
               heading_len*cos(active_rob.state(3)), ...
               heading_len*sin(active_rob.state(3)), 0, ...
               'LineWidth',1);

for k = 1:200
    [active_rob.u, active_rob.omega] = ROB_control(active_rob.state, target_pkg.state, r_goal, Kp_u, Kp_theta); % compute robot control at current timestamp
    active_rob.state = fun(active_rob.state, active_rob.u, active_rob.omega); % update robot state

    % aggiorna cronologia e grafica
    pathX(end+1) = active_rob.state(1);
    pathY(end+1) = active_rob.state(2);
    set(hPath, 'XData', pathX, 'YData', pathY);
    set(hRob,  'XData', active_rob.state(1), 'YData', active_rob.state(2));
    set(hHead, 'XData', active_rob.state(1), 'YData', active_rob.state(2), ...
               'UData', heading_len*cos(active_rob.state(3)), ...
               'VData', heading_len*sin(active_rob.state(3)));

    drawnow limitrate

    if arrived
        % segna stato se vuoi
        % target_pkg.picked = true;
        break
    end
end




