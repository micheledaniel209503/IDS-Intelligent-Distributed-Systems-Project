clc
close all
clear all


% --------DEFINITIONS---------
width = 100;
height = 80;


% ---------ROBOTS---------
% state: [x,y,theta]
% inputs: [u,omega]
rob1.state = [4,2,0];
rob2.state = [-4,2,pi/2];
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
    [Packages(k).x, Packages(k).y] = spawn_pkg(zones.inbound);
end

% robot list
Nrob = 3;
Robots(Nrob,1) = rob();
Robots(1).x = 50; Robots(1).y = 30; Robots(1).theta = 0; Robots(1).id = 1; % the rest is all set
Robots(2).x = 60; Robots(2).y = 50; Robots(2).theta = 0 ;Robots(2).id = 2; % the rest is all set
Robots(3).x = 20; Robots(3).y = 20; Robots(3).theta = 0; Robots(2).id = 3; % the rest is all set

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
    plot(Packages(k).x, Packages(k).y, 'ks', 'MarkerSize',6, 'MarkerFaceColor',[1 1 1]);
end
% robots
for k = 1:Nrob
    plot(Robots(k).x, Robots(k).y, 'ko', 'MarkerSize',6, 'MarkerFaceColor',[1 1 1]);
end

xlabel('X [m]')
ylabel('Y [m]')
