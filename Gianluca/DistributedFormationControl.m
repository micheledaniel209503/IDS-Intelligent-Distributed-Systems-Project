%% Formation Control
% The objective of this script is to make the sworm of robot to follow a
% desired trajectory, maintaining a circle formation around the object to
% carry.
% We assume to know:
% -desired path (can be also point-to-point)
% -robots already in the initial position around the object,
% The robots have to estimate their position through an EKF algorithm that
% updates the predicted estimate with a measure of distance between the
% robots and some fixed anchors (whose positions are known). The
% orientation is treated separately, using a specific sensor.

%% Control
% The control of each robot consists of the superposition of two separate
% control inputs. For simplicity we analise the control input in cartesian
% space: u = [ux, uy] and not the one of the unicycle: u = [v, omega]. In
% particular we define: u = [ux, uy] = u_form + u_target
% where u_form consists in a proportional controller that aims at
% maintaining a predefined distance between robot-i and every other robot.
% These distances are computed in order to maintain the circular formation.
% u_target is simply the proportional input that guide each robot toward
% the desired location.
% Since the robots are executing 2 tasks, we can control the relative
% importance of each task with the weight parameters assigned to them.
% It has to be noted that with this kind of control we don't need to
% specify the exact desired position for each robot.
% NB The initial formation of the robots at t=0 must be sufficiently close
% to the desired one, otherwise the formation control can lead to
% instabilities. For now the initial position of each robot is set with an
% uncertainty of 0.3m, so the controller behave quite robustely.

clc; clear; close all;


% Define the system dynamics
fun = @(x, y, theta, vel, omega, dT) [x + vel * cos(theta) * dT; y + vel * sin(theta) * dT; theta + omega * dT]; 

% Define the system matrix
A = @(x, y, theta, vel, omega, dT) [1, 0, -vel * sin(theta) * dT; 0, 1, vel * cos(theta) * dT; 0, 0, 1];

G = @(x, y, theta, vel, omega, dT)[ dT*cos(theta), 0; dT*sin(theta), 0; 0,  dT ];

%% Simulation Parameters
Ts = 0.01;            % timestep
Tfinal = 60;           % total sim time (s)
K = round(Tfinal/Ts);
form_R = 3;           % formation radius around package centroid
%% Uncertainties
noise_std = 0.2; % [m] std on the distance measures with the anchors
sigma_theta = 0.02; % [RAD] std on the orientation measure 
Q = 0.1*eye(2); % process noise covariance (set to zero if no uncertainty on the model)

%% Building the environment (only inbound zone)
width = 100;
height = 100;
map.W = width;
map.H = height;
map.polygon = [0 0; map.W 0; map.W map.H; 0 map.H];

%% Anchors
% Choose N of anchors
Nanchors = 5;
anchors = zeros(Nanchors,2);
anchors(:,1) = map.W * rand(Nanchors,1);
anchors(:,2) = map.H * rand(Nanchors,1);

%% Initialization of package position (initial centroid)
Package.c(1) = map.W/2;
Package.c(2) = map.H/3;

%% Desired trajectory for the package (figure-eight path using Lissajous curve)
c_center = [Package.c(1); Package.c(2)];
traj_Amp = 35;                % trajectory amplitude
omega_traj = 0.1;             % angular frequency

c_des = zeros(2,K);

for k = 1:K
    t = (k-1)*Ts;

    % Figure-eight trajectory (Lissajous curve) NB only ref. position, no
    % more ref. velocity and orientation
    c_des(:,k) = c_center + [traj_Amp * sin(omega_traj * t);
                             traj_Amp * sin(2 * omega_traj * t) / 2];
end


%% Plot desired trajectory and anchors (static)
fig = 0; %initialization for figure index
fig = fig+1;
figure(fig); clf; hold on; axis equal;
plot([0 map.W map.W 0 0],[0 0 map.H map.H 0],'k-'); % outline
plot(c_des(1,:),c_des(2,:),'b-','LineWidth',1.2);
plot(anchors(:,1),anchors(:,2),'ro','MarkerFaceColor','w','MarkerSize',8)
title('Desired Trajectory for the package centroid');
xlabel('x'); ylabel('y');
legend('Map boundary','Desired centroid path','Anchors');


%% Robots
Nrobots = 6;
Robots(Nrobots,1) = rob();

% Initial formation angles around the package
phi = linspace(0, 2*pi, Nrobots + 1); phi = phi(1:end-1);

%% Initial robot states (x,y,theta)
for i = 1:Nrobots
    Robots(i).id = i;
    Robots(i).x = zeros(3,1);
    Robots(i).x(1) = Package.c(1) + form_R * cos(phi(i)) + 0.3*randn(); % initial X
    Robots(i).x(2) = Package.c(2) + form_R * sin(phi(i)) + 0.3*randn(); % initial Y
    Robots(i).x(3) = 0; % random heading------------------------------------------------------
end

%% Initialisation for the estimation algorithm (recursive least square)
% We can't use EKF for the estimate's initialization (static situation)
for i = 1:Nrobots
    distances = sqrt(sum((anchors - [Robots(i).x(1), Robots(i).x(2)]).^2, 2));
    distances_noisy = distances + noise_std * randn(Nanchors, 1);

    [H,z,C] = trilateration(anchors, distances_noisy, 0.1);
    P = (H'*C^-1*H)^-1;
    x_ls = P*H'*C^-1*z;
    Robots(i).x_est(1) = x_ls(1);
    Robots(i).x_est(2) = x_ls(2);
    Robots(i).x_est(3) = Robots(i).x(3); % initial estimate on initial theta = exact initial theta
  
    P_values = [P, zeros(2, 1); zeros(1, 3)];
    P_values(3,3) = 0.1;
    Robots(i).P = P_values;
    
end

%% Initialisation for the control inputs
% initial zero control input
for i=1:Nrobots
    Robots(i).u = 0;
    Robots(i).omega = 0;
end


%% Controller gains and limits
% Different weights assigned to the different tasks: formation keeping and
% follow trajectory
%% WEIGHTS
w_form  = 0.1; % FORMATION TASK
w_target = 1.0; % TARGET TASK
%% GAINS
k_u   = 3.0;      % linear
k_omega = 2.0;      % angular
u_max = 3.0;         % m/s
omega_max = 1.0;     % rad/s

%% Initialize figure for animation
fig = fig+1;
figure(fig); clf; hold on;axis equal;
plot([0 map.W map.W 0 0],[0 0 map.H map.H 0],'k-'); % outline
plot(anchors(:,1),anchors(:,2),'ro','MarkerFaceColor','w','MarkerSize',8)
h_des_path = plot(NaN,NaN,'b-'); % desired centroid path shown progressively
h_robots = gobjects(Nrobots,1);
for i=1:Nrobots
    h_robots(i) = plot(Robots(i).x(1), Robots(i).x(2), 'bo', 'MarkerFaceColor','w');
end
h_package = plot(Package.c(1), Package.c(2), 'ro', 'MarkerFaceColor','g', 'MarkerSize',8);
title('Robots following desired centroid while maintaining formation');
xlabel('x'); ylabel('y'); grid on;
legend('','Anchors','Desired Path','Robots','Package centroid','Location','bestoutside');

% define vector for actual package traj
act_traj = zeros(2,K);

% initialisation of the trace of the covariance matrix
trace_P = zeros(K);
%% Simulation loop
dt = Ts;
for k = 2:K
    % Current desired centroid and orientation
    % NB We make the hyp that these quantities are known
    cd = c_des(:,k);
    
    % For plotting desired path up to current time
    set(h_des_path,'XData',c_des(1,1:k),'YData',c_des(2,1:k));
    
    % For each robot compute desired formation position and control
    for i = 1:Nrobots

        % Measure of distances wrt anchors and orientation
        distances = sqrt(sum((anchors - [Robots(i).x(1), Robots(i).x(2)]).^2, 2));
        distances_noisy = distances + noise_std * randn(Nanchors, 1);
        % trilateration
        [H_tril,z_tril,R_tril] = trilateration(anchors, distances, noise_std);
        z_theta = Robots(i).x(3) + sigma_theta*randn();
        % Extend H to include theta state and add orientation measurement row
        m = size(H_tril,1);
        H_extended = [H_tril, zeros(m,1)];  
        H_theta = [0, 0, 1];                   
    
        % Stack H and z
        H = [H_extended; H_theta];         
        z = [z_tril; z_theta];

        R = blkdiag(R_tril, sigma_theta^2);



        % Estimation of position and orientation using EKF (based on the last control input)
        [Robots(i).x_est, Robots(i).P] = EKF_function(Robots(i).x_est, Robots(i).u, Robots(i).omega, fun, A, G, z, H, Robots(i).P, Q, R, dt);
        % Update the trace of the covariance matrix
        trace_P(k) = trace_P(k) + trace(Robots(i).P);
    end

            %% CONTROL (all-pair formation):
    for i = 1:Nrobots

        pos_i = Robots(i).x_est(1:2);  % use estimated positions for control

        u_form = [0;0];
        for j = 1:Nrobots
            if j == i
                continue;
            end

            pos_j = Robots(j).x_est(1:2);   % use estimated position of robot j

            diff = pos_i - pos_j;
            dij = norm(diff);            % actual distance between i and j

            % compute circular-formation desired distance (chord length)
            k_sep = abs(i - j);
            k_sep = min(k_sep, Nrobots - k_sep);   % shortest wrap-around separation
            d_des_ij = norm(2 * form_R * sin(pi * k_sep / Nrobots));

            % accumulate formation term
            u_form = u_form - (dij^2 - d_des_ij^2) * diff;
        end

        % target drive toward desired centroid as before
        %u_target = (cd + form_R * [cos(phi(i)),sin(phi(i))]' - pos_i);
        u_target = (cd - pos_i);
        u_des = w_form*u_form + w_target*u_target;


        % conversione in unicycle
        theta_des = atan2(u_des(2), u_des(1));
        e_theta   = wrapToPi(theta_des - Robots(i).x_est(3));
        u     = k_u * norm(u_des) * cos(e_theta);
        omega = k_omega * e_theta;

        % apply saturation
        Robots(i).u     = max(min(u, u_max), -u_max);
        Robots(i).omega = max(min(omega, omega_max), -omega_max);
        
        % Integrate unicycle dynamics (simple Euler)-> real state
        % Current robot state
        x_i = Robots(i).x(1);
        y_i = Robots(i).x(2);
        th_i = Robots(i).x(3);

        x_i = x_i + Robots(i).u * cos(th_i) * dt;
        y_i = y_i + Robots(i).u * sin(th_i) * dt;
        th_i = wrapToPi(th_i + Robots(i).omega * dt);
        
        Robots(i).x = [x_i; y_i; th_i];
    end
    
    % Update package real centroid as the mean of robot positions
    avg_state = mean(reshape([Robots.x],3,[])',1)'; % Avg robots state
    Package.c = avg_state(1:2); % Extract only the x and y components
    act_traj(:,k) = Package.c;
    
    % Update plotting
    for i=1:Nrobots
        % delete(h_robots(i))   % remove old shape
        % h_robots(i) = drawRobot(Robots(i).x(1), Robots(i).x(2), Robots(i).x(3), 0.7, 'b');

        set(h_robots(i),'XData',Robots(i).x(1),'YData',Robots(i).x(2));
    end
    set(h_package,'XData',Package.c(1),'YData',Package.c(2));
    
    drawnow limitrate nocallbacks;
end


fig = fig+1;
figure(fig); clf; hold on;
plot(linspace(0,Tfinal,K),trace_P)
plot(linspace(0,Tfinal,K),trace_P);
ylabel('Trace of Matrix P ');
title('mean trace of Covariance Matrix over Time');
grid on;


%% TO DO
% 1)In this script, to estimate their relative distance, I have supposed that
% each robot have access to the estimated location of every other robot.
% This is reasonable, since we assumed full comunication between robots,
% but we can also use another set of sensors (UWB) in order to make more
% accurate measures.

% 2) Include package dynamics


%% Functions
function t = wrapToPi(t)
    t = atan2(sin(t), cos(t));
end


function h = drawRobot(x, y, theta, L, color)
    % Triangle vertices in robot frame
    verts = [ L   0;       % tip
             -L/2  L/3;    % back-left
             -L/2 -L/3];   % back-right
    
    % Rotation matrix
    R = [cos(theta) -sin(theta);
         sin(theta)  cos(theta)];
    
    % Rotate + translate
    verts = (R * verts')' + [x y];
    
    % Draw filled polygon
    h = fill(verts(:,1), verts(:,2), color, 'EdgeColor','k');
end
