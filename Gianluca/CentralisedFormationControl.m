%% Formation Control
% The objective of this script is to make the sworm of robot to follow a
% desired trajectory, maintaining a circle formation around the object to
% carry.
% We assume to know:
% -temporal law of the trajectory
% -robots already in the initial position perfectly around the object,
% pointing in the direction of the trajectory to follow.
% The robots have to estimate their position through an EKF algorithm that
% updates the predicted estimate with a measure of distance between the
% robots and some fixed anchors (whose positions are known). The
% orientation is treated separately, using a specific sensor.
clc; clear all; close all;

%% Simulation Parameters
Ts = 0.05;            % timestep
Tfinal = 70;           % total sim time (s)
K = round(Tfinal/Ts);
form_R = 3;           % formation radius around package centroid
%% Uncertainties
noise_std = 0.3; % [m] std on the distance measures with the anchors
sigma_theta = 0.02; % [RAD] std on the orientation measure 
Q = 0.1*eye(2); % process noise covariance (set to zero if no uncertainty on the model)

%% Building the environment (only inbound zone)
width = 100;
height = 100;
map.W = width;
map.H = height;
map.polygon = [0 0; map.W 0; map.W map.H; 0 map.H];

%% Anchors
Nanchors = 4;
anchors = zeros(Nanchors,2);
anchors(:,1) = map.W * rand(Nanchors,1);
anchors(:,2) = map.H * rand(Nanchors,1);

%% Initialization of package position (initial centroid)
Package.c(1) = map.W/1.5;
Package.c(2) = map.H/3;

%% Desired trajectory for the package (figure-eight path using Lissajous curve)
c_center = [Package.c(1); Package.c(2)];
traj_Amp = 25;                % trajectory amplitude
omega_traj = 0.1;             % angular frequency

for k = 1:K
    t = (k-1)*Ts;

    % Figure-eight trajectory (Lissajous curve)
    c_des(:,k) = c_center + [traj_Amp * sin(omega_traj * t);
                             traj_Amp * sin(2 * omega_traj * t) / 2];

    % Velocity (time derivative)
    c_des_dot(:,k) = [traj_Amp * omega_traj * cos(omega_traj * t);
                      traj_Amp * omega_traj * cos(2 * omega_traj * t)];

    % Orientation angle (continuous theta)
    if atan2(c_des_dot(2,k), c_des_dot(1,k)) >= 0
        theta_des(k) = atan2(c_des_dot(2,k), c_des_dot(1,k));
    else
        theta_des(k) = atan2(c_des_dot(2,k), c_des_dot(1,k)) + 2*pi;
    end
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

% fig = fig+1;
% figure(fig);
% plot(linspace(1,Tfinal,length(theta_des)), theta_des(:)*180/pi);


%% Robots
Nrobots = 5;
Robots(Nrobots,1) = rob();
% Initial formation angles around the package
phi = linspace(0, 2*pi, Nrobots + 1); phi = phi(1:end-1);

%% Initial robot states (x,y,theta)
for i = 1:Nrobots
    Robots(i).id = i;
    Robots(i).x = zeros(3,1);
    Robots(i).x(1) = Package.c(1) + form_R * cos(phi(i)) + 0.2*randn(); % X
    Robots(i).x(2) = Package.c(2) + form_R * sin(phi(i)) + 0.2*randn(); % Y
    Robots(i).x(3) = theta_des(1); % heading
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
k_rho   = 2.0;      % linear
k_alpha = 1.0;      % bearing
k_beta  = -1.0;     % final heading (must be negative)
v_max = 5.0;         % m/s
omega_max = 1.0;     % rad/s

%% Initialize figure for animation
fig = fig+1;
figure(fig); clf; hold on; axis equal;
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

%% Simulation loop
dt = Ts;
for k = 2:K
    % Current desired centroid and orientation
    % NB We make the hyp that these quantities are known
    cd = c_des(:,k);
    thd = theta_des(k);
    
    % For plotting desired path up to current time
    set(h_des_path,'XData',c_des(1,1:k),'YData',c_des(2,1:k));
    
    % For each robot compute desired formation position and control
    for i = 1:Nrobots
        % desired relative angle for this robot (rotate formation with theta_des)
        desired_theta = thd; %-> desired orientation of each robot is the same orientation of the package: to be revised!!
        desired_pos = cd + form_R * [cos(desired_theta + phi(i) - pi/2); sin(desired_theta + phi(i) - pi/2)];
        
        % Current robot state
        x_i = Robots(i).x(1);
        y_i = Robots(i).x(2);
        th_i = Robots(i).x(3);

        % Measure of distances wrt anchors and orientation
        distances = sqrt(sum((anchors - [Robots(i).x(1), Robots(i).x(2)]).^2, 2));
        distances_noisy = distances + noise_std * randn(Nanchors, 1);
        z_theta = th_i + sigma_theta*randn();

        % Estimation of position and orientation using EKF (based on the last control input)
        [Robots(i).x_est, Robots(i).P] = unicycle_ekf(Robots(i).x_est, Robots(i).P, Robots(i).u, Robots(i).omega, Ts, anchors, distances_noisy, noise_std, Q, sigma_theta, z_theta);
         % Compute (proportional) unicycle control to reach desired_pos and orientation aligned with theta_des
         % Choose type of controller:
         % 1) GPT proportional control
         % 2) Simple proportional controller

        % [Robots(i).u, Robots(i).omega] = unicycle_controller1(Robots(i).x_est, desired_pos, desired_theta, k_rho, k_alpha,k_beta, v_max, omega_max);
        [Robots(i).u, Robots(i).omega] = unicycle_controller1(Robots(i).x_est, cd, desired_theta, k_rho, k_alpha,k_beta, v_max, omega_max);
        
        % Integrate unicycle dynamics (simple Euler)
        x_i = x_i + Robots(i).u * cos(th_i) * dt;
        y_i = y_i + Robots(i).u * sin(th_i) * dt;
        th_i = wrapToPi(th_i + Robots(i).omega * dt);
        
        Robots(i).x = [x_i; y_i; th_i];
    end
    
    % Update package centroid as the mean of robot positions
    avg_state = mean(reshape([Robots.x],3,[])',1)'; % Avg robots state
    Package.c = avg_state(1:2); % Extract only the x and y components
    act_traj(:,k) = Package.c;
    
    % Update plotting
    for i=1:Nrobots
        set(h_robots(i),'XData',Robots(i).x(1),'YData',Robots(i).x(2));
    end
    set(h_package,'XData',Package.c(1),'YData',Package.c(2));
    
    drawnow limitrate nocallbacks;
end

% err_traj = sqrt((act_traj(1,:)-c_des(1,:)).^2+(act_traj(1,:)-c_des(2,:)).^2);
% fig = fig+1;
% figure(fig); hold on;
% plot(linspace(1,400,K),err_traj,'-')




%%  Controller function (proportional unicycle controller) ---

function [v, omega] = unicycle_controller1(state, p_des, th_des, ...
                                               k_rho, k_alpha, k_beta, ...
                                               v_max, omega_max)
% state : [x; y; theta]
% p_des : desired position [x; y]
% th_des: desired final heading
% gains : k_rho>0, k_alpha>0, k_beta<0 (with k_alpha > k_rho recommended)

    x = state(1); y = state(2); th = state(3);
    dx = p_des(1) - x;
    dy = p_des(2) - y;

    rho   = hypot(dx, dy);
    alpha = atan2(dy, dx) - th;         % bearing error
    alpha = atan2(sin(alpha), cos(alpha));
    beta  = th_des - th - alpha;        % heading-at-goal error
    beta  = atan2(sin(beta), cos(beta));

    % Core law (robust, Lyapunov-based)
    v     = k_rho  * rho * cos(alpha);
    omega = k_alpha* alpha + k_beta * beta;

    % Optional: slow down when turning hard (helps with noise / sharp corners)
    v = v / (1 + 0.5*abs(alpha));

    % Saturations
    v     = max(-v_max, min(v_max, v));
    omega = max(-omega_max, min(omega_max, omega));
end

function [v, omega] = unicycle_controller2(state, p_des, th_des, ...
                                               k_rho, k_alpha, k_beta, ...
                                               v_max, omega_max)
% state : [x; y; theta]
% p_des : desired position [x; y]
% th_des: desired final heading
% gains : k_rho>0, k_alpha>0 (with k_alpha > k_rho recommended)

    x = state(1); y = state(2); th = state(3);
    dx = p_des(1) - x;
    dy = p_des(2) - y;

    rho   = norm([dx, dy]);
    alpha = atan2(dy, dx) - th;         % bearing error
    alpha = atan2(sin(alpha), cos(alpha));

    % Core law (simple proportional)
    v     = k_rho  * rho;
    omega = k_alpha* alpha;
end

%% Function for angle conversion
function t = wrapToPi(t)
    t = atan2(sin(t), cos(t));
end

%% Extended Kalman filter
function [state_est, P_est] = unicycle_ekf( ...
    state, P, v, omega, dt, anchors, distances, noise_std, Q, sigma_theta, z_theta)
% state = [x; y; theta]
% P     = 3x3 covariance
% v, omega = control inputs (omega used for prediction only)
% dt = timestep
% anchors = Nx2 anchor positions
% distances = Nx1 measured ranges to anchors
% noise_std = std of individual range measurements (used by trilateration)
% Q = 3x3 process noise covariance
% sigma_theta = std of the orientation sensor (rad)
% z_theta = measured orientation (rad) from the orientation sensor

    % --- Prediction ---
    x = state(1);
    y = state(2);
    th = state(3);

    % Nonlinear motion model
    x_pred  = x + v*dt*cos(th);
    y_pred  = y + v*dt*sin(th);
    th_pred = th + omega*dt;

    state_pred = [x_pred; y_pred; th_pred];

    % Jacobian A = df/dx evaluated at current state
    A = [1, 0, -v*dt*sin(th);
         0, 1,  v*dt*cos(th);
         0, 0,  1];
    % Jacobian G (effect of disturbances on the model)
    G = [ dt*cos(th), 0; 
          dt*sin(th), 0;
                0      , dt ];
    % Predicted covariance
    P_pred = A * P * A' + G * Q * G';

    % --- Measurement (trilateration for x,y) ---
    % trilateration returns H_xy (m-1 x 2), z_tril ((m-1)x1), and C ((m-1)x(m-1))
    [H_xy, z_tril, C] = trilateration(anchors, distances, noise_std);

    % Number of pseudo-range equations
    m = size(H_xy,1);

    % Extend H to include theta state and add orientation measurement row
    H_pos_extended = [H_xy, zeros(m,1)];   % m x 3
    H_theta = [0, 0, 1];                   % 1 x 3

    % Stack H and z
    H = [H_pos_extended; H_theta];         % (m+1) x 3
    z_stack = [z_tril; z_theta];           % (m+1) x 1

    % Predicted measurement (h) for the stack:
    h_pos = H_xy * state_pred(1:2);        % m x 1
    h_theta = state_pred(3);               % 1 x 1
    h_stack = [h_pos; h_theta];            % (m+1) x 1

    % Innovation (measurement residual)
    y_innov = z_stack - h_stack;

    % Normalize angular innovation (last element) to [-pi,pi]
    y_innov(end) = atan2(sin(y_innov(end)), cos(y_innov(end)));

    % Measurement covariance (stacked): trilateration C and orientation sigma
    C_new = blkdiag(C, sigma_theta^2);         % (m+1) x (m+1)

    % Kalman gain
    S = H * P_pred * H' + C_new;
    K = P_pred * H' / S;

    % State update
    dx = K * y_innov;
    state_est = state_pred + dx;

    % Normalize estimated theta to [-pi,pi]
    state_est(3) = atan2(sin(state_est(3)), cos(state_est(3)));

    % Covariance update (Joseph form can be used for numerical stability)
    I = eye(size(P));
    P_est = (I - K*H) * P_pred * (I - K*H)' + K * C_new * K';
end

%% trilateration function 
function [H,z,C] = trilateration(anchors, distances, noise_std)
    % Number of anchors
    n = size(anchors, 1);
    
    % Initialize matrices
    H = zeros(n-1, 2);
    z = zeros(n-1, 1);
    C = zeros(n-1);
    
    % Iterate over all anchors
    for i = 1:n-1
        % Fill the matrices
        H(i, :) = 2*[anchors(i+1, 1) - anchors(i, 1), anchors(i+1, 2) - anchors(i, 2)];
        z(i) = - distances(i+1)^2  + distances(i)^2 + anchors(i+1, 1)^2 - anchors(i, 1)^2 + anchors(i+1, 2)^2 - anchors(i, 2)^2;
        % Fill the covariance matrix
        if i == 1
            C(i,i) = 4 * noise_std^2 * (distances(i+1)^2 + distances(i)^2);
            if n > 2
                C(i,i+1) = -4 * noise_std^2 * distances(i+1)^2;
            end
        elseif i < n-1
            C(i,i-1) = -4 * noise_std^2 * distances(i)^2;
            C(i,i) = 4 * noise_std^2 * (distances(i+1)^2 + distances(i)^2);
            C(i,i+1) = -4 * noise_std^2 * distances(i+1)^2;
        else
            C(i,i-1) = -4 * noise_std^2 * distances(i)^2;
            C(i,i) = 4 * noise_std^2 * (distances(i+1)^2 + distances(i)^2);
        end
    end
end
