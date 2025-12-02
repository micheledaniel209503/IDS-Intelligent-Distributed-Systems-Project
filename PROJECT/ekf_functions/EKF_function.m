function [state_k_1, P_k_1] = EKF_function(state, v, omega, Dyn_fun, A_fun, G_fun, ...
                                           H_xy, z_xy, R_xy, H_theta, z_theta, R_theta, ...
                                           P, Q, dt, k)
% Inputs:
%   state    - [3x1] current estimated state vector [x; y; theta]
%   v        - linear velocity input
%   omega    - angular velocity input
%   Dyn_fun  - function handle for nonlinear system dynamics:
%   A_fun    - function handle for Jacobian of dynamics w.r.t. state
%   G_fun    - function handle for Jacobian of dynamics w.r.t. process noise
%   H_xy     - measurement matrix for position
%   z_xy     - position measurement vector
%   R_xy     - measurement noise covariance for position
%   H_theta  - measurement matrix for orientation (1 x 3)
%   z_theta  - scalar orientation measurement
%   R_theta  - scalar orientation measurement noise variance
%   P        - [3x3] state covariance matrix
%   Q        - [2x2] process noise covariance matrix
%   dt       - time step
%   k        - current iteration number


    %% --- Update periods ---
    f_xy = 1; % Frequency for position updates
    f_theta = 5; % Frequency for orientation updates

    N_xy    = ceil(1 / (f_xy * dt));
    N_theta = ceil(1 / (f_theta * dt));

    %% --- Prediction step ---
    x = state(1);
    y = state(2);
    th = state(3);

    % Nonlinear motion model
    state_pred = Dyn_fun(x, y, th, v, omega, dt);

    % Jacobians
    A = A_fun(x, y, th, v, omega, dt);
    G = G_fun(x, y, th, v, omega, dt);

    % Predicted covariance
    P_pred = A * P * A' + G * Q * G';

    % Initialize for update phase
    state_k_1 = state_pred;
    P_k_1 = P_pred;

    %% Update position only every N_xy iterations
    if mod(k, N_xy) == 0

        innov_xy = z_xy - H_xy * state_k_1;          % innovation
        S_xy = H_xy * P_k_1 * H_xy' + R_xy;          % innovation covariance
        K_xy = P_k_1 * H_xy' * S_xy^-1;                 % Kalman gain

        state_k_1 = state_k_1 + K_xy * innov_xy;     % state correction
        state_k_1(3) = atan2(sin(state_k_1(3)), cos(state_k_1(3))); % normalize angle

        I = eye(3);
        P_k_1 = (I - K_xy * H_xy) * P_k_1;
    end

    %% Update position only every N_theta iterations
    if mod(k, N_theta) == 0
        innov_th = z_theta - H_theta * state_k_1;    % innovation
        innov_th = atan2(sin(innov_th), cos(innov_th)); % normalization

        S_th = H_theta * P_k_1 * H_theta' + R_theta; % innovation covariance
        K_th = (P_k_1 * H_theta') * S_th^-1;            % Kalman gain

        state_k_1 = state_k_1 + K_th * innov_th;     % state correction
        state_k_1(3) = atan2(sin(state_k_1(3)), cos(state_k_1(3))); % normalize angle

        I = eye(3);
        P_k_1 = (I - K_th * H_theta) * P_k_1;
    end

end