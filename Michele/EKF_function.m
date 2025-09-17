function [state_k_1, P_k_1] = EKF_function(state, v, omega, Dyn_fun, A_fun, G_fun, z, H, P, Q, R, dt)
% EKF_FUNCTION  Extended Kalman Filter for nonlinear state estimation
%
%   [state_k_1, P_k_1] = EKF_function(state, v, omega, Dyn_fun, A_fun, ...
%                                     G_fun, z, H, P, Q, R, dt)
%
%   This function performs one prediction and update cycle of the Extended
%   Kalman Filter (EKF) for a nonlinear system with state vector:
%
%       state = [x; y; theta]
%
%   Inputs:
%       state   - 3x1 current estimated state vector [x; y; theta]
%       v       - linear velocity control input
%       omega   - angular velocity control input
%       Dyn_fun - function handle for nonlinear system dynamics
%                 state_pred = Dyn_fun(x, y, theta, v, omega, dt)
%       A_fun   - function handle for Jacobian of dynamics w.r.t. state
%       G_fun   - function handle for Jacobian of dynamics w.r.t. process noise
%       z       - measurement vector
%       H       - measurement model matrix, dimension mx3 (m: N of anchors)
%       P       - 3x3 state covariance matrix
%       Q       - 3x3 process noise covariance matrix
%       R       - measurement noise covariance matrix
%       dt      - timestep duration
%
%   Outputs:
%       state_k_1 - 3x1 updated state estimate after EKF step
%       P_k_1     - 3x3 updated state covariance matrix

    %% --- Prediction ---
    x = state(1);
    y = state(2);
    th = state(3);

    % Nonlinear motion model
    state_pred = Dyn_fun(x, y, th, v, omega, dt);

    % Jacobian A = df/dx evaluated at current state
    A = A_fun(x, y, th, v, omega, dt);
    % Jacobian G (effect of disturbances on the model)
    G = G_fun(x, y, th, v, omega, dt);
    
    % Predicted covariance
    P_pred = A * P * A' + G * Q * G';

    %% --- Update ---

    y = (z - H * state_pred); % innovation

    % Normalize angular innovation (last element) to [-pi,pi]
    y(end) = atan2(sin(y(end)), cos(y(end)));

    S = H * P_pred * H' + R;
    K = P_pred * H' / S;
    % State update
    state_k_1 = state_pred + K * y;
    % Normalize estimated theta to [-pi,pi]
    state_k_1(3) = atan2(sin(state_k_1(3)), cos(state_k_1(3)));

    % Covariance update
    P_k_1 = (eye(size(P)) - K*H) * P_pred * (eye(size(P)) - K*H)' + K * R * K';

   
end
