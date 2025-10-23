function [u, omega] = conversion_to_unicycle(u_des, theta, k_u, k_omega, u_max, omega_max)

%
%   Inputs:
%       u_des      - desired 2D control vector [ux; uy]
%       theta      - current robot heading (radians)
%       k_u        - gain for linear velocity control
%       k_omega    - gain for angular velocity control
%       u_max      - maximum allowed linear velocity
%       omega_max  - maximum allowed angular velocity
%
%   Outputs:
%       u          - commanded linear velocity (saturated)
%       omega      - commanded angular velocity (saturated)
%

    % Handle zero vector case safely
    if norm(u_des) < 1e-8
        u = 0;
        omega = 0;
        return;
    end

    % Desired heading direction
    theta_des = atan2(u_des(2), u_des(1));

    % Orientation error (wrapped to [-pi, pi])
    e_theta = wrapToPi(theta_des - theta);

    % Compute control commands
    u     = k_u * norm(u_des) * cos(e_theta);
    omega = k_omega * e_theta;

    % Apply saturation
    u     = max(min(u, u_max), -u_max);
    omega = max(min(omega, omega_max), -omega_max);
end