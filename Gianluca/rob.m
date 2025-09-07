classdef rob
    properties
        id int8 = NaN; % [-] identification number

        x; % Real state of the robot [x, y, theta]
        x_est; % Estimated state of the robot [x, y, theta]

        u; % [m/s] forward velocity
        omega; % [rad/s] angular velocity

        P; % Covariance matrix for the estimation algorithm

        state; % [-] working state [loaded / free]
        item_id; % current package id that the robot is transporting

    end
end