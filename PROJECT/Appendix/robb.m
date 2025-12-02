classdef robb
    properties
        id int8 = NaN; % [-] identification number

        x; % Real state of the robot [x, y, theta]
        x_est; % Estimated state of the robot [x, y, theta]

        y; 

        u; % [m/s] forward velocity
        omega; % [rad/s] angular velocity

        P; % Covariance matrix for the position estimation algorithm

        state; % [-] working state [loaded / free]
        item_id; % current package id that the robot is transporting

        rel_dist = NaN; % relative distance with other (Nrobots-1) robots in formation
        P_rel_dist; % Cov. matrix for the estimate rel_dist

    end
end