classdef rob
    properties
        state  (1,3) double = [NaN NaN NaN] % [x y theta] real state
        state_est (1,3) double = [NaN NaN NaN] % [x y theta] state estimate
        id     (1,1) int32  = int32(0)      % [-] identification number, 0 = unassigned
        u      (1,1) double = 0.0           % [m/s]
        omega  (1,1) double = 0.0           % [rad/s]
        working_state (1,1) char = 'i'      % [-] working state: 'f' = free, 'l' = lineup, 'r' = ring, 't' = transportation, 'i' = inbound attraction
        target (1,2) double = [NaN NaN]     % [x y] target position
        P; % covariance matrix of the estimation algorithm
        item_id; % current package id that the robot is transporting
        sr (1,1) double = 18; % [m] robot sensing radius 
    end
end