classdef rob
    properties
        id int8 = NaN % [-] identification number
        x double = NaN % [m] x-coord
        y double = NaN % [m] y-coord
        theta double = NaN % [rad] orientation
        u double = 0.0 % [m/s] forward velocity
        omega double = 0.0 % [rad/s] angular velocity
        loaded logical = false % [-] working state
    end
end