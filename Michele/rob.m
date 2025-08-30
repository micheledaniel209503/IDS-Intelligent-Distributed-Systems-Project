classdef rob
    properties
        state  (1,3) double = [NaN NaN NaN] % [x y theta]
        id     (1,1) int32  = int32(0)      % [-] identification number, 0 = unassigned
        u      (1,1) double = 0.0           % [m/s]
        omega  (1,1) double = 0.0           % [rad/s]
        loaded (1,1) logical = false        % [-] working state
    end
end