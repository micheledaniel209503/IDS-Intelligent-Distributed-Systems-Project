classdef pkg
    properties
        state  (1,2) double = [NaN NaN] % [x y]
        state_est  (1,2) double = [NaN NaN] % [x y]
        id int32 = int32(0) % [-] identification number
        r double = NaN % [m] Radius of the package
        s double = NaN % [m^2] Surface dimension of the package
        delivered logical = false % [-] delivered/undelivered state
        target (1,2) double = [NaN NaN] % [x y]
    end
end