classdef pkg
    properties
        state  (1,2) double = [NaN NaN] % [x y]
        id int32 = int32(0) % [-] identification number
        PL int32 = NaN % [-] Priority Level
        s double = NaN % [m^2] Surface dimension of the package
        v double = NaN % [m^3] Volumetric dimension of the package
        picked logical = false % [-] picked/unpicked state
        delivered logical = false % [-] delivered/undelivered state
    end
end