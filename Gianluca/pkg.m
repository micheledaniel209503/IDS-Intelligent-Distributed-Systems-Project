classdef pkg
    properties
        id int8 = NaN % [-] identification number
        PL int8 = NaN % [-] Priority Level
        x double = NaN % [m] x-coord
        y double = NaN % [m] y-coord
        v double = NaN % Volumetric dimension of the package
        picked logical = false % [-] picked/unpicked state
        delivered logical = false % [-] delivered/undelivered state
    end
end