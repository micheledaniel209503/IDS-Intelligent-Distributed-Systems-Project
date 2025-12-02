classdef pkgg
    properties
        id int8 = NaN % [-] identification number
        PL int8 = NaN % [-] Priority Level
        x double = NaN
        y double = NaN
        c = zeros(2,1); % centroid coordinates [x, y, theta]
        v double = NaN % Volumetric dimension of the package
        picked logical = false % [-] picked/unpicked state
        delivered logical = false % [-] delivered/undelivered state
    end
end