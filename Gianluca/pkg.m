classdef pkg
    properties
        id int8 = NaN % [-] identification number
        PL int8 = NaN % [-] Priority Level
        c; % centroid coordinates [x, y, theta]
        v double = NaN % Volumetric dimension of the package
        picked logical = false % [-] picked/unpicked state
        delivered logical = false % [-] delivered/undelivered state
    end
end