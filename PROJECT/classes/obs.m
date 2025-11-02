classdef obs
    properties
        state  (1,2) double = [NaN NaN] % [x y]
        id int32 = int32(0) % [-] identification number
        type char = 's' % c = circle / 's' = square
        l double = NaN % [m] characteristic dimension of the obstacle (diameter if circle, edge if square)
        facecolor (1,3) double = [0.6 0.6 0.6] % colore riempimento per plot
    end
end