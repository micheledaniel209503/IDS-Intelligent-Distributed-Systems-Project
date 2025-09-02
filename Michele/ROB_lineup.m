function points = ROB_lineup(RN, trg_pos, radius)
% DESCRIPTION
% computes position of RN equidistant points on a circumference of radius
% "radius", then rigidly translate each point by trg_pos to surround
% trg_pos
% NOTE: trg_pos = pkg.state for initial line-up

if RN == 1 % only 1 rob needed
    points = trg_pos;
    return
end

points = zeros(RN, 2);
angle = 2*pi/RN; % angle increment
phi0 = pi/2; % initial shift (to put rob(1) at the top)

for i = 1:RN
    points(i,1) = trg_pos(1) + radius*cos((i-1)*angle + phi0); % x
    points(i,2) = trg_pos(2) + radius*sin((i-1)*angle + phi0); % y
end


end