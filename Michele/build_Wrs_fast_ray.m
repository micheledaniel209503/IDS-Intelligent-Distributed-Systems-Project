function Wrs = build_Wrs_fast_ray(pr, Vmask, obs, X, Y, rs, ds)
% BUILD_WRS_FAST_RAY  W_rs = (S ∩ V) ∩ B_rs(pr) via ray casting
% pr     = [x_r,y_r]
% Vmask  = Voronoi cell of robot
% obs    = ~free_mask (true = obstacle)
% X,Y    = meshgrid
% rs     = sensing radius (Inf or [] = no limit --> classic voronoi)
% ds     = step size along radius

    if ~any(Vmask(:)), Wrs = false(size(Vmask)); return; end % pretty unlikely
    
    % how many rotations
    %ntheta = 128;
    ntheta = 96;

    % precompute
    obsNum = double(obs);
    [m,n]  = size(X);
    xv = X(1,:);
    yv = Y(:,1);

    visible = false(m,n); % initialization

    % angles (-180 --> 180 deg)
    ang = linspace(-pi, pi, ntheta+1); ang(end) = [];

    for t = 1:numel(ang) % loop through each angle
        th = ang(t);
        % analyze the radius, from 0 up to rs (step = ds)
        for r = 0:ds:rs % check ONLY ON UP TO THE SENSING RADIUS
            % compute point on the radius
            x = pr(1) + r*cos(th);
            y = pr(2) + r*sin(th);

            % check for obstacles / out of bounds (default=1 --> obstacle)
            if interp2(X, Y, obsNum, x, y, 'nearest', 1) ~= 0
                break;  % interrupt (found an obstacle OR out of bound)
            end
            
            % if the cycle didn't break, means (x,y) is NOT an obstacle point 
            % project to the nearest point of the grid (the grid has a
            % finite resolution so..)
            [~, jx] = min(abs(xv - x));
            [~, iy] = min(abs(yv - y));

            % visibility mask
            visible(iy, jx) = true;
        end
    end

    % intersection (just for safety, not necessary if the loop is on rs)
    Rmask = isfinite(rs) & ((X-pr(1)).^2 + (Y-pr(2)).^2 > rs^2);
    if any(Rmask(:))
        visible(Rmask) = false;
    end

    % W = S inters V ; W_rs was already cut (from the for loop on rs)
    Wrs = visible & Vmask;

    % fallback: in case of radius too small --> use classic voronoi
    if ~any(Wrs(:)), Wrs = Vmask; end
end
