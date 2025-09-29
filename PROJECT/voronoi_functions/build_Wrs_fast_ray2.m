function Wrs = build_Wrs_fast_ray2(pr, Vmask, obs, X, Y, sr, ds)
% BUILD_WRS_FAST_RAY (grid-based) — W_rs = (S ∩ V) ∩ B_rs(pr) via ray
% casting USING INDEX COORDINATES
% pr     = [x_r,y_r]
% Vmask  = Voronoi cell of robot
% obs    = ~free_mask (true = obstacle)
% X,Y    = meshgrid
% sr     = sensing radius (Inf or [] = no limit --> classic voronoi)
% ds     = step size along radius

    if ~any(Vmask(:)), Wrs = false(size(Vmask)); return; end % pretty unlikely

    % how many rotations
    ntheta = 256;
    % grid parameters
    xv = X(1,:);
    yv = Y(:,1);
    dx = xv(2) - xv(1); % step size along x
    dy = yv(2) - yv(1); % step size along y
    [m,n] = size(X);

    % conversion world coordinates --> indices of the grid
    % jx = round((x - xv(1))/dx) + 1
    % iy = round((y - yv(1))/dy) + 1
    invdx = 1/dx; 
    invdy = 1/dy; 
    x0 = xv(1); 
    y0 = yv(1);

    % world coordinates --> indeces for ROBOT POSITION
    ixr = (pr(1) - x0)*invdx + 1;
    iyr = (pr(2) - y0)*invdy + 1;

    % initialization
    visible = false(m,n);

    % angles to be scanned (-180 --> 180 deg)
    ang = linspace(-pi, pi, ntheta+1); ang(end) = [];

    % number of steps on the scanned radius (scan all the radius up to sr,
    % with step size = ds)
    nsteps = floor(sr/ds);

    for t = 1:numel(ang) % loop through each angle
        th = ang(t);

        % step increment (along radius) in index space (increment index given coordinate
        % increment)
        % NOTE: if theta = 0 we would have dix = ds/dx and diy = 0 but we are
        % scanning 360° so
        dix = (cos(th) * ds) * invdx;   % column step
        diy = (sin(th) * ds) * invdy;   % row step

        % start from robot position (index coords)
        ix = ixr;
        iy = iyr;

        for s = 0:nsteps
            % closest cell point
            jx = round(ix);  iy_int = round(iy);

            % map boundaries check (point we are checking is outside map
            % boundaries)
            % note: jx<1 or iy<1 means negative points (outside map)
            if jx < 1 || jx > n || iy_int < 1 || iy_int > m
                break; % stop
            end

            % found an obstacle on the point (check point on the obstacle
            % mask)
            if obs(iy_int, jx)
                break;
            end

            % if the cycle didn't break, means point is NOT an obstacle point 
            % so mark it as visible on the visibility matrix (still index
            % coordinates)
            % REMINDER: meshgrid in matlab X(i,j) = xv(j), Y(i,j) = yv(i)
            visible(iy_int, jx) = true;

            % advance on the radius
            ix = ix + dix;
            iy = iy + diy;

        end
    end
    
    % disk intersection (just for safety, not necessary if the loop is on sr)
    Rmask = (X - pr(1)).^2 + (Y - pr(2)).^2 <= sr^2;

    % final intersection: W_rs = S inters V inters disk
    Wrs = visible & Vmask & Rmask;

    if ~any(Wrs(:)), Wrs = Vmask; end % in the very unlucky case of no visibility for example (robot on an obstacle maybe) --> give full voronoi
end
