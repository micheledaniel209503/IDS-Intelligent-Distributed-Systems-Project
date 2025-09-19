function Wrs = build_Wrs(pr, Vmask, obs, X, Y, rs, ds)

    if ~any(Vmask(:)) % pathological case (unlikely)
        Wrs = false(size(Vmask)); return
    end

    % Visibility set S (line-of-sight robot->point) only on V points
    idx = find(Vmask); % all points on V to be checked
    Smask = false(size(Vmask));  % initialization of the visibility mask (matrix same size as V but all FALSE elements)

    for k = 1:numel(idx) % go through all Vmask points and check visibility
        ii = idx(k);
        qx = X(ii); qy = Y(ii); % select point
        if is_line_free(pr(1), pr(2), qx, qy, obs, X, Y, ds) % if point is visible (1)
            Smask(ii) = true; % put it on the visibility set
        end
    end

    % INTERSECT W = S inters V 
    Wmask = Vmask & Smask;

    % SENSING RANGE W_rs
    Rmask = ( (X - pr(1)).^2 + (Y - pr(2)).^2 ) <= rs^2; % circle of radius rs
    % intersect W_rs = W inters R
    Wrs = Wmask & Rmask;


    if ~any(Wrs(:)) % very pathological case in which rs is really small, cutting all points of the cell
        Wrs = Wmask;
    end
end


% CHECKS FOR OBSTACLES IN THE LINE CONNECTING P1 and P2
function tf = is_line_free(x1,y1,x2,y2,obs,X,Y,ds)
    % sample the line connecting p1 (prob) and p2 (point) and check for obstacles
    L = hypot(x2-x1, y2-y1);

    if L < eps % means i'm checking the robot's position
        tf = true; 
        return; 
    end

    t = 0:ds:L;  % sample line
    xs = x1 + (x2-x1).*t/L; % line points
    ys = y1 + (y2-y1).*t/L; % line points

    obsNum = double(obs); % matrix of 0/1 --> 1=obstacle, 0=free
    v = interp2(X, Y, obsNum, xs, ys, 'nearest', 1); % VISIBILITY: interpolate the obstacle line on the grid 1=obstacle, 0=free; "1" is the out-of-grid points
    
    % tf will be: 0 if all points on line connecting p1 and p2 are free; 1 if there is any obstacle on line connecting p1 and p2
    tf = true;   % initialize
    for k = 1:numel(v) % check
        if v(k) ~= 0   % found an obstacle on the line
            tf = false;
            break; % exit the cicle (found an obstacle)
        end
    end 

end