function C_proj = project_centroid(C, Wmask, X, Y)
% PROJECT_CENTROID_TO_WRS
% Projects the point C onto the nearest point belonging to the mask Wmask
% (typically Wrs_set{k}), using the grid X,Y
%
% C      : [Cx, Cy] (coordinates in world)
% Wmask  : logical mask of the visible Voronoi Wrs for that robot
% X,Y    : meshgrid of the map

    % if centroid is undefined -> don't do anything
    if any(~isfinite(C))
        C_proj = C;
        return;
    end
    % grid
    [m,n] = size(X);
    dx = X(1,2) - X(1,1);
    dy = Y(2,1) - Y(1,1);
    x0 = X(1,1);
    y0 = Y(1,1);

    % grid coordinates of the centroid
    ix = round((C(1) - x0) / dx) + 1;
    iy = round((C(2) - y0) / dy) + 1;

    % if centroid is already on Wrs -> keep it and ignore any additional
    % computation
    if ix >= 1 && ix <= n && iy >= 1 && iy <= m && Wmask(iy, ix)
        C_proj = C;
        return;
    end

    % if centroid is NOT on Wrs -> look for the point nearest (euclidean distance) to C that is
    % in Wrs
    D2 = (X - C(1)).^2 + (Y - C(2)).^2; % euclidean distance between every point of the grid and centroid
    D2(~Wmask) = Inf; % ignore points that are not on the Wrs mask

    [min_val, idx_min] = min(D2(:)); % find the point of minimum distance, store the value of distance and its position

    if isinf(min_val) % very unfortunate case: if the distance is not finite, keep C
        C_proj = C;
    else
        [iy_min, ix_min] = ind2sub(size(D2), idx_min); % store the position
        C_proj = [X(iy_min, ix_min), Y(iy_min, ix_min)]; % return the projected centroid
    end
end
