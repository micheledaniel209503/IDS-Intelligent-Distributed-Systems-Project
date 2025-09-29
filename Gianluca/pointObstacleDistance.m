

function d = pointObstacleDistance(point, obstacle)
    % Distance between a point and an obstacle
    switch obstacle.type
        case 'c' % circle
            center = obstacle.state;
            r = obstacle.l/2;
            d = norm(point - center) - r;
        case 's' % square (axis-aligned)
            cx = obstacle.state(1); cy = obstacle.state(2);
            half = obstacle.l/2;
            dx = max(abs(point(1)-cx) - half, 0);
            dy = max(abs(point(2)-cy) - half, 0);
            d = sqrt(dx^2 + dy^2);
        otherwise
            d = inf;
    end
end