function trajectory = planTrajectory(mapSize, obstacles, startPoint, goalPoint, nPoints, minDistance)
    % Use fewer optimization variables (control points)
    nCtrl = 15; % number of control points (including start & goal)
    x0 = linspace(startPoint(1), goalPoint(1), nCtrl);
    y0 = linspace(startPoint(2), goalPoint(2), nCtrl);
    initGuess = [x0(2:end-1), y0(2:end-1)]'; 
    initGuess = initGuess(:);

    % Optimization setup (penalty-based, no constraints)
    options = optimoptions('fmincon','Display','none','Algorithm','sqp','MaxIterations', 50);
    problem.objective = @(z) costFunctionSoft(z, startPoint, goalPoint, nCtrl, obstacles, minDistance);
    problem.x0 = initGuess;
    problem.lb = repmat([0;0], nCtrl-2, 1);
    problem.ub = repmat(mapSize(:), nCtrl-2, 1);
    problem.solver = 'fmincon';
    problem.options = options;

    zOpt = fmincon(problem);

    % Rebuild trajectory with spline interpolation
    ctrlPoints = [startPoint; reshape(zOpt,[],2); goalPoint];
    t = linspace(0,1,nCtrl);
    tt = linspace(0,1,nPoints);
    xs = spline(t, ctrlPoints(:,1), tt);
    ys = spline(t, ctrlPoints(:,2), tt);

    traj_non_resampled = [xs(:), ys(:)];
    trajectory = resampleTrajectory(traj_non_resampled, nPoints);
end


function J = costFunctionSoft(z, startPoint, goalPoint, nCtrl, obstacles, minDistance)
    % Set weights for each task:
    w_pathLength = 100;
    w_distance = 10;

    % Build control points
    ctrlPoints = [startPoint; reshape(z,[],2); goalPoint];

    % Interpolated trajectory for evaluation
    t = linspace(0,1,nCtrl);
    tt = linspace(0,1,50); % fewer samples for cost
    xs = spline(t, ctrlPoints(:,1), tt);
    ys = spline(t, ctrlPoints(:,2), tt);
    waypoints = [xs(:), ys(:)];

    % Path length
    diffs = diff(waypoints,1,1);
    pathLength = sum(vecnorm(diffs,2,2));

    % Obstacle penalty
    penalty = 0;
    for i = 1:size(waypoints,1)
        for j = 1:numel(obstacles)
            d = pointObstacleDistance(waypoints(i,:), obstacles(j));
            penalty = penalty + max(0, minDistance - d)^2;
        end
    end

    % Weighted cost
    J = w_pathLength * pathLength + w_distance * penalty;
end

function trajResampled = resampleTrajectory(trajectory, nPoints)
    % trajectory: [N x 2] original path
    % nPoints: desired number of points in resampled trajectory

    % Compute cumulative arc length
    diffs = diff(trajectory,1,1);
    segLen = vecnorm(diffs,2,2);
    s = [0; cumsum(segLen)];

    % Normalize arc length
    sNorm = s / s(end);

    % Interpolate to uniform arc length spacing
    sNew = linspace(0,1,nPoints);
    xNew = interp1(sNorm, trajectory(:,1), sNew);
    yNew = interp1(sNorm, trajectory(:,2), sNew);

    trajResampled = [xNew(:), yNew(:)];
end