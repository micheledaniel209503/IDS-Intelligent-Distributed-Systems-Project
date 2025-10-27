function p_des = TRANSPORT_control(Robots, ids_group, i_robot, form_R, target, w_form, w_att, w_damp, u_max, dt)
    % INPUTS:
    %   Robots     - structure array with fields: state_est = [x, y, theta], u = speed
    %   ids_group  - list of robot IDs that belong to this formation
    %   i_robot    - ID of the current robot being controlled
    %   form_R     - desired formation radius
    %   target     - [x; y] target position
    %   w_form     - weight for formation term
    %   w_att      - weight for attraction term
    %   w_damp     - weight for damping term
    %   u_max      - maximum control input magnitude
    %   dt         - time step
    %
    % OUTPUT:
    %   p_des      - desired next position [1x2] for robot i

    %% Initialization
    Nrobots = length(ids_group);
    pos_i = Robots(i_robot).state_est(1:2);
    pos_i = pos_i(:);
    target = target(:);

    %% Dynamic angular ordering of robots

    positions = zeros(2, Nrobots);
    for k = 1:Nrobots
        j = ids_group(k);
        positions(:, k) = Robots(j).state_est(1:2);
    end

    % centroid of the formation
    center = mean(positions, 2);

    % Compute angles of each robot around the center
    angles = atan2(positions(2,:) - center(2), positions(1,:) - center(1));

    % Sort robots counterclockwise by angle
    [~, sort_idx] = sort(angles);
    ids_group = ids_group(sort_idx);

    %% --- Formation control initialization ---
    u_form = [0; 0];
    u_damp = [0; 0];

    % Velocity of current robot i
    vel_i = compute_velocity_components(Robots(i_robot))';

    %% --- Iterate over all other robots in the group ---
    for k = 1:Nrobots
        j = ids_group(k);
        if j == i_robot
            continue;
        end

        % Neighbor position and velocity
        pos_j = Robots(j).state_est(1:2);
        pos_j = pos_j(:);
        vel_j = compute_velocity_components(Robots(j))';

        % Formation control
        diff = pos_i - pos_j;
        dij = norm(diff);

        idx_i = find(ids_group == i_robot);
        idx_j = find(ids_group == j);
        k_sep = abs(idx_i - idx_j);
        k_sep = min(k_sep, Nrobots - k_sep);

        % Desired theorical distance for ring formation
        d_des_ij = 2 * form_R * sin(pi * k_sep / Nrobots);

        u_form = u_form - (dij^2 - d_des_ij^2) * diff;

        % Relative velocity damping
        rel_vel = vel_i - vel_j;
        u_damp = u_damp - rel_vel;
    end
    %% Attraction to target
    u_att = target - pos_i;
    u_att = u_att / (norm(u_att) + 0.01);

    %% Weighted control combination
    u = w_form * u_form + w_att * u_att + w_damp * u_damp;

    %% Position computation
    p_des = pos_i + (dt / u_max) * u;
    p_des = p_des(:)'; 

end

% COMPUTE_VELOCITY_COMPONENTS: Convert robot heading and speed to XY velocity
function vel = compute_velocity_components(robot)
    theta = robot.state_est(3);
    u = robot.u;
    vel = [u * cos(theta), u * sin(theta)];
end