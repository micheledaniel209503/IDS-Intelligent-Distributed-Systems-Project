function p_des = TRANSPORT_control(Robots, ids_group, i_robot, form_R, target, w_form, w_att, w_damp, u_max, dt)
%FORMATION_CONTROL_DAMPED_RELATIVE Formation control with relative velocity damping
%
%   Adds damping that penalizes relative velocity differences between robots
%   to prevent oscillations and improve stability.

    Nrobots = length(ids_group);
    pos_i = Robots(i_robot).state_est(1:2);
    pos_i = pos_i(:); % column vector
    target = target(:); % column vector

    %% --- Formation control term ---
    u_form = [0; 0];
    u_damp = [0; 0];   % Initialize relative damping accumulator

    % Get velocity of robot i
    vel_i = compute_velocity_components(Robots(i_robot))'; 

    for k = 1:Nrobots
        j = ids_group(k);
        if j == i_robot
            continue;
        end

        % Positions and velocities of neighbor j
        pos_j = Robots(j).state_est(1:2);
        pos_j = pos_j(:); % column vector
        vel_j = compute_velocity_components(Robots(j))'; 

        % --- Formation control part ---
        diff = pos_i - pos_j;
        dij = norm(diff);

        idx_i = find(ids_group == i_robot);
        idx_j = find(ids_group == j);
        k_sep = abs(idx_i - idx_j);
        k_sep = min(k_sep, Nrobots - k_sep);

        d_des_ij = 2 * form_R * sin(pi * k_sep / Nrobots);
        u_form = u_form - (dij^2 - d_des_ij^2) * diff;

        % --- Relative damping part ---
        rel_vel = vel_i - vel_j;   % velocity difference
        rel_speed = norm(rel_vel);

        if rel_speed > 0.01
            damping_gain = min(1, rel_speed / 0.5);
            u_damp = u_damp - damping_gain * rel_vel;
        end
    end

    %% --- Attraction control term ---
    u_att = target - pos_i;
    u_att = u_att / (norm(u_att) + 0.01);

    %% --- Weighted combination ---
    u = w_form * u_form + w_att * u_att + w_damp * u_damp;
    p_des = pos_i + (dt/u_max) * u;
    p_des = p_des(:)';  % row vector output

end

function vel = compute_velocity_components(robot)
    theta = robot.state_est(3);
    u = robot.u;
    vel = [u * cos(theta), u * sin(theta)];
end