function u = formation_control(Robots, ids_group, i_robot, form_R, target, w_form, w_att)

%   u = formation_control(Robots, ids_group, i_robot, form_R, target, w_form, w_att)
%
%   Inputs:
%       Robots      - struct array containing robot states (each with field .x_est)
%       ids_group   - vector of robot IDs involved in transporting the package
%       i_robot     - ID (index) of the specific robot to compute the control for
%       form_R      - desired radius of the circular formation
%       target      - 2D vector [x; y] of the desired package goal position
%       w_form      - weight for the formation control term
%       w_att       - weight for the attraction-to-target control term
%
%   Output:
%       u           - 2D normalized control vector to be applied to the robot
%
%   Notes:
%       - The function enforces a circular formation around the package center.
%       - Both formation and attraction terms are normalized and then combined.
%       - If a control term is zero, its direction is ignored in normalization.

    Nrobots = length(ids_group);

    % Get position of the robot i in 2D (estimated)
    pos_i = Robots(i_robot).state_est(1:2);

    %% Formation control term
    u_form = [0; 0];
    for k = 1:Nrobots
        j = ids_group(k);
        if j == i_robot
            continue;
        end

        pos_j = Robots(j).state_est(1:2);
        diff = pos_i - pos_j;
        dij = norm(diff);

        % Compute desired chord distance in circular formation
        % k_sep = circular separation index (shortest wrap-around)
        idx_i = find(ids_group == i_robot);
        idx_j = find(ids_group == j);
        k_sep = abs(idx_i - idx_j);
        k_sep = min(k_sep, Nrobots - k_sep);

        d_des_ij = 2 * form_R * sin(pi * k_sep / Nrobots);

        % Accumulate formation control contribution
        u_form = u_form - (dij^2 - d_des_ij^2) * diff;
    end

    % Normalize formation control (avoid NaN if zero)
    if norm(u_form) > 1e-8
        u_form = u_form / norm(u_form);
    else
        u_form = [0; 0];
    end

    %% Attraction control term
    % Vector from current robot position to target centroid
    u_att = target - pos_i;
    if norm(u_att) > 1e-8
        u_att = u_att / norm(u_att);
    else
        u_att = [0; 0];
    end

    %% Weighted combination
    u = w_form * u_form + w_att * u_att;

    % Normalize total control
    if norm(u) > 1e-8
        u = u / norm(u);
    end
end