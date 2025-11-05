function [p_des, center] = TRANSPORT_control(Robots, ids_group, i_robot, form_R, target, ...
                                   w_form, w_att, w_obs, u_max, dt, Obstacles)
% INPUTS:
%   Robots     - array of robot objects
%   ids_group  - IDs of robots in formation
%   i_robot    - index of current robot
%   form_R     - formation radius
%   target     - target position [x; y]
%   w_form, w_att, w_obs - control weights
%   u_max  - maximum linear velocity
%   dt         - time step
%   Obstacles  - array of obstacle objects



    %% Parameters
    
    prob_comm = 0.7;   % Probability of successful communication
    
    rep_range = 2;   % Start of repulsion
    
    rep_sigma = 0.5;   % Width of repulsion zone

    %% Initialization
    Nrobots = length(ids_group);
    pos_i = Robots(i_robot).state_est(1:2);
    pos_i = pos_i(:);
    target = target(:);

    %% Ordering of robots (angle around centroid)
    positions = zeros(2, Nrobots);
    for k = 1:Nrobots
        j = ids_group(k);
        positions(:, k) = Robots(j).state_est(1:2);
    end
    center = mean(positions, 2);
    angles = atan2(positions(2,:) - center(2), positions(1,:) - center(1));
    [~, sort_idx] = sort(angles);
    ids_group = ids_group(sort_idx);

    %% Initialize control terms
    u_form = [0; 0];
    u_obs  = [0; 0];

    %% Formation + Cooperative Obstacle Avoidance
    for k = 1:Nrobots
        j = ids_group(k);
        if j == i_robot
            continue;
        end
        
        if rand(1) <= prob_comm % Communication probability
            %% Formation control
            pos_j = Robots(j).state_est(1:2);
            pos_j = pos_j(:);

            diff = pos_i - pos_j;
            dij = norm(diff);

            idx_i = find(ids_group == i_robot);
            idx_j = find(ids_group == j);
            k_sep = abs(idx_i - idx_j);
            k_sep = min(k_sep, Nrobots - k_sep);

            d_des_ij = 2 * form_R * sin(pi * k_sep / Nrobots);
            u_form = u_form - (dij^2 - d_des_ij^2) * diff;

            %% Obstacle avoidance
            sr_j = Robots(j).sr;
            u_obs_j = [0; 0];

            for o = 1:length(Obstacles)
                pos_o = Obstacles(o).state(:);

                % Equivalent radius depending on obstacle shape
                if Obstacles(o).type == 'c'
                    r_eff = Obstacles(o).l / 2;   % circle
                else
                    r_eff = (Obstacles(o).l / 2) * sqrt(2);  % square ≈ circle of radius l/2 * sqrt(2)
                end

                diff_o = pos_j - pos_o;
                d = norm(diff_o) - r_eff;

                % Only consider obstacle if within sensing range
                if d < sr_j && d > 1e-3
                    % Repulsion intensity
                    phi = 10*exp(-((d - rep_range/2)^2) / (2*rep_sigma^2));
                    dir_obs = diff_o / (norm(diff_o) + 1e-6);

                    % Tangential component (in order to avoid motion stops)
                    dir_target = (target - pos_j);
                    dir_target = dir_target / (norm(dir_target) + 1e-6);
                    sign_tan = sign(det([dir_obs, dir_target]));   % direction of rotation
                    tang = sign_tan * [-dir_obs(2); dir_obs(1)];   % 90° rotated vector
                    phi_tan = 0.6 * phi;                           % tangential weight

                    % repulsive + tangential effects
                    u_obs_j = u_obs_j + phi * dir_obs + phi_tan * tang;
                end
            end

            % Accumulate
            u_obs = u_obs + u_obs_j;
        end
    end

    % Normalization
    if norm(u_obs) > 1e-3
        u_obs = u_obs / norm(u_obs);
    end

    %% Attraction term
    u_att = target - pos_i;
    u_att = u_att / (norm(u_att) + 0.01);

    %% Total input (sum of the different contributions)
    u = w_form * u_form + w_att * u_att + w_obs * u_obs;

    %% Reference position for next step 
    p_des = pos_i + (dt / u_max) * u; % do not exceed the maximum distance that the robot can cover in a single time step
    p_des = p_des(:)'; 

end
