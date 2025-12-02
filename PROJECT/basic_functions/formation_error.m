function err_avg = formation_error(Robots, ids_group, form_R)
% INPUTS:
%   Robots     - array of robot objects
%   ids_group  - array of robot IDs in the group
%   form_R     - desired formation radius
%
% OUTPUT:
%   err_avg    - average formation error (meters)

    %% --- Initialization ---
    Nrobots = length(ids_group);

    % Extract positions
    positions = zeros(2, Nrobots);
    for k = 1:Nrobots
        j = ids_group(k);
        positions(:, k) = Robots(j).state(1:2);
    end

    center = mean(positions, 2);

    %% --- Order robots by angle around center ---
    angles = atan2(positions(2,:) - center(2), positions(1,:) - center(1));
    [~, sort_idx] = sort(angles);
    positions  = positions(:, sort_idx);
    %% Initialization
    total_err = 0;
    count = 0;

    for a = 1:Nrobots
        pos_a = positions(:, a);

        for b = a+1:Nrobots
            pos_b = positions(:, b);

            % Actual distance
            dij = norm(pos_a - pos_b);

            % Index separation for desired inter-robot distance
            k_sep = abs(a - b);
            k_sep = min(k_sep, Nrobots - k_sep);

            % Theoretical distance for circular formation
            d_des = 2 * form_R * sin(pi * k_sep / Nrobots);

            % Accumulate absolute error
            total_err = total_err + abs(dij - d_des);
            count = count + 1;
        end
    end

    %% --- Average error---
    err_avg = total_err / count;

end
