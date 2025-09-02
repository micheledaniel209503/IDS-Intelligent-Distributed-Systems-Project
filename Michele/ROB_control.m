function [u, omega] = ROB_control(state_rob, point_trg, u_sat, omega_sat, Kp_u, Kp_theta, r_goal, theta_goal)
% DESCRIPTION: point-to-point control for a unicycle robot. Control is
% provided by the two velocities that are given as outputs of the function.
% CONTROL LOGIC: given the position of the target (point_trg [x,y]), first orient the
% robot to point to the target (within an angular tolerance), then proceed
% walking towards it, controlling the orientation continuously. Once
% reached (stops at r_eps distance from target), re-orient the robot
% pointing in the direction given by theta_goal (rad)

% util
wrap = @(a) atan2(sin(a), cos(a)); % wrap to choose the shortest rotation --> may need to unwrap the resulting rob.state(3)
sat  = @(v,lim) max(min(v, lim), -lim);

% tolerances/thresholds
r_eps = 0.15;
theta_eps_s = deg2rad(8); % [rad] tolerance on initial orientation of rob
theta_eps_e = deg2rad(1); % [rad] tolerance on final orientation of rob

% positions
x_r = state_rob(1); y_r = state_rob(2); theta_r = state_rob(3);
x_t = point_trg(1); y_t = point_trg(2);

% errors
dx = x_t - x_r; dy = y_t - y_r;
d = hypot(dx, dy);
ang_req = atan2(dy, dx);
e_trg = wrap(ang_req - theta_r); % ang error to trg 
e_goal  = wrap(theta_goal - theta_r); % ang error to final orientation

% stopping condition (if rob is close enough to trg)
if (d - r_goal < r_eps)
    u = 0;
    if abs(e_goal) < theta_eps_e % if already aligned with goal
        omega = 0; return;
    end
    % otherwise re-orient robot to theta_goal
    omega = sat(Kp_theta * e_goal, omega_sat); % angle control, objective: final orientation
    return;
end

% angle control, objective: trg (always active and corrective)
omega = sat(Kp_theta * e_trg, omega_sat);

% first: robot points to trg
if abs(e_trg) > theta_eps_s % not aligned
    u = 0; % don't start, keep aligning with trg
else % aligned --> proceed
    e_d = max(d - r_goal, 0); % deadzone on distance
    u = min(max(Kp_u*e_d, 0), u_sat); % distance control, objective: trg
end

end