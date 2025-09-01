function [u, omega] = ROB_control(state_rob, state_pkg, u_sat, omega_sat, Kp_u, Kp_theta)
r_goal = 0.5; % [m]
theta_th = deg2rad(8); % [rad] threshold on initial orientation of rob

% TO BE TESTED
x_r = state_rob(1); y_r = state_rob(2); theta_r = state_rob(3);
x_p = state_pkg(1); y_p = state_pkg(2);

dx = x_p - x_r; dy = y_p - y_r;
d = hypot(dx, dy);

ang_req = atan2(dy, dx);
dtheta = atan2(sin(ang_req - theta_r), cos(ang_req - theta_r));

% stopping condition (if rob is close enough to pkg)
if (d < r_goal)
    u = 0; omega = 0;
end

% angle control (always active and corrective)
omega = max(min(Kp_theta * dtheta, omega_sat), -omega_sat);

% first: robot points to pkg
if abs(dtheta) > theta_th % not aligned
    u = 0; % don't start
    omega = min(Kp_theta*dtheta, omega_sat);
else % aligned --> proceed
    e_d = max(d - r_goal, 0); % err
    u = min(max(Kp_u*e_d, 0), u_sat);
end

end