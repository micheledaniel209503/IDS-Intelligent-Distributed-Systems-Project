function [u, omega] = ROB_control(state_rob, state_pkg, r_goal, Kp_u, Kp_theta)
%UNTITLED Summary of this function goes here

% TO BE TESTED
x_r = state_rob(1); y_r = state_rob(2); theta_r = state_rob(3);
x_p = state_pkg(1); y_p = state_pkg(2);

dx = x_p - x_r; dy = y_p - y_r;
d = hypot(dx, dy);

ang_req = atan2(dy, dx);
dtheta = atan2(sin(ang_req - theta_r), cos(ang_req - theta_r)); 

if (d < r_goal)
    u = 0; omega = 0;
    arrived = true;
    return;
end

u = Kp_u*d;
omega = Kp_theta*dtheta;
arrived = false;

end