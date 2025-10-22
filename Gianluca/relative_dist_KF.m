function rob = relative_dist_KF(ri, rj, j, z_meas, R_meas, dt)
% UPDATE_RELATIVE_DISTANCE
% EKF scalar filter for the relative distance between robot i and robot j
%
% INPUTS:
%   robots - array of rob objects
%   i, j   - indices of the current robot and the other robot
%   z_meas - measured distance from i to j (scalar). If empty or NaN, only predict
%   dt     - timestep (scalar)
%
% OUTPUT:
%   ri  - i-th robot with updated estimate of the relative distance

% Previous estimate and covariance
r_prev = ri.rel_dist(j);
P_prev = ri.P_rel_dist(j);

% Compute relative vector and unit vector
d = rj.x_est(1:2) - ri.x_est(1:2);
r_hat_pos = norm(d);
if r_hat_pos < 1e-12
    e = [1;0]; % fallback if positions are identical
else
    e = d / r_hat_pos;
end

% Compute relative velocity vector
vi_vec = ri.u * [cos(ri.x_est(3)); sin(ri.x_est(3))];
vj_vec = rj.u * [cos(rj.x_est(3)); sin(rj.x_est(3))];
v_rel = vj_vec - vi_vec;

% --- Prediction step ---
alpha = e' * v_rel;
r_pred = r_prev + dt * alpha;
P_pred = P_prev; 

% --- Update step (if measurement available) ---
if isempty(z_meas) || isnan(z_meas)
    r_plus = r_pred;
    P_plus = P_pred;
else
    H = 1;  % measurement model: h(r) = r
    S = H * P_pred * H' + R_meas;
    K = (P_pred * H') / S;
    resid = z_meas - r_pred;
    
    r_plus = r_pred + K * resid;
    P_plus = (1 - K * H) * P_pred;
end

% Store updated values back in the robot object
ri.rel_dist(j) = r_plus;
ri.P_rel_dist(j) = P_plus;

rob = ri;

end
