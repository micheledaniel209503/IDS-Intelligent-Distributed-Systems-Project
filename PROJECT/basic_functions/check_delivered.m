function delivered = check_delivered(center_xy, target_xy, outboundPoly, tol_target)
% LOGIC: delivered = true if pkg is close enough to target (given a
% tolerance) AND it's inside the outbound zone
% center_xy = [cx cy] estimated position of pkg
% target_xy = [tx ty] target of the pkg
% outboundPoly = [N x 2] polygon
% tol_target = tolerance on target distance (es. max(pkgR, 1.0))

    cx = center_xy(1); cy = center_xy(2);
    close_to_trg = false;
    delivered = false; % initialize

    % criteria 1: close to target
    if norm(center_xy - target_xy) <= tol_target
        close_to_trg = true;
    end

    % criteria 2: inside outbound zone
    in = inpolygon(cx, cy, outboundPoly(:,1), outboundPoly(:,2));
    if in && close_to_trg
        delivered = true;
    end
end
