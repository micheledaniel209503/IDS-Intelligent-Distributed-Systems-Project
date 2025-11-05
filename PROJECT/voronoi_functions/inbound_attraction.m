function rho = inbound_attraction(X, Y, fraction, p)
% PDF that makes robot attracted to the lower portion of the map, where the
% inbound region is located.
% fraction: portion of the map that is considered "attractive"
    Ymin = min(Y(:));
    Ymax = max(Y(:));
    Yth  = Ymin + fraction*(Ymax - Ymin); % threshold y coord

    rho = 4.*ones(size(Y));

    above = (Y > Yth);
    % normalized vertical distance (coordinate) in the upper part
    t = (Y(above) - Yth) ./ max(eps, (Ymax - Yth));
    % negative slope with an exponent (1 - t)^p up to 0
    rho(above) = max(0, (1 - t).^p);
end
