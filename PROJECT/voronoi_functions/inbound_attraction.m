function rho = inbound_attraction(X, Y, fraction, p)
% PDF that makes robot attracted to the lower portion of the map, where the
% inbound region is located.
% fraction: portion of the map that is considered "attractive"
    Ymin = min(Y(:));
    Ymax = max(Y(:));
    Yth  = Ymin + fraction*(Ymax - Ymin);
    den  = max(eps, (Ymax - Yth));

    rho = zeros(size(Y));

    plateau = 10;
    base = 0;

    above = (Y > Yth);
    rho(above) = base + (plateau - base) .* (Ymax - Y(above)) ./ den;

    rho(~above) = plateau;
end
