function plot_phi_debug(Phi, X, Y, L, Robots, idx_use, centroids, opts)
% Phi: m×n×N
% X,Y: meshgrid
% L  : labels (1..N)
% Robots: array
% idx_use: indeces of used robot on voronoi
% centroids: N×2
% opts.logscale (bool), opts.mode ('auto'|'lineup'|'ring'), opts.mu_pkg (1×2), opts.R (scalar)

    if nargin < 8, opts = struct; end
    if ~isfield(opts,'logscale'), opts.logscale = false; end
    if ~isfield(opts,'mode'),     opts.mode     = 'auto'; end
    if ~isfield(opts,'mu_pkg'),   opts.mu_pkg   = [NaN NaN]; end
    if ~isfield(opts,'R'),        opts.R        = NaN; end

    xv = unique(X(1,:));  yv = unique(Y(:,1));
    dx = mean(diff(xv));  dy = mean(diff(yv));
    cell_area = dx*dy;

    N = size(Phi,3);
    tl = tiledlayout(ceil(sqrt(N)), ceil(sqrt(N)), "TileSpacing","compact","Padding","compact");

    for k = 1:N
        nexttile; hold on; axis equal tight; set(gca,'YDir','normal'); grid on
        W = Phi(:,:,k);
        if opts.logscale
            imagesc(xv, yv, log10(W + eps)); colorbar; titleStr = sprintf('log10 \\Phi_k, k=%d',k);
        else
            imagesc(xv, yv, W); colorbar; titleStr = sprintf('\\Phi_k, k=%d',k);
        end

        try
            contour(X, Y, double(L==k), [0.5 0.5], 'k', 'LineWidth', 1.0);
        catch
        end

        % rob + cent
        id = idx_use(k);
        pr = Robots(id).state(1:2);
        plot(pr(1), pr(2), 'wo', 'MarkerFaceColor','k', 'MarkerSize',4);
        if all(isfinite(centroids(k,:)))
            pc = centroids(k,:);
            plot(pc(1), pc(2), 'r+', 'MarkerSize',8,'LineWidth',1.2);
        end

        modek = opts.mode;
        if modek=="auto"
            ws = lower(string(Robots(id).working_state));
            if any(ws==["l","lineup"]), modek = "lineup";
            elseif any(ws==["r","ring"]), modek = "ring";
            else, modek = "free";
            end
        end

        switch modek
            case "lineup"
                mu = Robots(id).target;
                if all(isfinite(mu))
                    plot(mu(1), mu(2), 'bx', 'MarkerSize',8, 'LineWidth',1.2);
                    titleStr = titleStr + sprintf(' | lineup: \\mu=[%.1f %.1f]', mu(1),mu(2));
                end
            case "ring"
                mu = opts.mu_pkg;
                if all(isfinite(mu))
                end
        end

        % mass
        mask_k = (L==k);
        Mk = sum(W(mask_k), 'all') * cell_area;
        title(sprintf('%s | M=%.2e (id=%d)', titleStr, Mk, id));
    end
    title(tl, 'Debug \Phi per robot (bordo Voronoi, robot, centroide, target/anello)');
end
