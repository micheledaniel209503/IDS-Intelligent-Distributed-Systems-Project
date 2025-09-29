function rho = ring2d(X, Y, mu, R, sigma_r, eps_bg)
    % RING_PDF Function to compute the ring probability density function
    % Set default value for background noise (valley of the pdf) if not provided
    if nargin < 6 || isempty(eps_bg), eps_bg = 1e-15; end
    R = 0.97*R; % small compensation for the asimmetry of the pdf
    % Calculate the differences from the mean position
    DX = X - mu(1);
    DY = Y - mu(2);
    % Compute the radial distance from the mean position
    rr = hypot(DX, DY);
    % Calculate the ring probability density function
    rho = eps_bg + exp(- (rr - R).^2 / (2*sigma_r^2));
end