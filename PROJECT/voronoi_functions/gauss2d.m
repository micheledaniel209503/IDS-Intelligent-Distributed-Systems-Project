function rho = gauss2d(X, Y, mu, sigma_x2, sigma_y2)
% Works with no correlation rho
% X, Y : meshgrid matrices
% mu   : [mu_x, mu_y] location of the pdf

dx = X - mu(1);
dy = Y - mu(2);

% pdf expr
rho = (1/(2*pi*sqrt(sigma_x2 * sigma_y2))) .* ...
      exp( - ( dx.^2 ./ (2*sigma_x2) + dy.^2 ./ (2*sigma_y2) ) );

end
