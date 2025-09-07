function rho = gauss2d(X, Y, mu, Sigma)
% DESCRIPTION: computes expression of bivariate gaussian pdf given meshgrid
% (map grid) and 2 parameters. The output is the bivariate gaussian pdf
% expressed point by point on the given [X, Y] grid
% X,Y: meshgrid
% mu = [mx,my] : target point
% Sigma 2x2 : covariance matrix (use eye)
Sinv = inv(Sigma);
detS = det(Sigma);
dx = X - mu(1); dy = Y - mu(2);
q = Sinv(1,1)*dx.^2 + 2*Sinv(1,2).*dx.*dy + Sinv(2,2)*dy.^2;
rho = exp(-0.5*q) / (2*pi*sqrt(detS));
end
