function [H,z,C] = trilateration(anchors, distances, noise_std)
    % Number of anchors
    n = size(anchors, 1);
    
    % Initialize matrices
    H = zeros(n-1, 2);
    z = zeros(n-1, 1);
    C = zeros(n-1);
    
    % Iterate over all anchors
    for i = 1:n-1
        % Fill the matrices
        H(i, :) = 2*[anchors(i+1, 1) - anchors(i, 1), anchors(i+1, 2) - anchors(i, 2)];
        z(i) = - distances(i+1)^2  + distances(i)^2 + anchors(i+1, 1)^2 - anchors(i, 1)^2 + anchors(i+1, 2)^2 - anchors(i, 2)^2;
        % Fill the covariance matrix
        if i == 1
            C(i,i) = 4 * noise_std^2 * (distances(i+1)^2 + distances(i)^2);
            if n > 2
                C(i,i+1) = -4 * noise_std^2 * distances(i+1)^2;
            end
        elseif i < n-1
            C(i,i-1) = -4 * noise_std^2 * distances(i)^2;
            C(i,i) = 4 * noise_std^2 * (distances(i+1)^2 + distances(i)^2);
            C(i,i+1) = -4 * noise_std^2 * distances(i+1)^2;
        else
            C(i,i-1) = -4 * noise_std^2 * distances(i)^2;
            C(i,i) = 4 * noise_std^2 * (distances(i+1)^2 + distances(i)^2);
        end
    end
end