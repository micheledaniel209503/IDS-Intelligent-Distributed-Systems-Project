function O = spawn_obstacles(N, map, l, minDist)
% DESCRIPTION: spawns N obstacles of characteristic dimension l inside the
% map, ensuring a minimum distance between them
%
% Inputs:
%   N       : number of obstacles to spawn
%   map     : struct with fields .W and .H (map width and height)
%   l       : characteristic dimension of the obstacles
%   minDist : minimum allowed distance between obstacle centers
%
% Output:
%   O : array of obs objects

    O = repmat(obs, N, 1); % preallocate
    maxTries = 1000;       % max attempts to place each obstacle
    
    for k = 1:N
        placed = false;
        tries = 0;
        
        while ~placed && tries < maxTries
            tries = tries + 1;
            
            % Randomly choose type and position
            a = 0.3;
            if rand <= 0%0.55
                type = 'c';
                x = l/2 + (map.W - l) * rand;
                y = map.H*a + (map.H*(1-2*a) - l) * rand;
            else
                type = 's';
                half = l/2;
                x = half + (map.W - l) * rand;
                y = map.H*a + (map.H*(1-2*a) - l) * rand;
            end
            
            candidate = [x, y];
            
            % Check minimum distance with previously placed obstacles
            tooClose = false;
            for j = 1:k-1
                prev = O(j).state;
                d = norm(candidate - prev);
                if d < (minDist + max(l, O(j).l)/2) % conservative check
                    tooClose = true;
                    break;
                end
            end
            
            % If position is valid, save obstacle
            if ~tooClose
                O(k).state = candidate;
                O(k).id = k;
                O(k).type = type;
                O(k).l = l;
                placed = true;
            end
        end
        
        % If not placed after many tries
        if ~placed
            warning('Obstacle %d could not be placed with given constraints.', k);
            O(k) = []; % remove it
        end
    end
end
