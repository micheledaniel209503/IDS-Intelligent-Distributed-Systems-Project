function O = spawn_obstacles(N, map, l)
% DESCRIPTION: spawns N obstacles of characteristic dimension l inside the
% map
% N: number of obstacles to spawn
% mapW, mapH: map dimensions (0..W, 0..H)
% l: characteristic dimension
% Output: array of obstacles
    O = repmat(obs, N, 1);
    for k = 1:N
        if rand <= 0.55 % rand generates a random number between 0 and 1 so it's 50-50 between circle and square
            type = 'c';
            x = l/2 + (map.W - 2*l)*rand;
            y = l/2 + (map.H - 2*l)*rand;
        else
            type = 's';
            half = l/2;
            x = half + (map.W - 2*half)*rand;
            y = half + (map.H - 2*half)*rand;
        end
        O(k).state = [x, y];
        O(k).id = k;
        O(k).type = type;
        O(k).l = l;
    end
end
