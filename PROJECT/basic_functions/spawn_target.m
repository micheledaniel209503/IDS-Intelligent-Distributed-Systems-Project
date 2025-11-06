function [x,y] = spawn_target(outbound)
% SPANWS a package
% Given outbound polygon area, spawns a target point in a random position inside
% this area, giving as output x,y coordinates

xmin = min(outbound.polygon(:,1)); xmax = max(outbound.polygon(:,1));
ymin = min(outbound.polygon(:,2)) + 2; ymax = max(outbound.polygon(:,2) - 5);

% spawn package
x = xmin + rand*(xmax - xmin);
y = ymin + rand*(ymax - ymin);

end