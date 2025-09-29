function [x,y] = spawn_pkg(inbound)
% SPANWS a package
% Given inbound polygon area, spawns a package in a random position inside
% this area, giving as output x,y coordinates

xmin = min(inbound.polygon(:,1)); xmax = max(inbound.polygon(:,1));
ymin = min(inbound.polygon(:,2)) + 4; ymax = max(inbound.polygon(:,2));

% spawn package
x = xmin + rand*(xmax - xmin);
y = ymin + rand*(ymax - ymin);

end