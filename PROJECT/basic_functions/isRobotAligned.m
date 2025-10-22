function isAligned = isRobotAligned(robot, objectPosition, desiredRadius, tolerance)
% ISROBOTALIGNED checks if a single robot is correctly positioned
% on the circular formation around the object within a given tolerance.
%
% INPUTS:
%   robot          : instance of class 'rob'
%   objectPosition : 1x2 vector [x, y] of the object's position
%   desiredRadius  : scalar, desired radius of circular formation
%   tolerance      : scalar, acceptable deviation from desired radius
%
% OUTPUT:
%   isAligned      : boolean, true if the robot is within tolerance

    % Extract robot's estimated position
    robotPos = robot.state_est(1:2);  
    
    % If position is invalid (NaN), immediately return false
    if any(isnan(robotPos))
        isAligned = false;
        return;
    end
    
    % Compute distance from robot to object
    distance = norm(robotPos - objectPosition);
    
    % Compute deviation from desired radius
    deviation = abs(distance - desiredRadius);
    
    % Check alignment condition
    isAligned = deviation <= tolerance;
end