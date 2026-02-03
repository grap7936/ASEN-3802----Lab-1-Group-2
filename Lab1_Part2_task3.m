function [P,a] = Lab1_Part2_task3(F3D,P,I,A,c)

% outputs:
    % P, one of the forces
    % a, distance of force from edge

% inputs (left to right, only need F3D and F2):
    % F3D, inline stress gague
    % F2, data from standalone sensor
    % I, moment of inertia
    % A, % area of one rod
    % c, distance from center to edge of truss (h/2)

% default inputs:
arguments
    F3D
    P
    I = 5.9458; % in^4
    A = pi*(3/16)^2-pi*(2/16)^2; % in^2
    c = 250/25.4/2; % in
end

Fi = F3D - F3D(1);

a = Fi.*I./A./c./P;

end
