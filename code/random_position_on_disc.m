function [x, y] = random_position_on_disc(R, N)
% Places N points randomly on a disc of radius R

assert(nargin == 2, 'Wrong number of arguments!');

phi = 2 * pi * rand(N, 1); % equally distributed on the angle
r = R * sqrt(rand(N, 1)); % more for higher radius, because circumference ~ R
x = r .* cos(phi);
y = r .* sin(phi);

end