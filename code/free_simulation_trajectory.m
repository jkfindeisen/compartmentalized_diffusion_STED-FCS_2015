function [x, y] = free_simulation_trajectory(dr, R, Nt, x0, y0)
% Single molecule free diffusion trajectory.
%
%   dr  average displacement per time step [m]
%   R   sandbox radius [m]
%   Nt  number of time steps
%   x0, y0  starting point [m]
%
% All parameters must be real, scalar values. dr, R, and Nt must be
% positive. Nt must be an integer number sqrt(x0^2+y0^2) must be smaller
% than R.
%
% Author: Jan Keller-Findeisen (jan.keller@mpibpc.mpg.de)
% Dep. NanoBiophotonics, MPI for Biophysical Chemistry, Goettingen, Germany
% Date: 2013/06-2014/01

% parameter check
assert(nargin == 5, 'Wrong number of arguments!');

R2 = R^2;
if x0^2 + y0^2 >= R
    error('Initial point outside of simulation area');
end

% draw N random number normally distributed and scale with dr = steps
dx = dr * randn(Nt, 1);
dy = dr * randn(Nt, 1);
% sqrt(mean(dx.^2+dy.^2))/dr ~ sqrt(2) because two directions

% loop over all time steps
x = zeros(Nt, 1);
y = zeros(Nt, 1);
for ki = 1 : Nt
    
    % update position
    x0 = x0 + dx(ki);
    y0 = y0 + dy(ki);
    
    % check for being outside of sandbox
    if x0^2 + y0^2 >= R2
        % put back in sandbox on other side
        phi = atan2(y0, x0);
        x0 = x0 - 2 * R * cos(phi);
        y0 = y0 - 2 * R * sin(phi);
    end
    
    % store
    x(ki) = x0;
    y(ki) = y0;
end

end