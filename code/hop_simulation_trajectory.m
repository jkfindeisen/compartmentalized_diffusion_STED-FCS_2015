function [xt, yt, x1, y1, hops] = hop_simulation_trajectory(dr, R, Nt, x1, y1, hop_map, hop_pixelsize, hop_probability)
% Single molecule hopping diffusion trajectory given a hopping map.
%
%   dr  average displacement per time step [m]
%   R   sandbox radius [m]
%   Nt  number of time steps
%   x0, y0  starting point [m]
%   Hc hoping lattice constant
%   Hp hoping probability
%
% All parameters must be real, scalar values. dr, R, and Nt must be
% positive. Nt must be an integer number. sqrt(x0^2+y0^2) must be smaller
% than R.
%
% Author: Jan Keller-Findeisen (jan.keller@mpibpc.mpg.de)
% Dep. NanoBiophotonics, MPI for Biophysical Chemistry, Goettingen, Germany
% Date: 2013/06-2013/08

% parameter check
assert(nargin == 8, 'Wrong number of arguments!');

R2 = R^2;
if x1^2 + y1^2 >= R
    error('Initial point outside of simulation area');
end

% initial hop compartment

id1 = get_compartment(x1, y1, hop_map, hop_pixelsize);

% draw N random number normally distributed and scale with dr = steps
dx = dr * randn(Nt, 1);
dy = dr * randn(Nt, 1);
% sqrt(mean(dx.^2+dy.^2))/dr ~ sqrt(2) because two directions

% loop over all time steps
xt = zeros(Nt, 1);
yt = zeros(Nt, 1);
hops = false(Nt, 1);
for ki = 1 : Nt
    
    % update position
    x2 = x1 + dx(ki);
    y2 = y1 + dy(ki);
    
    % check for hoping
    id2 = get_compartment(x2, y2, hop_map, hop_pixelsize);
    
    % if we could hop but don't do it
    if id1 ~= id2 && rand() > hop_probability
        % find one inside the old compartment
        while id2 ~= id1
            x2 = x1 + dr * randn();
            y2 = y1 + dr * randn();
            id2 = get_compartment(x2, y2, hop_map, hop_pixelsize);
        end
    end
    
    % check for being outside of sandbox
    if x2^2 + y2^2 >= R2
        % put back in sandbox on other side
        phi = atan2(y2, x2);
        x1 = x2 - 2 * R * cos(phi);
        y1 = y2 - 2 * R * sin(phi);
        id1 = get_compartment(x1, y1, hop_map, hop_pixelsize);
    else
        x1 = x2;
        y1 = y2;
        id1 = id2;
    end
    
    % store
    xt(ki) = x1;
    yt(ki) = y1;
end

end

function id = get_compartment(x, y, map, pixelsize)
id = map(round(x / pixelsize) + (size(map, 1) + 1) / 2, round(y / pixelsize) + (size(map, 2) + 1) / 2);
end