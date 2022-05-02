function [x, y, t, xt, yt] = trap_simulation_trajectory(D, R, T, dt, x0, y0, Trap)
% Simulates a trajectory in the presence of tracking.
%
%   D - free diffusion coefficient [m^2/s]
%   R - sandbox radius [m]
%   T - total time of rajectory [s]
%   dt - time step [s]
%   x0, y0 - starting point [m]
%   Trap - trapping structure (see trap_simulation_run.m)
%
% Author: Jan Keller-Findeisen (jan.keller@mpibpc.mpg.de)
% Dep. NanoBiophotonics, MPI for Biophysical Chemistry, Goettingen, Germany
% Date: 2013/06-2014/01

assert(nargin == 7, 'Not enough arguments!');

% parameters
dr = sqrt(2 * D * dt); % average displacement (free diffusion)
Trap.dr = sqrt(2 * Trap.D * dt); % average displacement trapping centers
Nt = round(T / dt);
R2 = R^2;
Trap.R2 = Trap.R^2;

% random values
sx = dr * randn(Nt, 1);
sy = dr * randn(Nt, 1);
Trap.sx = Trap.dr * randn(Trap.N, Nt);
Trap.sy = Trap.dr * randn(Trap.N, Nt);

% output values
x = zeros(Nt, 1);
x(1) = x0;
y = zeros(Nt, 1);
y(1) = y0;
t = zeros(Nt, 1);
xt = zeros(Trap.N, Nt);
yt = zeros(Trap.N, Nt);
[xt(:, 1), yt(:, 1)] = random_position_on_disc(R, Trap.N);

% loop over time steps from second to last
for ki = 2 : Nt
    
    % diffuse traps
    xt(:, ki) = xt(:, ki - 1) + Trap.sx(:, ki);
    yt(:, ki) = yt(:, ki - 1) + Trap.sy(:, ki);
    % out of border check
    o = xt(:, ki).^2 + yt(:, ki).^2 > R2;
    phi = atan2(yt(o, ki), xt(o, ki));
    xt(o, ki) = xt(o, ki) - 2 * R * cos(phi);
    yt(o, ki) = yt(o, ki) - 2 * R * sin(phi);
    
    % is it trapped?
    if t(ki - 1) > 0
        % it is, try to untrap
        if rand() < Trap.poff
            % no need to set t(ki) to zero, is by default
            % diffuse freely
            x(ki) = x(ki - 1) + sx(ki);
            y(ki) = y(ki - 1) + sy(ki);
            % out of border check
            if x(ki)^2 + y(ki)^2 > R2
                phi = atan2(y(ki), x(ki));
                x(ki) = x(ki) - 2 * R * cos(phi);
                y(ki) = y(ki) - 2 * R * sin(phi);
            end
        else
            % it remains trapped, move with trap instead
            t(ki) = t(ki - 1);
            x(ki) = xt(t(ki), ki);
            y(ki) = yt(t(ki), ki);
        end
    else
        % it is not trapped, look if it could be trapped
        [dm, idx] = min((xt(:, ki) - x(ki - 1)).^2 + (yt(:, ki) - y(ki - 1)).^2);
        if dm < Trap.R2 && rand() < Trap.pon
            % it got trapped, set to traps position
            t(ki) = idx;
            x(ki) = xt(idx, ki);
            y(ki) = yt(idx, ki);
        else
            % it did not get trapped for whatever reason, move freely
            x(ki) = x(ki - 1) + sx(ki);
            y(ki) = y(ki - 1) + sy(ki);
            % out of border check
            if x(ki)^2 + y(ki)^2 > R2
                phi = atan2(y(ki), x(ki));
                x(ki) = x(ki) - 2 * R * cos(phi);
                y(ki) = y(ki) - 2 * R * sin(phi);
            end
        end
    end
end

end