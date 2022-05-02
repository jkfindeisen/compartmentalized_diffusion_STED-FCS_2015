function hop_trajectory_example()
% An example showing you a hop trajectory.
%
% Jan Keller-Findeisen (jan.keller@mpibpc.mpg.de)
%
% Part of software for "Cortical actin networks induce spatio-temporal
% confinement of phospholipids in the plasma membrane – a minimally
% invasive investigation by STED-FCS." by Andrade, D., Clausen, M.,
% Keller, J. et al. Sci Rep 5, 11454 (2015).

initialize();
fprintf('Show an example hopping diffusion trajectory.\n');

% parameters
Tmax = 100;     % maximal time of the trajectory [s]
dt = 50e-6;     % time step [s]
L  = 10e-6;     % sandbox length [m]
dL = 20e-9;     % compartment map pixel size [m]
Df = 8e-13;     % free diffusion coefficient [m^2/s]
HL = 160e-9;    % average compartment diameter [m]
HP = 0.01;      % hopping probability [0-1]
seed = 0;       % random generator seed

% make trajectory
[trajectory, grid] = hop_trajectory(Tmax, dt, L, dL, Df, HL, HP, seed);

% show trajectory
h = figure;
plot([0, 0, L, L, 0], [0, L, L, 0, 0], 'k'); % rectangle (our sandbox)
hold on;

% show compartments
for kc = 1 : length(grid.C)
    v = grid.V(grid.C{kc}, :);
    if all(isfinite(v(:)))
        plot(v(:, 1), v(:, 2), 'b-');
        plot(v([end, 1], 1), v([end, 1], 2), 'b-');
    end
end

% show compartment centers
% plot(grid.centers(:, 1), grid.centers(:, 2), 'k.');

% show trajectory with certain time resolution
dtDisplay = 5e-3;     % display time step [s]
step = max(1, round(dtDisplay / dt));
p1 = trajectory(1, 2:3);
for kt = 1  + step : step : size(trajectory, 1)
    p2 = trajectory(kt, 2:3);
    plot([p1(1), p2(1)], [p1(2), p2(2)], 'r');
    p1 = p2;
end

text(trajectory(1, 2), trajectory(1, 3), 'start');
text(trajectory(end, 2), trajectory(end, 3), 'end');

% correct boundaries and everything
xlim([0, L]);
ylim([0, L]);
pbaspect([1 1 1]);
axis off;
title(sprintf('trajectory, compartment borders and centers [length %.1fs]', trajectory(end, 1)));
saveas(h, 'hop_trajectory_example.eps', 'epsc');

end

function [trajectory, grid] = hop_trajectory(Tmax, dt, L, dL, Df, HL, HP, seed)
% Simulates a hopping diffusion trajectory of a single molecule for a
% hopping diffusion model (i.e. free diffusion inside of compartments but
% changing from one compartment to the next only with a certain
% probability)
%
% Syntax
%   function [trajectory, grid] = hopping_trajectories(Tmax, dt, L, Df, HL, HP, seed)
%
% Input parameters
%   Tmax      - maximal total length [s]
%   dt        - time step [s]
%   L         - sandbox length [m]
%   dL        - compartment map pixel size [m]
%   Df        - free diffusion coefficient [m^2/s]
%   HL        - average compartment diameter/length [m]
%   HP        - hopping probability [0-1]
%   seed      - random generator seed (nonnegative integer)
%
% Output parameters
%   trajectory  - vector Nx4 with order (time [s], x [m], y [m], compartment [id])
%   grid        - characterization of the simulated
%
% Author
%   Jan Keller-Findeisen (jkeller1@gwdg.de)
%   Department of NanoBiophotonics
%   Max Planck Institute of Biophysical Chemistry
%   Am Fassberg 11, 3077 Goettingen, Germany
%
% Last change
%   2013-08

assert(nargin == 8, 'Wrong number of arguments!');

%% initializations
rng(seed);

% mean displacement in two directions (2D)
md = sqrt(4 * Df * dt);
if md > HL / 2
    warning('mean displacement (2D) is more than half of average compartment length, results might be compromised');
end
d = md / sqrt(2); % we need the mean displacement in 1D later

% number of compartments
num = round((L / HL)^2);

%% create hopping map
fprintf(' create hopping map\n');
centers = L * rand(num, 2);

% compute voronoi diagram
[V, C] = voronoin(centers);

% fill map with ids
[xm, ym] = ndgrid(0 : dL : L, 0 : dL : L);
map = zeros(size(xm));
id = 0;
nc = length(C);
for kc = 1 : nc
    v = V(C{kc}, :);
    if all(isfinite(v(:)))
        v = v([1:end, 1], :); % add first line again
        in = inpolygon(xm, ym, v(:, 1), v(:, 2));
        id = id + 1;
        map(in) = id;
    end
end

%% create hopping trajectory
fprintf(' simulate hopping diffusion\n');

% starting values, in the center of the map
x = L / 2;
y = L / 2;
id = map(floor(x / dL) + 1, floor(y / dL) + 1);

% the indefinitely long loop for each time step
store = cell(100, 1);
store_counter = 1;
buffer = zeros(1e5, 3);
buffer_counter = 1;
t = 0;
while 1
    % update position and compartment
    x0 = x + d * randn();
    y0 = y + d * randn();
    
    % outside of area, stop it here
    if  x0 < 0 || x0 > L || y0 < 0 || y0 > L
        break;
    end
    
    id0 = map(floor(x0 / dL) + 1, floor(y0 / dL) + 1);
    
    % different compartment and we do not hop
    if id ~= id0 && rand() > HP
        % retry until compartment is the same again
        while id0 ~= id
            x0 = x + d * randn();
            y0 = y + d * randn();
            
            id0 = map(floor(x0 / dL) + 1, floor(y0 / dL) + 1);
        end
    end
    
    % update position and compartment
    x = x0;
    y = y0;
    id = id0;
    
    % store position and compartment
    buffer(buffer_counter, :) = [x, y, id];
    buffer_counter = buffer_counter + 1;
    
    if buffer_counter > size(buffer, 1)
        store{store_counter} = buffer;
        store_counter = store_counter + 1;
        buffer = zeros(1e5, 3);
        buffer_counter = 1;
    end
    
    % increase timer
    t = t + dt;
    
    % time out, stop here
    if t > Tmax
        break;
    end
end

% trim buffer
buffer = buffer(1 : buffer_counter - 1, :);
store{store_counter} = buffer;

%% save
grid.map = map;
grid.V = V;
grid.C = C;
grid.centers = centers;
trajectory = cat(1, store{:});
trajectory = [(1 : size(trajectory, 1))' * dt, trajectory];

end