function [Ds, Di, TauD] = hop_simulation_run(Tp, Np, Dp, dt, SR, Lp, Hp, fwhms, seed)
% FCS simulation run for hopped diffusion
%
% Simulates trajectories and photon traces, computes correlations, fits
% anomalous diffusion model.
%
% Syntax
%   function [Ds, Di, TauD] = hop_simulation_run(Tp, Np, Dp, dt, Lp, Hp, fwhms, seed)
%
% Input parameters
%   Tp - total time [s]
%   Np - number of molecules in the sandbox volume simultaneously
%   Dp - Free diffusion coefficient [m^2/s]
%   dt - simulation time step [s]
%   SR - sand box radius [m]
%   Lp - Average size of compartments (sqrt of average area) [m]
%   Hp - Hopping probability [0-1]
%   fwhms - Vector of all FWHMs of the microscope that you used [m]
%           Gaussian radially symmetric 2D
%   seed - random number generator seed value
%
% Output parameters
%   Ds - mean of Di in second dimension (averaging over all foci in the sandbox)
%        that means a single diffusion time for each FWHM of the PSF
%   Di apparent diffusion times [µm^2/s]
%   TauD apparent tau d times [s]
%
% Comment:
% The comparments are created by throwing random points into the sandbox
% (seed points for a Voronoi mesh) and creating a mesh around them so that
% irregular (randomly) shaped polgonal areas are created. Then a hopping
% map is constructed (a big matrix) where for each cell the corresponding
% compartment (they are numbered) is written, i.e. the compartments are
% rastered.
%
% For large FWHMs the apparent diffusion coefficient might go up slightly
% because of finite size of sandbox. Try with larger sandbox values and
% larger distance between neighbored foci. Also make sure the time step is
% small enough to not jump over two compartments!
%
% Be aware that too large time step (mean displacement should be smaller
% than smallest FWHM), too small sandbox volume, not enough repetitions,
% can lead to deviations from the inserted free diffusion coefficient Dp in
% the outcome. This is a statistics effect with random and systematic
% errors.
%
% Author: Jan Keller-Findeisen (jan.keller@mpibpc.mpg.de)
% Dep. NanoBiophotonics, MPI for Biophysical Chemistry, Goettingen, Germany
% Date: 2013/06-2013/08

assert(nargin == 9, 'Wrong number of parameter!');

close all; % close all figures that are open

%% simulation paremeters

rng(seed);     % random number seed
% rng('shuffle');

% save physical parameters
sim.T = Tp;                % Trace duration [s]
sim.time_bin = Tp;           % Trace (compute traces in this time bin)
sim.N = Np;                 % Number of particles (lipids, ..)
sim.R = SR;               % sandbox radius [m]
sim.D = Dp;                % diffusion coefficient [m^2s^-1]
sim.dt = dt;            % time step dt (default: 100 µs) % Pixel dwell time [s]
sim.mean_displacement = sqrt(2 * sim.D * sim.dt); % in one direction (1D)
fprintf(' mean displacement during one step (free diffusion) %.2fnm\n', sim.mean_displacement * sqrt(2) / 1e-9); % in 2D

% additonal hop mode parameters
sim.hop.R = 1.1 * sim.R + Lp; % hopping map box radius (is a bit larger than the sandbox in order to catch particles that diffuse outside)
sim.hop.N = round(4 * (sim.hop.R / Lp)^2); % number of seed points needed to achieve exactly the desired length (square root of average area)
sim.hop.pixel_size = 4e-9;  % pixel size of hoping map [m]

% Diffusion Mode
correlator.nbr_octaves = 20;
correlator.nbr_bins_per_octaves = 5;

% excitation
psf.FWHMs = fwhms;
psf.number_FWHMs = length(psf.FWHMs);
psf.R = 2 * sim.R; % calculation distance for the PSF calculation
psf.pixel_size = 2e-9; % step size on which we calculate the PSF [m]
psf.psf_distance = 600e-9; % distance between two neighbored foci
psf.Nr = 3; % results in (2*Nr+1)^2 foci

if sqrt(2) * psf.psf_distance * psf.Nr > sim.R
    error('PSF foci outside of sandbox radius, increase sandbox');
end

% calculate excitation psf
psf.ri = ndgrid(0:psf.pixel_size:psf.R);

[psf.cxi, psf.cyi] = ndgrid((-psf.Nr:psf.Nr) * psf.psf_distance, (-psf.Nr:psf.Nr) * psf.psf_distance);
psf.number_centers = numel(psf.cxi);

psf.V = cell(psf.number_FWHMs, 1);
for ki = 1 : psf.number_FWHMs
    psf.V{ki} = power(2., - psf.ri.^2 / (psf.FWHMs(ki) / 2)^2); % gaussian distribution
end

% fit parameters
fit.nbr_component = 1; % one component
% order: N, tauD, alpha, DC
fit.x0 = [0, 0, 1, 0]; % alpha0 and DC0
fit.lb = [0, 0, 0.2, -0.05];
fit.ub = [Inf, Inf, 1.5, 0.05];

%% simulation
fprintf(' calculate compartment distribution\n');

%% make hoping map

% centers of compartments
L = sim.hop.R;
% randomly oriented seeds
seeds = L * (rand(sim.hop.N, 2) - 0.5) * 2;
% nn = omex_nearest_neighbour(seeds);
% fprintf('mean distance to nearest neighbor of compartment centers %.2f nm\n', mean(nn(:, 1)) / 1e-9);
fprintf(' average sqrt(area) %.2f nm\n', (L * 2) / sqrt(sim.hop.N) / 1e-9);
% compute voronoi diagram
[V, C] = voronoin(seeds);
% make figure
%h = figure;
phi = 0: 0.1 : 2 * pi + 0.1;
%plot(sin(phi) * sim.R, cos(phi) * sim.R, 'k');
% hold on;
for kc = 1 : length(C)
    v = V(C{kc}, :);
    if all(isfinite(v(:)))
        %plot(v(:, 1), v(:, 2), 'b-');
        %plot(v([end, 1], 1), v([end, 1], 2), 'b-');
    else
        
    end
end

% plot(seeds(:, 1), seeds(:, 2), 'k.');

W = 1.2 * sim.R;
%         xlim([-W, W]);
%         ylim([-W, W]);
%         pbaspect([1 1 1])
%         axis off;

% saveas(h, [output_name,'.hopping.map.eps'], 'epsc');
% close(h);

% determine if rastered (hopping map) pixel is in polgyon (time consuming
% and not very elegant, but working)
L = 1.2 * sim.R;
[xm, ym] = ndgrid(-L : sim.hop.pixel_size : L, -L : sim.hop.pixel_size : L);
map = zeros(size(xm));
id = 0;
ui_progressbar(0.);
for kc = 1 : length(C)
    v = V(C{kc}, :);
    if all(isfinite(v(:)))
        v = v([1:end, 1], :); % add first line again
        in = inpolygon(xm, ym, v(:, 1), v(:, 2));
        id = id + 1;
        map(in) = id;
    end
    ui_progressbar(kc / length(C));
end
ui_progressbar(1.);
sim.hop.map = map;


Nsim = ceil(sim.T / sim.time_bin); % number of simulation steps
%% simulate trajectories and photon traces
fprintf(' simulate photon traces\n');

% for all molecules (we add up over all molecules)
trace = cell(psf.number_FWHMs, psf.number_centers);
[trace{:}] = deal(0); % set all to zero
ui_progressbar(0.);
for kn = 1 : sim.N
    
    % place molecule randomly within sandbox
    [x0, y0] = random_position_on_disc(sim.R, 1);
    
    % for all simulation steps
    trajectory = cell(Nsim, 1);
    for ks = 1 : Nsim
        
        % T determines actual bin size (usually sim.time_bin, except for the last bin which can be smaller
        T = min(sim.time_bin, sim.T - sim.time_bin * (ks - 1));
        
        % number of samples in this step
        Nsam = round(T / sim.dt); % number of samples
        
        % get trajectory
        [x, y] = hop_simulation_trajectory(sim.mean_displacement, sim.R, Nsam, x0, y0, sim.hop.map, sim.hop.pixel_size, Hp);
        x0 = x(end);
        y0 = y(end);
        trajectory{ks} = [x, y];
    end
    
    % merge trajectory
    trajectory = cat(1, trajectory{:});
    x = trajectory(:, 1);
    y = trajectory(:, 2);
    
    % compute photon traces
    
    % for each psf position
    for kj = 1 : psf.number_centers
        
        % get radius in pixel of psf
        ri = sqrt((x - psf.cxi(kj)).^2 + (y - psf.cyi(kj)).^2);
        ri = round(ri / psf.pixel_size) + 1;
        
        % for each FWHM
        for ki = 1 : psf.number_FWHMs
            V = psf.V{ki};
            trace{ki, kj} = trace{ki, kj} + V(ri); % add all molecules
        end
    end
    ui_progressbar(kn / sim.N);
end
ui_progressbar(1.);

%% correlation
fprintf(' compute correlation curve\n');
correlator.nbr_octaves = floor(log(numel(trace{1})) + 2);
correlator.nbr_octaves = 15;

% do the correlations
fcs_curve = cell(psf.number_FWHMs, psf.number_centers);
ui_progressbar(0.);
for kn = 1 : numel(trace);
    % fprintf('Run %d of %d ...\n', kn, psf.number_FWHMs);
    handle = fcs_corr_init('type', 'multi', 'n_compression', 1, 'f_time', 1, 'n_octaves', correlator.nbr_octaves,'n_bins_per_octave', correlator.nbr_bins_per_octaves);
    fcs_corr_add(handle, trace{kn}'); % needs line vector
    [tau, fcs_curve{kn}] =  fcs_corr_curve(handle, 'n_norm', 2);
    fcs_corr_free(handle);
    ui_progressbar(kn / numel(trace));
end
ui_progressbar(1.);
tau = tau * sim.dt;

% normalize by the mean first three values of the average curve for each FWHM
normalization = zeros(psf.number_FWHMs, 1);
fcs_average = cell(psf.number_FWHMs, 1);
for ki = 1 : psf.number_FWHMs
    
    avg_curve = 0;
    for kj = 1 : psf.number_centers
        avg_curve = avg_curve + fcs_curve{ki, kj};
    end
    avg_curve = avg_curve / psf.number_centers;
    
    % normalization(ki) = mean(avg_curve(1:3)); % not sure what is better here
    normalization(ki) = avg_curve(1);
    fcs_average{ki} = avg_curve / normalization(ki);
    
    for kj = 1 : psf.number_centers
        fcs_curve{ki, kj} = fcs_curve{ki, kj} / normalization(ki);
    end
end

%% fit the fcs curves
fprintf(' fit correlation curves\n');
fits.N = zeros(psf.number_FWHMs, psf.number_centers);
fits.TauD = zeros(psf.number_FWHMs, psf.number_centers);
fits.alpha = zeros(psf.number_FWHMs, psf.number_centers);
fits.curves = cell(psf.number_FWHMs, psf.number_centers);

% opt = optimset('Display', 'iter', 'Algorithm', 'interior-point'); % if you want to see what is going on
opt = optimset('Display', 'off', 'Algorithm', 'interior-point');

% for all FWHMs
ui_progressbar(0.);
for ki = 1 : psf.number_FWHMs
    
    average = fcs_average{ki};
    
    % meaningful starting values
    fit.x0(1) = 1 / average(1);
    fit.x0(2) = psf.FWHMs(ki)^2 / (8 * log(2) * sim.D); % is this a good starting value
    fit.x0(end) = mean(average(end-3:end));
    fit.lb(end) = -0.05 * average(1);
    fit.ub(end) = 0.05 * average(1);
    
    % for all PSF centers
    for kj = 1 : psf.number_centers
        
        data = fcs_curve{ki, kj};
        fit.x0(1) = 1 / data(1); % start without the average
        % do the fitting
        [x, ~] = fmincon(@minizer, fit.x0, [], [], [], [], fit.lb, fit.ub, [], opt);
        
        fits.curves{ki, kj} = anomFCSCurve(x, tau, fit.nbr_component);
        fits.N(ki, kj) = x(1);
        fits.TauD(ki, kj) = x(2);
        fits.alpha(ki, kj) = x(3);
    end
    
    ui_progressbar(ki / psf.number_FWHMs);
end
ui_progressbar(1.);

    function y = minizer(x)
        G = anomFCSCurve(x, tau, fit.nbr_component);
        y = (G - data).^2;
        % y = y ./ sqrt(sqrt((1 : length(y))'));
        y = mean(y(:)) * 1e6;
    end


%% plot results
figure

% plot correlation curve and residuals
for kn = 1 : psf.number_FWHMs
    subplot(4, 2, [3,5]);
    semilogx(tau, fits.curves{kn, 1}, 'r-'); % fitted
    hold on;
    semilogx(tau, fcs_curve{kn, 1},'k.'); % real data
    ylim([-0.02, 1.05]);
    xlabel('Tau[s]')
    ylabel('G(Tau)');
    
    subplot(4,2,7);
    semilogx(tau, fcs_curve{kn, 1} - fits.curves{kn, 1}, '-r');
    hold on;
    xlabel('Tau[s]')
    ylabel('Residuals');
end

% plot alpha values
subplot(4, 2, 2);
fwhms = psf.FWHMs / 1e-9;
mAlpha = mean(fits.alpha, 2);
eAlpha = std(fits.alpha, 0, 2) / sqrt(size(fits.alpha, 2));
% plot(fwhms, mean(fits.alpha, 2), 'o--b');
errorbar(fwhms, mAlpha, eAlpha, 'o-b');
xlim([min(fwhms) - 20, max(fwhms) + 20]);
xlabel('FWHM [nm]');
ylabel('alpha');

% plot TauD values
subplot(4, 2, 4);
mTauD = mean(fits.TauD, 2) / 1e-3; % [ms]
eTauD = std(fits.TauD, 0, 2) / sqrt(size(fits.TauD, 2));
errorbar(fwhms.^2, mTauD, eTauD, 'o-b');
% xlim([min(fwhms) - 20, max(fwhms) + 20]);
xlabel('FWHM^2 [nm^2]');
ylabel('TauD [ms]');

% plot D
subplot(4, 2, 6);
Di = (repmat(psf.FWHMs, 1, psf.number_centers) / 1e-6).^2 ./ (8 * log(2) * fits.TauD); % [µm^2/s]
TauD = fits.TauD;
Ds = mean(Di, 2);
eD = std(Di, 0, 2) / sqrt(size(Di, 2));
errorbar(fwhms, Ds, eD, 'o-b');
xlim([min(fwhms) - 20, max(fwhms) + 20]);
xlabel('FWHM [nm]');
ylabel('D [µm^2/s]');

% plot example trace
% example_time = (1 : 1000 : numel(trace{1})) * sim.dt;
% example_trace = trace(1 : 1000 : end, :);
%
% subplot(4, 2, 1);
% plot(example_time, example_trace);
% xlabel('Time[s]');
% ylabel('Photon Counts');
%
% % plot also in table
% cnames={'Tau [ms]','D [1e13]','N','alpha'};
% rnames = cell(psf.number_FWHMs, 1);
% for ki = 1 : psf.number_FWHMs
%     rnames{ki} = sprintf('%d nm', round(psf.FWHMs(ki) / 1e-9));
% end
%
% data = [fits.TauD * 1e3, psf.FWHMs.^2 ./ 4 ./ fitresults.TauD * 1e13, fitresults.N, fitresults.alpha];
% uitable('Parent', gcf, 'Data', data, 'ColumnName', cnames, 'RowName', rnames, 'Units', 'normalized', 'Position', [0.53,0.07,0.4,0.2]);

% saveas(gcf, [output_name,'_Figure.fig'], 'fig');
end

%%
function G = anomFCSCurve(fitparameter, tau, component)
% generates FCS theory curve for fitting G(tau) = FCScurve(a, tau, TauD, s)
% N = vector of amplitudes (column)
% tau = time base (row vector)
% TauD = vector of diffusion times with amplitudes a (column)
% s = structure parameter (scalar)

N = fitparameter(1:component);
TauD = fitparameter(component + 1 : 2 * component);
alpha = fitparameter(end-1);
DC = fitparameter(end);

a = 1 ./ N;
td = 1 ./ TauD;   %for outer product multiplication

G = zeros(size(tau));
for ki = 1 : length(tau)
    g = 1 ./ (1 + (td .* tau(ki)).^alpha);
    G(ki) = a * g' + DC;
end

end