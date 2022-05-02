function hop_example()
% Examples of hopping FCS simulation
%
% Jan Keller-Findeisen (jan.keller@mpibpc.mpg.de)
%
% Part of software for "Cortical actin networks induce spatio-temporal
% confinement of phospholipids in the plasma membrane – a minimally
% invasive investigation by STED-FCS." by Andrade, D., Clausen, M.,
% Keller, J. et al. Sci Rep 5, 11454 (2015).

initialize();
fprintf('Run a hopping diffusion example.\n');

% producing some data
runSimulation('hop.example.mat');

end

function runSimulation(filename)
% performs one simulation run and stores results in a file

% parameters
fwhms = [20 : 20 : 100, 140 : 40 : 300]' * 1e-9; % FWHMS of gaussian PSFs [m]
Tp = 1e1; % total simulation time [s]
Np = 50; % number of molecules at the same time
Dp = 8e-13; % free diffusion coefficient [m^2/s]
dt = 50e-6; % simulation time step [s]
SR = 4e-6; % sand box radius

% hop parameter
Lp = 100e-9; % average size of compartments (sqrt of average area) [m]
Hp = 0.1; % hopping probability [0-1]

% number of repetitions of the whole simulation
repetitions = 3;

% loop
data = cell(repetitions, 2);
for ki = 1 : repetitions
    % call to simulation run (trapping simulation) with different seeds
    % each time (ki = 1, 2, 3, ..)
    [~, Di, TauD] = hop_simulation_run(Tp, Np, Dp, dt, SR, Lp, Hp, fwhms, ki);
    % store fitted diffusion coefficients and TauD times
    data{ki, 1} = Di;
    data{ki, 2} = TauD;
    % better to store after each repetition all the data obtained so far
    % (the file is overwritten)
    save(filename);
end

end