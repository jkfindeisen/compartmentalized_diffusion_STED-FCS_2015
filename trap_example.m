function trap_example()
% Examples of trapping FCS simulation.
%
% Jan Keller-Findeisen (jan.keller@mpibpc.mpg.de)
%
% Part of software for "Cortical actin networks induce spatio-temporal
% confinement of phospholipids in the plasma membrane – a minimally
% invasive investigation by STED-FCS." by Andrade, D., Clausen, M.,
% Keller, J. et al. Sci Rep 5, 11454 (2015).

initialize();
fprintf('Run a trapping diffusion example.\n');

% producing some data
runSimulation('trap.example.mat');

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

% trap parameters
Trap.N = 50; % number of trepping centers in sandbox
Trap.D = 1e-14; % diffusion of trapping centers  [m^2/s]
Trap.R = 5e-9; % trapping radius around center
Trap.pon = 1; % probability for trappings (distance < Trap.R)
Trap.poff = 0.0005; % propabyility for release (distance  < Trap.R and trapped)

% number of repetitions of the whole simulation
repetitions = 3;

% loop
data = cell(repetitions, 2);
for ki = 1 : repetitions
    % call to simulation run (trapping simulation) with different seeds
    % each time (ki = 1, 2, 3, ..)
    [~, Di, TauD] = trap_simulation_run(Tp, Np, Dp, dt, SR, Trap, fwhms, ki);
    % store fitted diffusion coefficients and TauD times
    data{ki, 1} = Di;
    data{ki, 2} = TauD;
    % better to store after each repetition all the data obtained so far
    % (the file is overwritten)
    save(filename);
end

end