function hop_fit_experimental_data()
% Hop trajectory fit of experimental data, uses the simulation model and
% Matlab function fmincon. The experimental data has already been fitted
% before with an anomalous
%
% Jan Keller-Findeisen (jan.keller@mpibpc.mpg.de)
%
% Part of software for "Cortical actin networks induce spatio-temporal
% confinement of phospholipids in the plasma membrane – a minimally
% invasive investigation by STED-FCS." by Andrade, D., Clausen, M.,
% Keller, J. et al. Sci Rep 5, 11454 (2015).

initialize();
fprintf('Fit experimental data to hop diffusion models.\n');

% experimental data of IA32 PE8
fwhms = [240, 190, 142, 90.9, 62.4, 50.7, 42]' * 1e-9; % (m)
Dcoeff = [0.382, 0.457, 0.487, 0.512, 0.581, 0.636, 0.737]; % µm/cm²

% boundaries and initial values for parameter
x0 = [Dcoeff(end), 100, 0.1];
lb = [0.3, 50, 0.001];
ub = [1.2, 400, 1];

opt = optimset('Display', 'iter', 'Algorithm', 'interior-point');
x = fmincon(@fitfunc, x0, [], [], [], [], lb, ub, [], opt);

    function y = fitfunc(x)
        dt = 20e-6;
        SR = 3e-6;
        Di = hop_simulation_run(1e2, 10, x(1) * 1e-12, dt, SR, x(2) * 1e-9, x(3), fwhms);
        y = sum((Di(:) - Dcoeff(:)).^2); % least squares minimization
    end

end