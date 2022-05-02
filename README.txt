Example software package for 

"Cortical actin networks induce spatio-temporal confinement of phospholipids in the plasma membrane
– a minimally invasive investigation by STED-FCS."

by Andrade, D., Clausen, M., Keller, J. et al. Sci Rep 5, 11454 (2015).

https://doi.org/10.1038/srep11454

Note: Contains a fast FCS curve calculation, which worked with Matlab-2014b, but does not work anymore with Matlab-2021b.
It might be simple to fix, but we cannot fix it right now.

1. Examples of free, trapped and hopping diffusion

Call free_example.m, hop_example.m and trap_example.m. Will compute correlation curves (on an old Matlab version).

2. Free, trappend and hopping diffusion simulations

See code/free_simulation_trajectory.m, hop_simulation_trajectory.m, trap_simulation_trajectory.m

3. FCS curve computation

See code/free_simulation_run, hop_simulation_run.m, trap_simulation_run.m and in particular fcs_corr_xxx.m.

Note: Requires compiled C code. Will not run on newer Matlabs (old Matlab 2014 will probably work), newer will not.

4. Fitting of FCS curves to experimental data

See hop_fit_experimental_data.m and code/free_simulation_run, hop_simulation_run.m, trap_simulation_run.m and there in particular subfunction anomFCSCurve, the curve

License: MIT license (see LICENSE.txt) except for ui_progressbar.m (2-clause BSD)
Code is provided as is, no guarantees for anything.

