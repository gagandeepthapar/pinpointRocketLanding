# Pinpoint Rocket Landing

## Background
Trajectory generation of an arbitrary rocket following the path laid out in *Convex Optimization for Trajectory Optimization: A Tutorial on Generating Dynamically Feasible Trajectories Reliably and Efficiently*, D. Malyuta et al.
 
## Approach and Misc. Notes
Dynamics of the rocket were encoded as individual constraints at each time step constrained by the amount of propellant onboard the vehicle. A standard state propagator was used, dissimilar to the approach in the paper, to define how the state changes over time.

A feasible trajectory was generated using simplified dynamics of the rocket and mass modeling. This trajectory was found to be similar to that in the paper, although requires additional work to achieve full reconciliation.

![feasible_3d](media/feasible_3d.png)
![feasible_2d](media/feasible_downrange.png)
![feasible_dim](media/feasible_traj.png)

Some interesting trajectories were generated such as a subsurface flight (violation of the glideslope constraints)

![subsurface_2d](media/subsurface_2dlanding.png)
![subsurface_3d](media/subsurface_3dlanding.png)
![subsurface_dims](media/subsurface_dimensionlanding.png)

and a technically feasible trajectory that approaches the ground vertically instead of following the glideslope closer. This is very clearly a suboptimal trajectory.

![feasible_2d](media/2dlanding.png)
![feasible_3d](media/3dlanding.png)
![feasible_dim](media/dimensionlanding.png)

Several suboptimal solutions were found (data can be found in **mat_files** directory). While the focus of the paper was Convex Optimization and Second Order Cone Programming, MATLAB is more suited to use *fmincon*, another approach to minimizing functions. *fmincon* is substantially slower than other SOCP solvers and, consequentially, takes ~30-45min to generate a 150 second trajectory given the constraints.

## Further Reading
A quick [presentation deck](https://github.com/gagandeepthapar/pinpointRocketLanding/blob/main/RocketLanding_PPT.pdf) and [summary](https://github.com/gagandeepthapar/pinpointRocketLanding/blob/main/RocketLanding_Summary.pdf) of the project can be found in this repository as well.

## Future Work
Future work includes modifying constraints and possibly porting the code to a lighter language e.g., Python to improve performance. It should be noted that the SOCP solver and constraint handler used in the paper, [ECOS](https://github.com/embotech/ecos) and [JuMP](https://github.com/jump-dev/JuMP.jl) respectively, were not used as there was no obvious support for MATLAB. ECOS has better Python support, although JuMP appears to be [Julia](https://julialang.org/) specific.
