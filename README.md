# Blade-element-theory-solver
Solver for computing propeller performance characteristics using Blade Element Momentum Theory, which is the combination of the tradiational blade element theory and the actuator disc theory. This solver can be used as a first order propeller design tool and was developed as an auxilliary project to my undergraduate thesis project.
The solver takes the propeller geometric characteristics (sectional chord length, radial station distance from the centre, hub diameter, sectional geometrical pitch, sectional airfoil) as input. The code for 2 propellers has been provided:

1). APC 11x4.7, which is a market standard propeller by the and has a diameter of 11 inches and a representative section pitch of 4.7 inches. The inward airfoil is Eppeler 63 and the outward airfoil is Clark Y.

2). A self-designed propeller having the ARAD-10 as the main airfoil and MH 32 8.7% as the tip airfoil. This propeller was the final propeller which was the result of the 1st order design in our undergraduate thesis project and according to the 1st order tool QPROP (developed by Dr. Mark Drela), showed better efficiency and higher thrust for low advance ratios than APC 11x4.7.

This solver is validated by comparing the QPROP results and the solver results for APC 11x4.7. Technically, QPROP is a much more sophisticated tool compared to this solver as QPROP makes use of the vortex theory which considers 3D flow of air around the propeller, while the blade element momentum theory is just a 2D concept (where the axial and lateral air flow is considered but the air flow across the radius of propeller blade is ignored).

A lot of sophistication and automation of this solver can still be done (for eg. the reynolds number or the flight velocity range still cant be controlled by the user)
