Water simulation started as a homework for a class. Combines the grid-based approach with particles.

Method used here closely follows the method from the paper ["Animating Sand as a Fluid" (2005)](http://dl.acm.org/citation.cfm?id=1073298). Particles are represented as small black circles, the velocities at each grid cell are lines pointing out of the cell with a size reletive to the magnatude, and pressure is visualized with blue for negative pressure and red for positive pressure.

Important
- This simulation requires freeglut or glut for proper compilation.
