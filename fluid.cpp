#include "fluid.hpp"

#include "opengl.hpp"
#include "sparse.hpp"
#include <iostream>
#include <cstdio>
#include <fstream>

using namespace std;


int framecount = 0;

// random double uniformly distributed between 0 and 1
/// double urand() {
    /// return (double)rand()/RAND_MAX;
/// }

void Fluid::init(int m, int n, Vector2d x0, double dx) {
    this->m = m;
    this->n = n;
    this->x0 = x0;
    this->dx = dx;

    showPressure = showVelocity = false;
    particleSize = 10;

    Vector2d xGrid = x0 + Vector2d(dx/2,dx/2); // center of grid cell (0,0)
    vel = StaggeredGrid(m, n, xGrid, dx);
    velOld = StaggeredGrid(m, n, xGrid, dx);
    velWeight = StaggeredGrid(m, n, xGrid, dx);
    rho = StaggeredGrid(m, n, xGrid, dx);
    rhoWeights = StaggeredGrid(m, n, xGrid, dx);
    flags = Grid<char>(m, n, xGrid, dx);
    indices = Grid<int>(m, n, xGrid, dx);
    pressure = Grid<double>(m, n, xGrid, dx);
    vel.assign(Vector2d(0,0));
    pressure.assign(0);
    // create particles
    for (int i = 0; i < m/2; i++) {
        for (int j = 0; j < 3*n/4; j++) {
            Vector2d xij = x0 + Vector2d(i, j)*dx;
            for (int i1 = 0; i1 < 4; i1++) {
                for (int j1 = 0; j1 < 4; j1++) {
                    Vector2d x = xij + Vector2d(i1, j1)*(dx/2)
                        + Vector2d(urand(), urand())*(dx/2);
                    particles.push_back(new Particle(x, Vector2d(0,0), Vector3d(0, 0, 1), 1));
                }
            }
        }
    }
    for(int i = 3*m/4; i < m; i++) {
        for(int j = 0; j < n/2; j++) {
            Vector2d xij = x0 + Vector2d(i, j)*dx;
            for(int i1 = 0; i1 < 4; i1++) {
                for(int j1 = 0; j1 < 4; j1++) {
                    Vector2d x = xij + Vector2d(i1, j1)*(dx/2) + Vector2d(urand(), urand())*(dx/2);
                    particles.push_back(new Particle(x, Vector2d(0, 0), Vector3d(0, 0, 0), 0.5));
                }
            }
        }
    }
    forces.push_back(new GravityForce(Vector2d(0.0, -9.8)));
}

void Fluid::draw() {
    double arrowScale = 1;
    double pressureScale = 1;
    if(showPressure)
    {
        glBegin(GL_QUADS);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                Vector2d x = x0 + Vector2d(i, j)*dx;
                double p = pressureScale*pressure(i,j);
                glColor3f(1 + p, 1 - abs(p), 1 - p);
                glVertex2f(x(0), x(1));
                glVertex2f(x(0)+dx, x(1));
                glVertex2f(x(0)+dx, x(1)+dx);
                glVertex2f(x(0), x(1)+dx);
            }
        }
        glEnd();
    }
    glColor3f(0.8, 0.8, 0.8);
    glPointSize(10);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            Vector2d x = x0 + Vector2d(i+0.5, j+0.5)*dx;
            Vector2d v = vel.interpolate(x);
            glBegin(GL_POINTS);
            glVertex2f(x(0), x(1));
            glEnd();
            if(showVelocity)
            {
                glBegin(GL_LINES);
                glVertex2f(x(0), x(1));
                glVertex2f(x(0) + arrowScale*v(0), x(1) + arrowScale*v(1));
                glEnd();
            }
        }
    }
    
    glPointSize(particleSize);
    glBegin(GL_POINTS);
    for (int p = 0; p < particles.size(); p++) {
        Particle *par = particles[p];
        glColor3f(par->color(0), par->color(1), par->color(2));
        Vector2d x = particles[p]->x;
        glVertex2f(x(0), x(1));
    }
    glEnd();
}

void Fluid::step(double dt) {
    advection(dt);
    particlesToGrid();
    velOld = vel;
    flagCells();
    addExternalForces(dt);
    applyBoundaryConditions();
    pressureSolve(dt);
    gridToParticles();
}

void Fluid::advection(double dt) {
    // TODO: for each particle, move it according to the 'vel' grid
    // (not its own velocity).
    // Use 2nd-order Runge-Kutta, i.e. move it according to u(x + u(x) dt/2).
    // If it ends up outside the grid, project it back.
    for(int i = 0; i < particles.size(); i++) {
		Particle* p = particles[i];
		Vector2d newpos = p->x + vel.interpolate(p->x + 0.5 * vel.interpolate(p->x) * dt) * dt;
		if(newpos(0) < x0(0)) {
			newpos(0) = x0(0) + 0.0001; // add a bit so it doesn't get stuck due to error
		}
		else if(newpos(0) > x0(0)+dx*m) {
			newpos(0) = x0(0) + dx*m - 0.0001;
		}
		if(newpos(1) < x0(1)) {
			newpos(1) = x0(1) + 0.0001;
		}
		else if(newpos(1) > x0(1)+dx*n) {
			newpos(1) = x0(1) + dx*n - 0.0001;
		}
		p->x = newpos;
	}
}

void Fluid::particlesToGrid() {
    vel.assign(Vector2d(0,0));
    velWeight.assign(Vector2d(1e-3,1e-3)); // avoid division by zero
    // TODO: for each particle, add its velocity to its nearby grid cells.
    // Use velWeight to determine the weighted average denominator
    for(int i = 0; i < particles.size(); i++) {
		Particle* p = particles[i];
		vel.addInterpolated(p->x, p->v);
		velWeight.addInterpolated(p->x, Vector2d(1, 1));
	}
	for(int i = 0; i < vel.u.m; i++) {
		for(int j = 0; j < vel.u.n; j++) {
			vel.u(i, j) /= velWeight.u.get(i, j);
		}
	}
	for(int i = 0; i < vel.v.m; i++) {
		for(int j = 0; j < vel.v.n; j++) {
			vel.v(i, j) /= velWeight.v.get(i, j);
		}
	}
}

void Fluid::flagCells() {
    flags.assign(AIR); // this will change if obstacles or free surfaces
    for(int i = 0; i < particles.size(); i++) {
        Particle* p = particles[i];
        Vector2d relP = p->x - x0;
        //printf("P: (%f, %f), x0: (%f, %f)\n", p->x(0), p->x(1), x0(0), x0(1));
        int x = floor(relP(0)/dx);
        int y = floor(relP(1)/dx);
        //printf("Rel Pos: (%f, %f), (%d, %d)\n", relP(0), relP(1), x, y);
        flags(x, y) = FLUID;
    }
}

void Fluid::addExternalForces(double dt) {
    for (int f = 0; f < forces.size(); f++)
        forces[f]->addForces(this, dt);
}

void Fluid::applyBoundaryConditions() {
    // TODO: zero out normal velocities at obstacle boundaries.
    // (Right now, the only obstacle boundaries are the ends of the grid.)

	for(int i = 0; i < n; i++) {
		vel.u.get(0, i) = 0;
		vel.u.get(n, i) = 0;
	}

	for(int i = 0; i < m; i++) {
		vel.v.get(i, 0) = 0;
		vel.v.get(i, m) = 0;
	}
}

void Fluid::pressureSolve(double dt) {
	indices.assign(-1.0);
    pressure.assign(0);
    addDensity();
    // 1. Enumerate cells flagged as fluid, store in 'indices' grid
    int num = 0;
    for(int i = 0; i < m; i++) {
		for(int j = 0; j < n; j++) {
			if(flags.get(i, j) == FLUID) {
				indices(i, j) = num;
				num++;
			}
		}
	}
    
    // 2. Build a linear system, solve it, put results into 'pressure' grid
    SpMatrix A(num, num);
    VectorXd b = VectorXd(num);
    double scale = 1.0 / (dx * dx);
    for(int i = 0; i < m; i++) {
		for(int j = 0; j < n; j++) {
			int cellIndex = indices.get(i, j);
			if(cellIndex != -1) { //fluid
				// Build A ( Del^2 p)
				double center = 0.0;
				double d;               //density at the face between the cells
                //Check if any neighboring cells are outside boundary or are obstacles
                if(i+1 < m && flags.get(i+1, j) != OBS) {
                    d = rho.u(i+1, j);
                    center += -1.0/d;
				}
				if(j+1 < n && flags.get(i, j+1) != OBS) {
                    d = rho.v(i, j+1);
                    center += -1.0/d;
				}
				if(i-1 > -1 && flags.get(i-1, j) != OBS) {
                    d = rho.u(i, j);
                    center += -1.0/d;
				}
				if(j-1 > -1 && flags.get(i, j-1) != OBS) {
                    d = rho.v(i, j);
                    center += -1.0/d;
				}
                if(isinf(center)) {
                    printf("Center Infinity\n");
                }
                center *= scale * dt;
                if(isinf(center)) {
                    printf("Center Infinity: %f, %f\n", scale, dt);
                    std::exit(0);
                }
                A.add(cellIndex, cellIndex, center);
				//Add values for neighboring cells
				if(i+1 < m && flags.get(i+1, j) == FLUID) {
                    d = rho.u(i+1, j);
                    //center -= -1.0/d;
					A.add(cellIndex, indices.get(i+1, j), scale*dt/d);
				}
				if(j+1 < n && flags.get(i, j+1) == FLUID) {
                    d = rho.v(i, j+1);
                    //center -= 1.0/d;
					A.add(cellIndex, indices.get(i, j+1), scale*dt/d);
				}
				if(i-1 > -1 && flags.get(i-1, j) == FLUID) {
                    d = rho.u(i, j);
                    //center -= 1.0/d;
					A.add(cellIndex, indices.get(i-1, j), scale*dt/d);
				}
				if(j-1 > -1 && flags.get(i, j-1) == FLUID) {
                    d = rho.v(i, j);
                    //center -= 1.0/d;
					A.add(cellIndex, indices.get(i, j-1), scale*dt/d);
				}
                
				//Build b (Del . u)
				//Del . u = (u1 - u0)/dx + (v1 - v0)/dx
				double du = (vel.u.get(i+1, j) - vel.u.get(i, j)) / dx;
				double dv = (vel.v.get(i, j+1) - vel.v.get(i, j)) / dx;
				b(cellIndex) = (du + dv);
			}
		}
	}
    //std::ofstream out("Matrix.txt");
    //out << "A:\n";
    //SparseMatrix<double> print(num, num);
    //print.setFromTriplets(A.triplets.begin(), A.triplets.end());
    //out << print;
    //out << "\n\nP:\n";
    //out.close();
    //std::exit(0);
    //Solve
    VectorXd p = A.solve(b);
    //Put in pressure grid
    //bool stop = false;
    for(int i = 0; i < m; i++) {
		for(int j = 0; j < n; j++) {
			int index = indices.get(i, j);
			if(index != -1) {
				pressure(i, j) = p(index);
                //if(isnan(pressure(i, j))) {
                //    stop = true;
                //}
                //out << pressure(i, j) << "\n";
			}
		}
	}
    //out.close();
    //if(stop)
    //    std::exit(0);
    // 3. Apply pressure forces to the velocity field 'vel'
	//u = u - dt * Del p / rho
    for(int i = 0; i < m; i++) {
		for(int j = 0; j < n; j++) {
            if(i+1 < m) {
                double Dpu = (pressure.get(i+1, j) - pressure.get(i, j)) / dx;
                double oldu = vel.u(i+1, j);
                if(rho.u(i+1, j) > 1e-6) {
                    vel.u(i+1, j) = vel.u(i+1, j) - dt * Dpu / rho.u(i+1, j);
                    if(isnan(vel.u(i+1, j))) {
                        printf("How? NaN u\n");
                        printf("Old U: %f\nDpu: %f\nRho: %.15g\n", oldu, Dpu, rho.u(i+1, j));
                        printf("Value: %f\n", dt * Dpu);
                        std::exit(0);
                    }
                }
                else {
                    vel.u(i+1, j) = vel.u(i+1, j) - dt * Dpu;
                }
                if(isnan(vel.u(i+1, j))) {
                    printf("NaN u\n");
                    printf("Old U: %f\nDpu: %f\nRho: %.15g\n", oldu, Dpu, rho.u(i+1, j));
                    printf("Value: %f\n", dt * Dpu);
                    std::exit(0);
                }
                vel.u(i+1, j) *= 0.995;
                //vel.u(i+1, j) = vel.u(i+1, j) - dt * Dpu;
            }
            if(j+1 < n) {
                double Dpv = (pressure.get(i, j+1) - pressure.get(i, j)) / dx;
                double oldv = vel.v(i, j+1);
                if(rho.v(i, j+1) > 1e-6) {
                    vel.v(i, j+1) = vel.v(i, j+1) - dt * Dpv / rho.v(i, j+1);
                    //printf("Rho: %.15g\n", rho.v(i, j+1));
                    if(isnan(vel.v(i, j+1))) {
                        printf("How? NaN v\n");
                        printf("Old V: %f\nDpv: %f\nRho: %.15g\n", oldv, Dpv, rho.v(i, j+1));
                        printf("Value: %f\n", dt * Dpv);
                        std::exit(0);
                    }
                }
                else {
                    vel.v(i, j+1) = vel.v(i, j+1) - dt * Dpv;
                }
                if(isnan(vel.v(i, j+1))) {
                    printf("NaN v\n");
                    printf("Old V: %f\nDpv: %f\nRho: %.15g\n", oldv, Dpv, rho.v(i, j+1));
                    printf("Value: %f\n", dt * Dpv);
                    std::exit(0);
                }
                vel.v(i, j+1) *= 0.995;
            }
        }
	}
}

void Fluid::gridToParticles() {
    // TODO: update particle velocities using vel and velOld grids
    for(int i = 0; i < particles.size(); i++) {
		Particle* p = particles[i];
		Vector2d pic = vel.interpolate(p->x);
		Vector2d flip = p->v + (vel.interpolate(p->x) - velOld.interpolate(p->x));
		p->v = (0.95 * flip) + (0.05 * pic);
	}
}

void Fluid::addDensity() {
    //Simpler method for accumulating the densities that gave better results
    rho.assign(Vector2d(0,0));
    rhoWeights.assign(Vector2d(1e-3,1e-3)); // avoid division by zero
    // TODO: for each particle, add its velocity to its nearby grid cells.
    // Use velWeight to determine the weighted average denominator
    for(int i = 0; i < particles.size(); i++) {
		Particle* p = particles[i];
		rho.addInterpolated(p->x, Vector2d(p->m, p->m));
		rhoWeights.addInterpolated(p->x, Vector2d(1, 1));
	}
	for(int i = 0; i < rho.u.m; i++) {
		for(int j = 0; j < rho.u.n; j++) {
			rho.u(i, j) /= rhoWeights.u.get(i, j);
		}
	}
	for(int i = 0; i < rho.v.m; i++) {
		for(int j = 0; j < rho.v.n; j++) {
			rho.v(i, j) /= rhoWeights.v.get(i, j);
		}
	}
    /*
    // Accumulate the particle masses on to the density of the cell face
    // density = Sum for j in Neighbors of mj*W(|r-rj|, h)
    //              mj - mass of particle j
    //              r  - position of the center of the cell face
    //              rj - position of particle j
    //              h  - radius of the kernel (max distance of the neighbors we want to consider)
    double h = 3.0*dx;
    //std::ofstream out("Matrix.txt");
    // 1. Loop over all cells
    //out << "Rho u:\n";
    for(int i = 0; i < rho.u.m; i++) {
		for(int j = 0; j < rho.u.n; j++) {
            // 2. Find the center of the cell
            Vector2d center = rho.u.x0 + Vector2d(i*rho.u.dx, j*rho.u.dx);
            //printf("Center u: (%f, %f)\n", center(0), center(1));
            double density = 0.0;
            // 3. We don't have any acceleration structure so we'll just need to loop over all particles
            for(int p = 0; p < particles.size(); p++) {
                Particle *par = particles[p];
                // 4. Get distance of particle to center
                double dist = (center - par->x).norm();
                // 5. Add weighted contribution mj*W(|r-rj|, h)
                density += par->m * weight(dist, h);
            }
            // 6. Put the total in that cell face
            rho.u(i, j) = density;
            //out << density << ", ";
		}
        //out << "\n";
	}
    //out << "\nRho v:\n";
	for(int i = 0; i < rho.v.m; i++) {
		for(int j = 0; j < rho.v.n; j++) {
            // 2. Find the center of the cell
            Vector2d center = rho.v.x0 + Vector2d(i*rho.v.dx, j*rho.v.dx-rho.v.dx/2.0);
            //printf("Center v: (%f, %f)\n", center(0), center(1));
            double density = 0.0;
            // 3. We don't have any acceleration structure so we'll just need to loop over all particles
            for(int p = 0; p < particles.size(); p++) {
                Particle *par = particles[p];
                // 4. Get distance of particle to center
                double dist = (center - par->x).norm();
                // 5. Add weighted contribution mj*W(|r-rj|, h)
                density += par->m * weight(dist, h);
            }
            // 6. Put the total in that cell face
            rho.v(i, j) = density;
            //out << density << ", ";
		}
        //out << "\n";
	}
    //out.close();
    //std::exit(0);
    */
}

double Fluid::weight(double r, double h) {
    if(r >= 0 && r <= h) {
        return (315.0 / (64.0*M_PI*std::pow(h, 9))) * (std::pow(h*h - r*r, 3));
        //return (1.0 - r/h) * (1.0 - r/h);
    }
    return 0;
}

void MouseForce::addForces(Fluid *fluid, double dt) {
    if (!enabled)
        return;
    for (int i = 0; i <= fluid->m; i++) {
        for (int j = 0; j < fluid->n; j++) {
            Vector2d x = fluid->x0 + Vector2d(i, j+0.5)*fluid->dx;
            if ((this->x - x).norm() <= r)
                fluid->vel.u(i,j) = v(0);
        }
    }
    for (int i = 0; i < fluid->m; i++) {
        for (int j = 0; j <= fluid->n; j++) {
            Vector2d x = fluid->x0 + Vector2d(i+0.5, j)*fluid->dx;
            if ((this->x - x).norm() <= r)
                fluid->vel.v(i,j) = v(1);
        }
    }
}

void GravityForce::addForces(Fluid *fluid, double dt) {
    if(!enabled) {
        return;
    }
    for(int i = 0; i <= fluid->m; i++) {
        for(int j = 0; j <= fluid->n; j++) {
            if(j != fluid->n) {
                fluid->vel.u(i, j) += g(0) * dt;
            }
            if(i != fluid->m) {
                fluid->vel.v(i, j) += g(1) * dt;
            }
        }
    }
}
