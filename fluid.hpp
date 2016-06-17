#ifndef FLUID_HPP
#define FLUID_HPP

#include "grid.hpp"

class Particle;
class Force;

inline double urand() {
    return (double)rand()/RAND_MAX;
}

class Fluid {
public:
    static const char FLUID = 'f', AIR = 'a', OBS = 'o'; // flags
    int m, n;
    Vector2d x0;
    double dx;
    //double rho;                  // density

    int particleSize;

    // debug info flags
    bool showVelocity;
    bool showPressure;

    StaggeredGrid vel, velOld;   // velocity field(s)
    StaggeredGrid rho;           // density between the cell faces
    StaggeredGrid rhoWeights;
    vector<Particle*> particles; // particles for advection
    vector<Force*> forces;       // external forces
    Grid<char> flags;            // whether cell is fluid, air, or obstacle
    // intermediate variables
    StaggeredGrid velWeight;     // used in particlesToGrid()
    Grid<int> indices;           // used in pressureSolve()
    Grid<double> pressure;       // used in pressureSolve()
    void init(int m, int n, Vector2d x0, double dx);
    void draw();
    void step(double dt);
    void advection(double dt);
    void particlesToGrid();
    void flagCells();
    void addExternalForces(double dt);
    void applyBoundaryConditions();
    void pressureSolve(double dt);
    void gridToParticles();
    void addDensity();
    double weight(double r, double h);
};

class Particle {
public:
    Vector2d x, v;
    Vector3d color;
    double m;
    Particle(Vector2d x, Vector2d v, Vector3d color, double m): x(x), v(v), color(color), m(m) {}
};

class Force {
public:
    virtual void addForces(Fluid *fluid, double dt) = 0;
};

class MouseForce: public Force {
public:
    Vector2d x, v;
    double r;
    bool enabled;
    MouseForce(double r): x(0,0), v(0,0), r(r), enabled(false) {}
    void addForces(Fluid *fluid, double dt);
};

class GravityForce: public Force {
public: 
    Vector2d g;
    bool enabled;
    GravityForce(Vector2d g): g(g), enabled(true) {}
    void addForces(Fluid *fluid, double dt);
};

#endif
