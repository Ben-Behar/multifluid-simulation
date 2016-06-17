#include "fluid.hpp"
#include "opengl.hpp"

#include <iostream>
#include <cstdio>
using namespace std;

int w = 800, h = 800;                          // size of window in pixels
double xmin = 0, xmax = 1, ymin = 0, ymax = 1; // range of coordinates drawn
double dt = 0.01;                              // time step
int lastTime = 0, prevTime = 0, frame = 0;

int n = 50;
Fluid fluid;
MouseForce mouseForce(0.1);

bool stopPhysics = false;

// random double uniformly distributed between 0 and 1
/// double orand() {
    /// return (double)rand()/RAND_MAX;
/// }


void stage(int i)
{
    stopPhysics = true;
    fluid.vel.assign(Vector2d(0,0)); // reset vel
    fluid.pressure.assign(0);        // reset pressure
    fluid.particles.clear();
    if(i == 1)
    {   
        for (int i = 0; i < fluid.m/2; i++) {
            for (int j = 0; j < 3*fluid.n/4; j++) {
                Vector2d xij = fluid.x0 + Vector2d(i, j)*fluid.dx;
                for (int i1 = 0; i1 < 4; i1++) {
                    for (int j1 = 0; j1 < 4; j1++) {
                        Vector2d x = xij + Vector2d(i1, j1)*(fluid.dx/2)
                            + Vector2d(urand(), urand())*(fluid.dx/2);
                        fluid.particles.push_back(new Particle(x, Vector2d(0,0), Vector3d(0, 0, 1), 1));
                    }
                }
            }
        }
        for(int i = 3*fluid.m/4; i < fluid.m; i++) {
            for(int j = 0; j < fluid.n/2; j++) {
                Vector2d xij = fluid.x0 + Vector2d(i, j)*fluid.dx;
                for(int i1 = 0; i1 < 4; i1++) {
                    for(int j1 = 0; j1 < 4; j1++) {
                        Vector2d x = xij + Vector2d(i1, j1)*(fluid.dx/2) + Vector2d(urand(), urand())*(fluid.dx/2);
                        fluid.particles.push_back(new Particle(x, Vector2d(0, 0), Vector3d(0, 0, 0), 0.5));
                    }
                }
            }
        }
    }
    else if(i == 2)
    {
        for (int i = 0; i < fluid.m/2; i++) {
            for (int j = 0; j < 3*fluid.n/4; j++) {
                Vector2d xij = fluid.x0 + Vector2d(i, j)*fluid.dx;
                for (int i1 = 0; i1 < 4; i1++) {
                    for (int j1 = 0; j1 < 4; j1++) {
                        Vector2d x = xij + Vector2d(i1, j1)*(fluid.dx/2)
                            + Vector2d(urand(), urand())*(fluid.dx/2);
                        fluid.particles.push_back(new Particle(x, Vector2d(0,0), Vector3d(0, 0, 1), 5));
                    }
                }
            }
        }
        for(int i = 3*fluid.m/4; i < fluid.m; i++) {
            for(int j = 0; j < fluid.n/2; j++) {
                Vector2d xij = fluid.x0 + Vector2d(i, j)*fluid.dx;
                for(int i1 = 0; i1 < 4; i1++) {
                    for(int j1 = 0; j1 < 4; j1++) {
                        Vector2d x = xij + Vector2d(i1, j1)*(fluid.dx/2) + Vector2d(urand(), urand())*(fluid.dx/2);
                        fluid.particles.push_back(new Particle(x, Vector2d(0, 0), Vector3d(0, 0, 0), 0.5));
                    }
                }
            }
        }
    }
    else if (i == 3) {
        fluid.particleSize = 5;
        for(int i = 0; i < fluid.m; i++) {
            for(int j = 0; j < fluid.n; j++) {
                Vector2d xij = fluid.x0 + Vector2d(i, j)*fluid.dx;
                double pmass = (i < fluid.m/2 && j < 3*fluid.n/4) ? 5.0 : 0.05;
                Vector3d pcolor = (i < fluid.m/2 && j < 3*fluid.n/4) ? Vector3d(0, 0, 1) : Vector3d(0.9, 0.9, 0.9);
                for (int i1 = 0; i1 < 4; i1++) {
                    for (int j1 = 0; j1 < 4; j1++) {
                        Vector2d x = xij + Vector2d(i1, j1)*(fluid.dx/2)
                            + Vector2d(urand(), urand())*(fluid.dx/2);
                        fluid.particles.push_back(new Particle(x, Vector2d(0,0), pcolor, pmass));
                    }
                }
            }
        }
    }
    stopPhysics = false;
}

void display() {


    glClearColor(1,1,1,1);
    glClear(GL_COLOR_BUFFER_BIT);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(::xmin, ::xmax, ::ymin, ::ymax);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    fluid.draw();
    glutSwapBuffers();
}

void reshape(int w, int h) {
    ::w = w;
    ::h = h;
    glViewport(0, 0, w, h);
}

void idle() {
    //FPS
    frame++;
    int time = glutGet(GLUT_ELAPSED_TIME);
    if(time - lastTime > 1000) {
        double fps = frame * 1000.0 / (time-lastTime);
        lastTime = time;
        frame = 0;
        //printf("FPS: %f\n", fps);
        char fpsString[15];
        snprintf(fpsString, 15, "%f", fps);
        glutSetWindowTitle(fpsString);
    }
    double timeStep = (time-prevTime) / 1000.0;
    if(!stopPhysics)  // allows us to make modifications to the fluid without causing concurrent issues
    {
        int substeps = 1;
        for (int i = 0; i < substeps; i++) {
            fluid.step(::dt);
            //fluid.step(timeStep/substeps);
        }
    }
    prevTime = time;

    glutPostRedisplay();
}

void mouse(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON) {
        if (state == GLUT_DOWN) {
            ::mouseForce.enabled = true;
            ::mouseForce.x = Vector2d(xmin + (double)x/w * (xmax - xmin),
                                      ymax - (double)y/h * (ymax - ymin));
            ::mouseForce.v = Vector2d(0,0);
        } else if (state == GLUT_UP) {
            ::mouseForce.enabled = false;
        }
    }
}

void motion(int x, int y) {
    Vector2d xOld = ::mouseForce.x;
    ::mouseForce.x = Vector2d(xmin + (double)x/w * (xmax - xmin),
                              ymax - (double)y/h * (ymax - ymin));
    ::mouseForce.v = (::mouseForce.x - xOld)/dt;
}

void onKey(unsigned char key, int x, int y) {
    printf("%i\n",key);
	switch(key) {
		case 27: // escape
			exit(0);
        case 32: // space
            fluid.particles.push_back(
                new Particle(
                    Vector2d((double)rand()/RAND_MAX,.9), // position
                    Vector2d(0, 0),                       // velocity
                    Vector3d(1, 0, 0),                    // color
                    1000));                               // Mass
            break;
        case '=': // +
            fluid.particleSize ++;
            break;
        case '-': // -
            fluid.particleSize --;
            break;
        case 'q': // q
            fluid.showPressure = !fluid.showPressure;
            break;
        case 'w': // w
            fluid.showVelocity = !fluid.showVelocity;
            break;
        case '1': // 1
            stage(1);
            break;
        case '2': // 2
            stage(2);
            break;
        case '3': // 2
            stage(3);
            break;
	}
}

int main(int argc, char **argv) {
    fluid.init(n, n, Vector2d(0, 0), 1./n);
    fluid.forces.push_back(&mouseForce);
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_MULTISAMPLE);
    glutInitWindowSize(::w, ::h);
    glutCreateWindow("Animation");
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutIdleFunc(idle);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutKeyboardFunc(onKey);
    glutMainLoop();
}
