#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define PARTICLE_NUMBER 100
#define BOX_WIDTH 20.0
#define BOX_HEIGHT 20.0
#define DELTA_T 0.05
#define KB_T 1.0  // increasing temperature will increase velocity
#define PI 3.14159265358979323846
#define FIILENAMELENGTH 100
#define DIM 2
/*Read: In this first iteration we are not implementing any particle interactions. Data processing is done in python.
*
*Note that the following parameters for the simulation we have not configured yet: temperatureinitial,  size of the time interval used in the integration algorithm, deltat,
* is the time at which the simulation ends, having started at 0.0,
*is an integer giving the number of integration steps between data writes,
*/

// Structure to hold particle data
typedef struct {
    double x, y;  // Position
    double v_x, v_y;  // Velocity
    
} Particle;

Particle particles[PARTICLE_NUMBER];



void adjustCenterOfMassVelocity() {

    double totalVx = 0.0, totalVy = 0.0;

    for (int i = 0; i < PARTICLE_NUMBER; i++) {
        
        totalVx += particles[i].v_x;
        totalVy += particles[i].v_y;
    }

    double centerMassVx = totalVx / PARTICLE_NUMBER;
    double centerMassVy = totalVy / PARTICLE_NUMBER;

    // Subtracting the center of mass velocity from each particle's velocity
    for (int i = 0; i < PARTICLE_NUMBER; i++) {
        particles[i].v_x -= centerMassVx;
        particles[i].v_y -= centerMassVy;
    }

}

//method samples from gaussian distribution
double samplegaussian(){
    
    double u1, u2, z0;
    
    u1 = drand48();
    u2 = drand48();
    
    z0 = sqrt(-2.0 * log(u1)) * cos(2*PI * u2);
    
    return (sqrt(KB_T) * z0);
    
}


void initializeVelocities() {
    for (int i = 0; i < PARTICLE_NUMBER; i++) {

        for (int d=0; d<DIM; d++) {
            
            particles[i].v_x = samplegaussian();
            particles[i].v_y = samplegaussian();
        }
    }

}


// Function to initialize particles
void initializeParticles() {
    int i;
    for (i = 0; i < PARTICLE_NUMBER; i++) {
        particles[i].x = drand48() * BOX_WIDTH;
        particles[i].y = drand48() * BOX_HEIGHT;
        
       
        particles[i].v_x = 0;  
        particles[i].v_y = 0;  
    }

    initializeVelocities();
    adjustCenterOfMassVelocity();

}


// Function to update particle positions
void updatePositions() {    
    int i;
    for (i = 0; i < PARTICLE_NUMBER; i++) {

        particles[i].x += particles[i].v_x * DELTA_T;
        particles[i].y += particles[i].v_y * DELTA_T;

        // Handle boundary conditions
        if (particles[i].x > BOX_WIDTH) particles[i].x -= BOX_WIDTH;
        if (particles[i].y > BOX_HEIGHT) particles[i].y -= BOX_HEIGHT;
        if (particles[i].x < 0) particles[i].x += BOX_WIDTH;
        if (particles[i].y < 0) particles[i].y += BOX_HEIGHT;
    }
}

//method that handles writing data to csv file
void writepositions(char filename[FIILENAMELENGTH]){
    
    int i;
    char filenameupdated[FIILENAMELENGTH];
    
    FILE *posfile;
    
    strcpy(filenameupdated, filename);
    strcat(filenameupdated, "_particle_positions.csv");
    
    posfile = fopen(filenameupdated, "a");

    
    
    for (i=0; i<(PARTICLE_NUMBER-1); i++) {
        fprintf(posfile, ",%f,%f", particles[i].x, particles[i].y);
    }
    
    fprintf(posfile, ",%f,%f", particles[PARTICLE_NUMBER-1].x, particles[PARTICLE_NUMBER-1].y);
    fprintf(posfile,"\n");
   
    fclose(posfile);
    
}







// MAIN
int main(int argc, const char * argv[]) {

    char outputfilename[FIILENAMELENGTH]; // array for file name(s) to which output is saved.

    strcpy(outputfilename, argv[0]);


    srand((long)time(NULL));  // Seed for random number generation
    initializeParticles();

    int steps = 1000;  // Number of time steps to simulate !!This can be changed to a while 
    int i;
    for (i = 0; i < steps; i++) {

        updatePositions();
        writepositions (outputfilename);
       // writeenergies (kineticenergy, currenttime, outputfilename);

    }       


    return 0;
}
