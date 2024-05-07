#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define PARTICLE_NUMBER 5
#define BOX_WIDTH 10.0
#define BOX_HEIGHT 10.0
#define DELTA_T 0.01
#define KB_T 1.0  // increasing temperature will increase velocity
#define PI 3.14159265358979323846
#define FIILENAMELENGTH 100
#define DIM 2
#define MASS 1.0
#define SIGMA 1.0
#define EPSILON 1.0 
#define DELTA 0.01
#define KB 1.38e-23

//Ercolessi "epsilon = 1 and sigma = 1"

//Sigma: dictates where the potential is zero
//Epsilon: the depth of the potential well ==> strength of the interactions


//QUESTIONS 
/* What should sigma and epsilon be?
* Temperature?
* Are we expecting the system to reach an equilibrium? --> Clustering?
*/ 


/*Read: 
*
*Note that the following parameters for the simulation we have not configured yet: temperatureinitial,  size of the time interval used in the integration algorithm, deltat,
* is the time at which the simulation ends, having started at 0.0,
*is an integer giving the number of integration steps between data writes,
*/

/* Phase 2 Simulation Methods: 
*
*/

/*1) want to neglect forces from particles that are far away  (DONE)
* 2) don't calculate pairs twice                       
* 3) calculate potential energy and kinetic energy  (DONE)
*/



// Structure to hold particle data
typedef struct {
    double x, y;  // Position
    double v_x, v_y;  // Velocity
    double a_x,a_y;
} Particle;


//GLOBAL VARIABLES
Particle particles[PARTICLE_NUMBER];

//better convention to have variables declared in main method??
double fx[PARTICLE_NUMBER];
double fy[PARTICLE_NUMBER];
double rc = 2.5*SIGMA;
double kinetic_E = 0.0;
double potential_E = 0.0;




/*
*This function Ensures systemic drift is removed
*/
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


//Initializing particles on a lattice
void initializeParticles(){
 
    int nrows, ncols;
    double minDiff = HUGE_VAL;
    
    // Find the best nrows and ncols
    for (int i = 1; i <= sqrt(PARTICLE_NUMBER); i++) {
        if (PARTICLE_NUMBER % i == 0) {
            int j = PARTICLE_NUMBER / i;
            double diff = fabs(i - j);
            if (diff < minDiff) {
                minDiff = diff;
                nrows = i;
                ncols = j;
            }
        }
    }

    // Calculate spacing
    double spacingX = BOX_WIDTH / ncols;
    double spacingY = BOX_HEIGHT / nrows;

    // Initialize particle positions
    int particleIndex = 0;
    for (int i = 0; i < nrows; i++) {
        for (int j = 0; j < ncols; j++) {
            if (particleIndex < PARTICLE_NUMBER) {
                particles[particleIndex].x = j * spacingX + spacingX / 2;
                particles[particleIndex].y = i * spacingY + spacingY / 2;
               
                particleIndex++;
            }
        }
    }

    initializeVelocities();
    adjustCenterOfMassVelocity();
}


/*Writes total kinetic energy to a specified file.
*
*/
void write_KE( double current_time, char filename[FIILENAMELENGTH]){

    char filenameupdated[FIILENAMELENGTH];
    
    FILE *KEfile;
    
    strcpy(filenameupdated, filename);
    strcat(filenameupdated, "_KE_list.txt"); // data written to filename that is pure text
    
    KEfile = fopen(filenameupdated, "a");

    fprintf(KEfile, "%f %f\n", current_time, kinetic_E);
    
    fclose(KEfile);

}


void write_PE( double current_time, char filename[FIILENAMELENGTH]){

    char filenameupdated[FIILENAMELENGTH];
    
    FILE *PEfile;
    
    strcpy(filenameupdated, filename);
    strcat(filenameupdated, "_PE_list.txt"); // data written to filename that is pure text
    
    PEfile = fopen(filenameupdated, "a");

    fprintf(PEfile, "%f %f\n", current_time, potential_E);
    
    fclose(PEfile);

}


void writepositions(char filename[FIILENAMELENGTH]){
    
    int i;
    char filenameupdated[FIILENAMELENGTH];
    
    FILE *posfile;
    
    strcpy(filenameupdated, filename);
    strcat(filenameupdated, "_particle_positions.csv"); // data written to filename that is formated and designated a mathematica notebook
    
    posfile = fopen(filenameupdated, "a");

    
    
    for (i=0; i<(PARTICLE_NUMBER-1); i++) {
        fprintf(posfile, ",%f,%f", particles[i].x, particles[i].y);
    }
    
    fprintf(posfile, ",%f,%f", particles[PARTICLE_NUMBER-1].x, particles[PARTICLE_NUMBER-1].y);
    fprintf(posfile,"\n");
    
    fclose(posfile);
    
}




//implementing the periodic images via the minimum image distance idea.
double minimumImageDistance(double dx, double boxDimension) {
    dx -= round(dx / boxDimension) * boxDimension;
    return dx;
}



double lennard_jones(double r, double rc){

    double sr = SIGMA / r;
    double sr2 = sr * sr;
    double sr6 = sr2 * sr2 * sr2;
    double sr12 = sr6 * sr6;


    potential_E += 4* EPSILON *(sr12 - sr6) - 4* EPSILON * (pow((SIGMA/rc),12) - pow((SIGMA/rc),6));

    double forceMagnitude = - 24 * EPSILON * (2 * sr12 / r - sr6 / r); //TODO: check this expression
    forceMagnitude += 24 * EPSILON * (2 * pow((SIGMA/rc),12) / rc -  pow((SIGMA/rc),6) / rc);   //trunctation expression

    //TODO: I think this is probably incredibly inefficent, calculating once is probably only necessary

    return forceMagnitude;
}



void apply_boundary_conditions(int i) {
    if (particles[i].x >= BOX_WIDTH) particles[i].x -= BOX_WIDTH;
    else if (particles[i].x < 0) particles[i].x += BOX_WIDTH;

    if (particles[i].y >= BOX_HEIGHT) particles[i].y -= BOX_HEIGHT;
    else if (particles[i].y < 0) particles[i].y += BOX_HEIGHT;
}



void calculate_new_forces_and_accelerations(int i) {

    fx[i] = 0.0;
    fy[i] = 0.0;
    for (int j = 0; j < PARTICLE_NUMBER; j++) {
        if (i != j) {

            double dx = particles[j].x - particles[i].x;
            double dy = particles[j].y - particles[i].y;

            dx = minimumImageDistance(dx, BOX_WIDTH);
            dy = minimumImageDistance(dy, BOX_HEIGHT);

            double r = sqrt(dx * dx + dy * dy);
            if(r<rc){
                double force_mag = lennard_jones(r, 2.5 * SIGMA) / MASS;
                fx[i] += force_mag * dx/r; 
                fy[i] += force_mag * dy/r;
            }
        }
    }
    particles[i].a_x = fx[i] / MASS;
    particles[i].a_y = fy[i] / MASS;
}


//Implementation of the velocity-verlet algorithm
void verlet(int i){

    //temperature by equipartition function
    double T = 0.0;
    for (int j = 0; j < PARTICLE_NUMBER; j++)
    {
        T += MASS * ((particles[j].v_x)*(particles[j].v_x) + (particles[j].v_y)*(particles[j].v_y));
    }
    
    T /= (PARTICLE_NUMBER *DIM); //KB removed
    printf("T: %f \n",T);
    //beta for temperature correction   
    double beta = sqrt(KB_T/T);

    


    // 1/2 v step
    particles[i].v_x += 0.5 * particles[i].a_x * DELTA_T;
    particles[i].v_y += 0.5 * particles[i].a_y * DELTA_T;
    //particles[i].v_x = particles[i].v_x * beta + 0.5 * particles[i].a_x * DELTA_T;
    //particles[i].v_y = particles[i].v_y * beta  + 0.5 * particles[i].a_y * DELTA_T;
    //particles[i].v_x = beta * (particles[i].v_x + 0.5 * particles[i].a_x * DELTA_T);
    //particles[i].v_y = beta * (particles[i].v_y + 0.5 * particles[i].a_y * DELTA_T);
    // Update positions
    particles[i].x += particles[i].v_x * DELTA_T + 0.5 * particles[i].a_x * DELTA_T * DELTA_T;
    particles[i].y += particles[i].v_y * DELTA_T + 0.5 * particles[i].a_y * DELTA_T * DELTA_T;

    apply_boundary_conditions(i);
    calculate_new_forces_and_accelerations(i);
   
    // Update velocities again with new accelerations
     particles[i].v_x += 0.5 * particles[i].a_y * DELTA_T;
     particles[i].v_y += 0.5 * particles[i].a_y  * DELTA_T;
    //particles[i].v_x = beta*particles[i].v_x + 0.5 * particles[i].a_y * DELTA_T;
    //particles[i].v_y = beta*particles[i].v_y + 0.5 * particles[i].a_y  * DELTA_T;
     //particles[i].v_x = beta * (particles[i].v_x + 0.5 * particles[i].a_x * DELTA_T);
     //particles[i].v_y = beta * (particles[i].v_y + 0.5 * particles[i].a_y * DELTA_T);

    particles[i].v_x *= beta;
    particles[i].v_y *= beta;
    
    //adjustCenterOfMassVelocity();
    kinetic_E += 0.5*MASS* (particles[i].v_x * particles[i].v_x  + particles[i].v_y * particles[i].v_y);


}




void update_step(){

    
    for (int i = 0; i < PARTICLE_NUMBER; i++)
    {
        verlet(i);
        //printf("Fy: %f \t Fx: %f",fx[i],fy[i]);
    }
    
    
}








// MAIN
int main(int argc, const char * argv[]) {

    char outputfilename[FIILENAMELENGTH]; // array for file name(s) to which output is saved.

    

    strcpy(outputfilename, argv[0]);


    srand(time(NULL));  // Seed for random number generation
    initializeParticles();

   

    int steps = 3000;  // Number of time steps to simulate !!This can be changed to a while 
    int i;
    for (i = 0; i < steps; i++) {
        
        kinetic_E = 0.0;
        potential_E = 0.0;

        update_step();
        
        writepositions (outputfilename);
        write_KE ( i*DELTA_T, outputfilename);
        write_PE( i*DELTA_T,outputfilename);
        
    }

    return 0;
}