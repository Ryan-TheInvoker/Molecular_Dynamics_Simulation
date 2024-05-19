#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define PARTICLE_NUMBER 550
#define BOX_WIDTH 40.0
#define BOX_HEIGHT 40.0
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

double virial = 0.0;



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


    double trunc_prime = -24.0 * EPSILON * (2.0 * pow(SIGMA/ rc, 12) - pow(SIGMA / rc, 6)) / (rc * rc); 
    potential_E += 4* EPSILON *(sr12 - sr6 - 4* EPSILON * (pow((SIGMA/rc),12) - pow((SIGMA/rc),6)) - r * trunc_prime);

    //printf("trunc prime:%f");
    //double forceMagnitude =  24 * EPSILON * (2 * sr12 / r - sr6 / r); //TODO: check this expression
    //forceMagnitude += 24 * EPSILON * (2 * pow((SIGMA/rc),12) / rc -  pow((SIGMA/rc),6) / rc);   //trunctation expression

    double forceMagnitude = +24 * EPSILON* sr2*(2.0 * sr12 + trunc_prime * sqrt(1/(r*r)));
    //TODO: I think this is probably incredibly inefficent, calculating once is probably only necessary

    //virial -= forceMagnitude* r*r;
    return forceMagnitude;
}

double lennard_jones_potential(double r, double rc){

    

    double sr = SIGMA / r;
    double sr2 = sr * sr;
    double sr6 = sr2 * sr2 * sr2;
    double sr12 = sr6 * sr6;


    double trunc_prime = -24.0 * EPSILON * (2.0 * pow(SIGMA/ rc, 12) - pow(SIGMA / rc, 6)) / (rc * rc); 
    potential_E += 4* EPSILON *(sr12 - sr6 - 4* EPSILON * (pow((SIGMA/rc),12) - pow((SIGMA/rc),6)) - r * trunc_prime);

    //printf("trunc prime:%f");
    //double forceMagnitude =  24 * EPSILON * (2 * sr12 / r - sr6 / r); //TODO: check this expression
    //forceMagnitude += 24 * EPSILON * (2 * pow((SIGMA/rc),12) / rc -  pow((SIGMA/rc),6) / rc);   //trunctation expression



}


void apply_boundary_conditions(int i) {
    if (particles[i].x >= BOX_WIDTH) particles[i].x -= BOX_WIDTH;
    else if (particles[i].x < 0) particles[i].x += BOX_WIDTH;

    if (particles[i].y >= BOX_HEIGHT) particles[i].y -= BOX_HEIGHT;
    else if (particles[i].y < 0) particles[i].y += BOX_HEIGHT;
}



void calculate_new_forces_and_accelerations(int i) {
    virial =0.0;
    potential_E = 0.0;
    fx[i] = 0.0;
    fy[i] = 0.0;

    for (int j = 0; j < PARTICLE_NUMBER; j++) {
        if (i != j) {

            double dx = particles[i].x - particles[j].x;
            double dy = particles[i].y - particles[j].y;

            dx = minimumImageDistance(dx, BOX_WIDTH);
            dy = minimumImageDistance(dy, BOX_HEIGHT);

            double r = sqrt(dx * dx + dy * dy);
            if(r<rc){
                double force_mag = lennard_jones(r, 2.5 * SIGMA);
                fx[i] += force_mag * dx; 
                fy[i] += force_mag * dy; // think we must remove /r here

                virial += fx[i]*dx + fy[i]*dy;
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
    
    //printf("T: %f \n",T);
    //beta for temperature correction   
    double beta = sqrt(KB_T/T);

    


    // // 1/2 v step
    // particles[i].v_x += 0.5 * particles[i].a_x * DELTA_T;
    // particles[i].v_y += 0.5 * particles[i].a_y * DELTA_T;

    particles[i].x += particles[i].v_x * DELTA_T+ 0.5 * particles[i].a_x * DELTA_T*DELTA_T;
    particles[i].y +=  particles[i].v_y * DELTA_T+ 0.5 * particles[i].a_y * DELTA_T*DELTA_T;
    
    particles[i].v_x += 0.5*(particles[i].a_x )*DELTA_T;
    particles[i].v_y += 0.5*(particles[i].a_y )*DELTA_T;
    // // Update positions
    // particles[i].x += particles[i].v_x * DELTA_T + 0.5 * particles[i].a_x * DELTA_T * DELTA_T;
    // particles[i].y += particles[i].v_y * DELTA_T + 0.5 * particles[i].a_y * DELTA_T * DELTA_T;

    
    calculate_new_forces_and_accelerations(i);

    particles[i].v_x += 0.5*(particles[i].a_x )*DELTA_T;
    particles[i].v_y += 0.5*(particles[i].a_y )*DELTA_T;

    apply_boundary_conditions(i);

    particles[i].v_x *= beta;
    particles[i].v_y *= beta;
    kinetic_E += 0.5*MASS* (particles[i].v_x * particles[i].v_x  + particles[i].v_y * particles[i].v_y);
    // // Update velocities again with new accelerations
    //  particles[i].v_x += 0.5 * particles[i].a_y * DELTA_T;
    //  particles[i].v_y += 0.5 * particles[i].a_y  * DELTA_T;


    

    
    //adjustCenterOfMassVelocity();
    


}




void update_step(){

    
    for (int i = 0; i < PARTICLE_NUMBER; i++)
    {
        verlet(i);
        //printf("Fy: %f \t Fx: %f",fx[i],fy[i]);
    }
    
    
}





//SECTION 3 METHODS//
void delete_old_files(const char *outputfilename) {
    char filename[FIILENAMELENGTH];

    snprintf(filename, sizeof(filename), "%s_KE_list.txt", outputfilename);
    remove(filename);

    snprintf(filename, sizeof(filename), "%s_PE_list.txt", outputfilename);
    remove(filename);

    snprintf(filename, sizeof(filename), "%s_particle_positions.csv", outputfilename);
    remove(filename);

    snprintf(filename, sizeof(filename), "%s_TP_data.txt", outputfilename);
    remove(filename);

    snprintf(filename, sizeof(filename), "%s_Potential_Energy_Plot", outputfilename);
    remove(filename);

    snprintf(filename, sizeof(filename), "%s_Kinetic_Energy_Plot.png", outputfilename);
    remove(filename);

    snprintf(filename, sizeof(filename), "%s_Total_Energy_Plot.png", outputfilename);
    remove(filename);

    snprintf(filename, sizeof(filename), "%s_TP_Plot.png", outputfilename);
    remove(filename);

    snprintf(filename, sizeof(filename), "particle_motion_2_ass2.mp4");
    remove(filename);
}

void print_timing_summary() {
    
    // printf("initializeParticles: %f seconds\n", (double)total_time_initializeParticles / CLOCKS_PER_SEC);
    // printf("initializeVelocities: %f seconds\n", (double)total_time_initializeVelocities / CLOCKS_PER_SEC);
    // printf("adjustCenterOfMassVelocity: %f seconds\n", (double)total_time_adjustCenterOfMassVelocity / CLOCKS_PER_SEC);
   // printf("build_neighbor_list: %f seconds\n", (double)total_time_build_neighbor_list / CLOCKS_PER_SEC);
    // printf("calculate_new_forces_and_accelerations: %f seconds\n", (double)total_time_calculate_new_forces_and_accelerations / CLOCKS_PER_SEC);
    // printf("verlet: %f seconds\n", (double)total_time_verlet / CLOCKS_PER_SEC);
    // //printf("update_step: %f seconds\n", (double)total_time_update_step / CLOCKS_PER_SEC);
    //  printf("temp_calc: %f seconds\n", (double)total_calc_temp / CLOCKS_PER_SEC);
    //   printf("len_jon: %f seconds\n", (double)total_len_jon / CLOCKS_PER_SEC);
    //    printf("len_jon_pot: %f seconds\n", (double)total_len_jon_pot / CLOCKS_PER_SEC);
    //     printf("MD: %f seconds\n", (double)total_MD / CLOCKS_PER_SEC);
    //      printf("WKE: %f seconds\n", (double)total_WKE / CLOCKS_PER_SEC);
    //       printf("WPE: %f seconds\n", (double)total_WPE / CLOCKS_PER_SEC);
}




/*
*Solely for T v P data
*/
void record_TP_data_point(double temperature, double pressure, char filename[FIILENAMELENGTH]) {
    
    char filenameupdated[FIILENAMELENGTH];
    FILE *datafile;
    strcpy(filenameupdated, filename);
    strcat(filenameupdated, "_TP_data.txt");
    datafile = fopen(filenameupdated, "a");
    fprintf(datafile, "%f %f\n", temperature, pressure);
    fclose(datafile);
}


/*
*
*/
double calculate_temperature() {
    clock_t start = clock();
    double total_kinetic_energy = 0.0;
    for (int i = 0; i < PARTICLE_NUMBER; i++) {
        double speed_squared = particles[i].v_x * particles[i].v_x + particles[i].v_y * particles[i].v_y;
        total_kinetic_energy += 0.5 * MASS * speed_squared;
    }
    return (2.0 * total_kinetic_energy) / (DIM * PARTICLE_NUMBER * KB_T);

    //total_calc_temp += clock() - start;
}


/*Pressure calculation
*
*/
double compute_pressure( double V) {

    // double virial_sum = 0.0;

    // // Compute the virial term
    // for (int i = 0; i < PARTICLE_NUMBER; i++) {
    //     for (int j = i + 1; j < PARTICLE_NUMBER; j++) {

    //         double dx = particles[j].x - particles[i].x;
    //         double dy = particles[j].y - particles[i].y;

    //         double r = sqrt(dx * dx + dy * dy);

    //         if (r < rc) {
    //             double fx_ij = (fx[i] - fx[j]);
    //             double fy_ij = (fy[i] - fy[j]);
    //             double Fij_dot_rij = fx_ij * dx + fy_ij * dy;
    //             virial_sum += Fij_dot_rij / r;
    //         }
    //     }
    // }

    // Calculate pressure using the virial equation of state
    printf("Virial = %f", virial);
    double pressure = (PARTICLE_NUMBER / V) * calculate_temperature()  + (virial / (DIM * V)) ;

    return pressure;
}




/*
*Function to adjust the velocities in order to simulate different temperatures
*/
void adjust_velocities( double factor) {
    for (int i = 0; i < PARTICLE_NUMBER; i++) {

            particles[i].v_x += particles[i].v_x * factor;
            particles[i].v_y += particles[i].v_y * factor;
        
    }
}

void equilibrate_system(int steps) {
    for (int k = 0; k < steps; k++) {
        // Perform one MD step (function not provided, placeholder for your MD step implementation)
        update_step();
    }
}



// CALORIC DATA FUNCTIONS//
void write_caloric( char filename[FIILENAMELENGTH]){



    char filenameupdated[FIILENAMELENGTH];
    
    FILE *TEfile;
    
    strcpy(filenameupdated, filename);
    strcat(filenameupdated, "_TE_list.txt"); // data written to filename that is pure text
    
    TEfile = fopen(filenameupdated, "a");

    fprintf(TEfile, "%f %f\n", calculate_temperature(), potential_E+kinetic_E);
    
    fclose(TEfile);


}

void caloric_data(int num_runs,int steps,char filename[FIILENAMELENGTH]){
    
    for (int i = 0; i < num_runs; i++)
    {
        equilibrate_system(500);
        write_caloric(filename);
        adjust_velocities(0.01);
    }
    


}
// MAIN
int main(int argc, const char * argv[]) {

    char outputfilename[FIILENAMELENGTH]; // array for file name(s) to which output is saved.

    
   
    strcpy(outputfilename, argv[0]);
    delete_old_files(outputfilename);

    srand(time(NULL));  // Seed for random number generation
    initializeParticles();

   

    int steps = 4000;  // Number of time steps to simulate !!This can be changed to a while 
    int i;
    for (i = 0; i < steps; i++) {
        
        kinetic_E = 0.0;
        

        update_step();
        
        writepositions (outputfilename);
        write_KE ( i*DELTA_T, outputfilename);
        write_PE( i*DELTA_T,outputfilename);
        
    }

    // int num_readings = 10;
    
    // for (int j = 0; j < num_readings; j++)
    // {
    
    //     equilibrate_system(500);
    //     record_TP_data_point(calculate_temperature(), compute_pressure(BOX_HEIGHT*BOX_WIDTH),outputfilename);

    //     adjust_velocities(0.01);
    // }
    
    //print_timing_summary(); 

    //caloric_data(10,500,outputfilename);

    return 0;
}