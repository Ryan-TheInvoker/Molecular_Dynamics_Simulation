//
//  main.c
//  mol-dyn-2D
//
//
//

/* After compiling the program it is called with 7 command-line arguments.
 In this example code the program will crash if they are not precisely correctly entered here,
 as no attempt is made at error checking.
 
 you call the code in the compiled code as follows
 
 "./a.out 0.0001 0.5 100 1.1 1.2 1.0 datafile"
 
 Here's an example I shall show you: "./a.out 0.01 1.0 1 1.0 1.0 1.0 test1"
 
 the first argument is the size of the time interval used in the integration algorithm, deltat,
 the second argument, after a space, is the time at which the simulation ends, having started at 0.0,
 the third argument is an integer giving the number of integration steps between data writes,
 the box dimensions x and y are the fourht and fifth arguments, respectively,
 the sixth argument is the initial temperature,
 and the final argument is a string which names the file(s) which are created with data.
 
 
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define PARTICLES 10  // define the number of particles in the simulation
#define DIM 2 // We restrict ourselves to two dimensions

#define TWOPI 6.283185307 // 2*pi

#define FIILENAMELENGTH 100 // Maximum characters for file name length string


/* ------------------------------------------------------------  */

/* This is an implementation of box-Muller approach to sampling from a normal distribution.
 The desired temperature for the Maxwell distribution is passed to the function which returns an
 appropriately normal-disttibruted random variable. This can be positive or negative.
 
 Read up about the box-Muller approach on this. And look at your third-year text-book
 for a reminder of why we sample velocities from the Maxwell distribution.
 
 Eventually, the precise way in which particle velocities are selected is not so important
 for the simulation. The temperature of a gas with interacting particles will generally
 not equal the value that we use here. Nevertheless an approach such as this can more-or-less
 start off with something in a meaningful way.
 
 Later you shall also learn a few things you could do to adjust the temperature.
 
 */

double samplegaussian(double temp){
    
    double u1, u2, z0;
    
    u1 = drand48();
    u2 = drand48();
    
    z0 = sqrt(-2.0 * log(u1)) * cos(TWOPI * u2);
    
    return (sqrt(temp) * z0);
    
}


/* This function initialises the position and velocities in of the particles in the box.
 A function used in this way (and called in the way it is called, changes the entries in the array
 in the function that calls it.  */

void initialiseparticles (double pos[PARTICLES][DIM], double vel[PARTICLES][DIM], double temperature, double boxdims[DIM]){
    
    int i, d;
    
    srand48((long)time(NULL)); // seeding peudorandom number, using inbuilt function - you're wlecome to use other random number generators
    
    for (i=0; i<PARTICLES; i++) {
        
        for (d=0; d<DIM; d++) {
            
            pos[i][d] = drand48()*boxdims[d];
            vel[i][d] = samplegaussian(temperature);
        }
    }
}

/* ------------------------------------------------------------  */

/* This function should eventually use the positions to calculate the total
 force on each particle, which gives the acceleration. The acceleration matrix
 is updated.
 
 This initial code sets the accelerations to 0.0. This is where you would need to
 put in a version of the Lennard-Jones force that calculates the forces between
 pairs of particles. */

void calculateaccelerations(double pos[PARTICLES][DIM], double acc[PARTICLES][DIM]){
    
    int i, j, d; // i and j would enerumerate particle pairs and d is the Cartesian component index.
    
    for (i=0; i<PARTICLES; i++) {
        for (d=0; d<DIM; d++) {
            acc[i][d] = 0.0;
        }
    }
    
}


/* ------------------------------------------------------------  */

/* This function performs the so-called integration step, which calculates new
 positions and new velocities for all particles, given the accelerations.
 
 Simultaneously is also computes the total kinetic energy. One can calculate this separately at the
 computational cost of looping through all particles and coordinates another time.
 
 This initial code sets the accelerations to 0.0. This is where you would implement the
 Verlet, velocity-Velet, leapfrog, or other integration algorithms. */

void moveparticles(double pos[PARTICLES][DIM], double vel[PARTICLES][DIM], double acc[PARTICLES][DIM], double dt, double *t, double *KE, double boxdims[DIM]){
    
    int i, d;
    double newcoord, newvcoord, unused;
    
    *KE = 0.0;
    
    for (i=0; i<PARTICLES; i++) {
        for (d=0; d<DIM; d++) {
            
            newcoord = pos[i][d]+ dt*vel[i][d] + 0.5 * dt * dt * acc[i][d]; // This is what you need to replace!
            newvcoord = vel[i][d];
            
            /* Our box is periodic, a particle leaving on one side is placed back into the box at the other.
             Does the following do this?   The function modf gives the part of the first argument after the decimal. */
            
            newcoord = modf(newcoord/boxdims[d]+1.0, &unused) * boxdims[d];
            
            /* Now place the new value into the matrices */
            pos[i][d] = newcoord;
            vel[i][d] = newvcoord;
            
            /* and update the kinetic energy */
            *KE += 0.5*newvcoord*newvcoord;
            
        }
    }
    *t += dt;
}


/* This functions appends the coordinates of the particles to a specific file.
 Here the data will be written in a way that you might use simply in mathematica,
 but there are definitely better ways to read position data (in easier format) if
 you wish to animate their motions.
 
 Please note, however, that particle positions are not really useful data.
 
 In long and large simulations one could fill whole hard disks with such data.
 
 We shall mostly need only energies, pressure, and such, with, quite seldomly,
 a snapshot of positions to illustrate something about the system.
 
 */

void writepositions(double pos[PARTICLES][DIM], char filename[FIILENAMELENGTH], double boxdims[DIM]){
    
    int i;
    char filenameupdated[FIILENAMELENGTH];
    
    FILE *posfile;
    
    strcpy(filenameupdated, filename);
    strcat(filenameupdated, "_mathematica.nb"); // data written to filename that is formated and designated a mathematica notebook
    
    posfile = fopen(filenameupdated, "a");

    fprintf(posfile, "ListPlot[{");
    
    for (i=0; i<(PARTICLES-1); i++) {
        fprintf(posfile, "{ %f, %f },\n", pos[i][0], pos[i][1]);
    }
    
    fprintf(posfile, "{ %f, %f }\n", pos[PARTICLES-1][0], pos[PARTICLES-1][1]);
    
    fprintf(posfile, "}, PlotStyle -> PointSize[Large], \n PlotRange -> {{0, %f}, {0, %f}}, AspectRatio -> 1],\n", boxdims[0], boxdims[1]);
    fclose(posfile);
    
}

/* The energies are written to a simple text file: timestamp and kinetic energy. */

void writeenergies(double KE, double t, char filename[FIILENAMELENGTH]){
    
    char filenameupdated[FIILENAMELENGTH];
    
    FILE *KEfile;
    
    strcpy(filenameupdated, filename);
    strcat(filenameupdated, "_KE_list.txt"); // data written to filename that is pure text
    
    KEfile = fopen(filenameupdated, "a");

    fprintf(KEfile, "%f %f\n", t, KE);
    
    fclose(KEfile);
    
}


/* The next two functions are get the mathematica right.[Running] python -u "/home/ryan/University/Honours/Semester 1/Statistical Physics B 721/Term 2 Simulations/data_analysis.py"
QApplication: invalid style override 'gtk2' passed, ignoring it.
	Available styles: Windows, Fusion

void tidyupmathematicafile(double pos[PARTICLES][DIM], char filename[FIILENAMELENGTH], double boxdims[DIM]){
    
    int i;
    char filenameupdated[FIILENAMELENGTH];
    
    FILE *posfile;
    
    strcpy(filenameupdated, filename);
    strcat(filenameupdated, "_mathematica.nb"); // data written to filename that is formated and designated a mathematica notebook
    
    posfile = fopen(filenameupdated, "a");

    fprintf(posfile, "ListPlot[{");
    
    for (i=0; i<(PARTICLES-1); i++) {
        fprintf(posfile, "{ %f, %f },\n", pos[i][0], pos[i][1]);
    }
    
    fprintf(posfile, "{ %f, %f }\n", pos[PARTICLES-1][0], pos[PARTICLES-1][1]);
    
    fprintf(posfile, "}, PlotStyle -> PointSize[Large], \n PlotRange -> {{0, %f}, {0, %f}}, AspectRatio -> 1]}]\n", boxdims[0], boxdims[1]);

    fclose(posfile);
    
}

void startmathematica(char filename[FIILENAMELENGTH]){
    
    char filenameupdated[FIILENAMELENGTH];
    
    FILE *posfile;
    
    strcpy(filenameupdated, filename);
    strcat(filenameupdated, "_mathematica.nb"); // data written to filename that is formated and designated a mathematica notebook
    
    posfile = fopen(filenameupdated, "w");

    fprintf(posfile, "ListAnimate[{");
    
    fclose(posfile);
    
}

/* ------------------------------------------------------------  */
/* ------------------------------------------------------------  */


int main(int argc, const char * argv[]) {
    
    double boxdims[DIM]; // array storing x and y dimensions of the simulation box
    double temperatureinitial; // temperature for distribution to sample particle speeds
    
    double currenttime=0.0, deltat, endtime; // intialise clock, declare variable for time step size, and final time
    double kineticenergy=0.0; // variable that will temporarily store the kinetic energy of all particles at a time step
    
    double positions[PARTICLES][DIM]; // 2D array to store x and y components of particle positions
    double velocities[PARTICLES][DIM]; // 2D array to store the x and y components of particle velocities
    double accelerations[PARTICLES][DIM]; // 2D array to store the accelerations.
    
    char outputfilename[FIILENAMELENGTH]; // array for file name(s) to which output is saved.
    
    int outputinterval, stepssinceoutput=0;
    
    
    
    
    /* The following steps read relevant quantities and filenames from the command line of the code */
    deltat = atof(argv[1]);
    endtime = atof(argv[2]);
    outputinterval = atof(argv[3]);
    boxdims[0] = atof(argv[4]);
    boxdims[1] = atof(argv[5]);
    temperatureinitial = atof(argv[6]);
    strcpy(outputfilename, argv[7]);
    // end: parsing the command line
    
    // TEST on command line START
    printf("code has read in delta t= %.4e end time=%.4e interval for data=%d\n box x=%.4e box y=%.4e initial t=%.4e file=%s\n",deltat, endtime, outputinterval, boxdims[0], boxdims[1], temperatureinitial, outputfilename);
    
    
    initialiseparticles(positions, velocities, temperatureinitial, boxdims);
    
    startmathematica(outputfilename);
    
    while (currenttime < endtime) {
        
        calculateaccelerations(positions, accelerations);
        
        moveparticles(positions, velocities, accelerations, deltat, &currenttime, &kineticenergy, boxdims);
        
        stepssinceoutput += 1;
        
        if (stepssinceoutput == outputinterval) {
            stepssinceoutput=0;
            
            writepositions (positions, outputfilename, boxdims);
            writeenergies (kineticenergy, currenttime, outputfilename);
            
        }
        
    }
    
    tidyupmathematicafile(positions, outputfilename, boxdims);
    
    
}
