#ifndef GRIDGADGET_H
#define GRIDGADGET_H

#define float double   //converts everything to double precision


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct particle_data 
{
    float Pos[3];
    float Vel[3];
    float Mass;
    // int Type;
    float Rho, U, Temp, Ne, Nh, Hsml;
};

struct gridStruct
{
    float Rho;
    float T;
    // Assuming 'vel' is intended to be part of 'gridStruct' as your code snippet suggests modifications to 'grid[ll].vel'
    float vel[3]; 
};

// Function prototypes

void gridSPH(struct gridStruct *grid, struct particle_data *part_struct, double BoxSize, int Na, int Ngas);
void set_sph_kernel(void);
double getSeparation(double x, double y, double BoxSize);

#endif // GRIDGADGET_H
 
