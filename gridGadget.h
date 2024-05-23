#ifndef GRIDGADGET_H
#define GRIDGADGET_H


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


struct particle_data_DM 
{
    float Pos[3];
    float Vel[3];
    float Mass;
};


struct gridStruct
{
    float Rho;
    float T;
    // Assuming 'vel' is intended to be part of 'gridStruct' as your code snippet suggests modifications to 'grid[ll].vel'
    float vel[3]; 
};

enum DensityKernelType {
    DENSITY_KERNEL_CUBIC_SPLINE = 1,
    DENSITY_KERNEL_QUINTIC_SPLINE = 2,
    DENSITY_KERNEL_QUARTIC_SPLINE = 4,
};

typedef struct {
    double H;
    double HH;
    double Hinv; /* convert from r to u*/
    int type;
    double support;
    const char * name;
    /* private: */
    double Wknorm;
    double dWknorm;
} DensityKernel;



void gridSPH(struct gridStruct *grid, struct particle_data *part_struct, double BoxSize, int Na, int Ngas);
void gridCIC_DM(struct gridStruct *grid, struct particle_data_DM *part_struct_DM, double BoxSize, int Nmesh, int NDM);
void set_sph_kernel(void);
double getSeparation(double x, double y, double BoxSize);

double density_kernel_dwk(DensityKernel * kernel, double u);
double density_kernel_wk(DensityKernel * kernel, double u);




#endif // GRIDGADGET_H
 
