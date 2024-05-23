
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gridGadget.h"

#include <math.h>


/********************************************/
//old kernel stuff
#define KERNEL_TABLE 10000
double  Kernel[KERNEL_TABLE+1],
               Kernel3D[KERNEL_TABLE+1],
               KernelRad[KERNEL_TABLE+1];


/**********************************************/


/**
 *
 * We use Price 1012.1885 kernels
 * sml in Gadget is the support big H in Price,
 *
 * u = r / H
 * q = r / h
 *
 * luckily, wk = 1 / H ** 3 W_volker(u)
 *             = 1 / h ** 3 W_price(q)
 * and     dwk = 1 / H ** 4 dw_volker/du
 *             = 1 / h ** 4 dw_price/dq
 *
 * wk_xx is Price eq 6 , 7, 8, without sigma
 *
 * the function density_kernel_wk and _dwk takes u to maintain compatibility
 * with volker's gadget.
 */
double wk_cs(DensityKernel * kernel, double q) {
    if(q < 1.0) {
        return 0.25 * pow(2 - q, 3) - pow(1 - q, 3);
    }
    if(q < 2.0) {
        return 0.25 * pow(2 - q, 3);
    }
    return 0.0;
}
double dwk_cs(DensityKernel * kernel, double q) {
    if(q < 1.0) {
        return - 0.25 * 3 * pow(2 - q, 2) + 3 * pow(1 - q, 2);
    }
    if(q < 2.0) {
        return -0.25 * 3 * pow(2 - q, 2);
    }
    return 0.0;
}
static double wk_qus(DensityKernel * kernel, double q) {
    if(q < 0.5) {
        return pow(2.5 - q, 4) - 5 * pow(1.5 - q, 4) + 10 * pow(0.5 - q, 4);
    }
    if(q < 1.5) {
        return pow(2.5 - q, 4) - 5 * pow(1.5 - q, 4);
    }
    if(q < 2.5) {
        return pow(2.5 - q, 4);
    }
    return 0.0;
}
static double dwk_qus(DensityKernel * kernel, double q) {
    if(q < 0.5) {
        return -4 * pow(2.5 - q, 3) + 20 * pow(1.5 - q, 3) - 40 * pow(0.5 - q, 3);
    }
    if(q < 1.5) {
        return -4 * pow(2.5 - q, 3) + 20 * pow(1.5 - q, 3);
    }
    if(q < 2.5) {
        return -4 * pow(2.5 - q, 3);
    }
    return 0.0;
}

static double wk_qs(DensityKernel * kernel, double q) {
    if(q < 1.0) {
        return pow(3 - q, 5) - 6 * pow(2 - q, 5) + 15 * pow(1 - q, 5);
    }
    if(q < 2.0) {
        return pow(3 - q, 5)- 6 * pow(2 - q, 5);
    }
    if(q < 3.0) {
        return pow(3 - q, 5);
    }
    return 0.0;
}
static double dwk_qs(DensityKernel * kernel, double q) {
    if(q < 1.0) {
        return -5 * pow(3 - q, 4) + 30 * pow(2 - q, 4)
             - 75 * pow (1 - q, 4);
    }
    if(q < 2.0) {
        return -5 * pow(3 - q, 4) + 30 * pow(2 - q, 4);
    }
    if(q < 3.0) {
        return -5 * pow(3 - q, 4);
    }
    return 0.0;
}
    

static struct {
    const char * name;
    double (*wk)(DensityKernel * kernel, double q);
    double (*dwk)(DensityKernel * kernel, double q);
    double support; /* H / h, see Price 2011: arxiv 1012.1885*/
    double sigma[3];
} KERNELS[] = {
    { "CubicSpline", wk_cs, dwk_cs, 2.,
        {2 / 3., 10 / (7 * M_PI), 1 / M_PI} },
    { "QuinticSpline", wk_qs, dwk_qs, 3.,
        {1 / 120., 7 / (478 * M_PI), 1 / (120 * M_PI)} },
    { "QuarticSpline", wk_qus, dwk_qus, 2.5,
        {1 / 24., 96 / (1199 * M_PI), 1 / (20 * M_PI)} },
};


double
density_kernel_dwk(DensityKernel * kernel, double u)
{
    double support = KERNELS[kernel->type].support;
    return kernel->dWknorm *
        KERNELS[kernel->type].dwk(kernel, u * support);
}

double
density_kernel_wk(DensityKernel * kernel, double u)
{
    double support = KERNELS[kernel->type].support;
    return kernel->Wknorm *
        KERNELS[kernel->type].wk(kernel, u * support);
}



/*

// Example modification of gridSPH to use the new kernel function
void gridSPH(struct gridStruct *grid, struct particle_data *part_struct, double BoxSize, int Na, int Ngas) {
    int ip, jp, kp;
    float dx[3], r;
    double u, wk, weight, vel;
    float h, Mass, Temp;
    int cell[3], num_sur;
    double cell_size = BoxSize / Na;
    set_sph_kernel(); //sets the old kernel table

    // Get the kernel weight using the new function
    DensityKernel kernel;
    kernel.type =  1; // I believe type 2 is the correct kernel if we are using the sophisticated SPH (the quintic kerrnel)
    kernel.Wknorm = 1.0 / (h * h * h); // Normalization might need adjustment
    
    printf("Name: %s\n", KERNELS[1].name);
    printf("Support: %f\n", KERNELS[1].support);
    

    for (int i = 0; i < Ngas; i++) {
        struct particle_data *p = &part_struct[i];
        ip = (int)(Na * (p->Pos[0] + 0.5) / BoxSize);
        jp = (int)(Na * (p->Pos[1] + 0.5) / BoxSize);
        kp = (int)(Na * (p->Pos[2] + 0.5) / BoxSize);

        Mass = p->Mass * 1e10; // Convert from gadget unit
        Temp = p->Temp;
        h = p->Hsml;
        num_sur = (int)(h / cell_size + 0.5);

        for (cell[0] = ip - num_sur; cell[0] <= ip + num_sur; cell[0]++)
            for (cell[1] = jp - num_sur; cell[1] <= jp + num_sur; cell[1]++)
                for (cell[2] = kp - num_sur; cell[2] <= kp + num_sur; cell[2]++) {
                    int ll = ((cell[0] + Na) % Na) * Na * Na + ((cell[1] + Na) % Na) * Na + (cell[2] + Na) % Na;

                    for (int j = 0; j < 3; j++) {
                        dx[j] = getSeparation(p->Pos[j], cell[j] * cell_size, BoxSize);
                    }

                    r = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);

                    if (r < h) {
                        u = r / h;


                        //wk = density_kernel_wk(&kernel, u);

                        int ii = (int) (u * KERNEL_TABLE);
                        wk = 1.0 / (h*h*h) * (Kernel3D[ii] +(Kernel3D[ii + 1] - Kernel3D[ii]) 
                            * (u - KernelRad[ii])* KERNEL_TABLE);
                        

                        grid[ll].Rho += Mass * wk;
                        weight = Mass * wk;

                        //probably need to update these with correct kernel function?
                        grid[ll].T += Temp * weight;

                        for (int kvel = 0; kvel < 3; kvel++) {
                            grid[ll].vel[kvel] += p->Vel[kvel] * weight;
                        }
                    }
                }
    }
}

*/

//taken from Lya forest code
//Na is nnumber of cells across
void gridSPH(struct gridStruct *grid, struct particle_data *part_struct, double BoxSize, int Na, int Ngas)
{

  int j, ii, ll, ip, jp, kp;
  float dx[3], r;
  float u, wk, weight, vel;
  float h, Mass, Temp;
  int cell[3], num_sur;
  double cell_size = BoxSize/Na;
  set_sph_kernel();


  for(int i=0; i<Ngas; i++)
  {
    struct particle_data *p = &part_struct[i];
    ip = (int) (Na*(p->Pos[0]+.5)/BoxSize);
    jp=  (int) (Na*(p->Pos[1]+.5)/BoxSize);
    kp = (int) (Na*(p->Pos[2]+.5)/BoxSize);

    Mass = p->Mass*1e10; //convert from gadget unit (although h's are still there)
    Temp = p->Temp;
    h = p->Hsml;

    num_sur = (int) (h/cell_size + .5);


    //printf("%d num_sur = %d", i, num_sur);  

    for(cell[0] = ip-num_sur; cell[0] <= ip+num_sur; cell[0]++)
      for(cell[1] = jp-num_sur; cell[1] <= jp+num_sur; cell[1]++)
        for(cell[2] = kp-num_sur; cell[2] <= kp+num_sur; cell[2]++)
    { 
      ll = ((cell[0]+Na)%Na)*Na*Na+((cell[1]+Na)%Na)*Na
        +(cell[2]+Na)%Na;

      // fprintf(stdout, "%d %d %d\n", cell[0], cell[1], cell[2]); 
      for(j=0; j<3; j++) 
        dx[j] = getSeparation(p->Pos[j], (cell[j])*cell_size, BoxSize);

      r = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
        
      if(r < h)
        {

          u = r / h;
          ii = (int) (u * KERNEL_TABLE);

          wk = 1.0 / (h*h*h) * (Kernel3D[ii] +(Kernel3D[ii + 1] - Kernel3D[ii]) 
              * (u - KernelRad[ii])* KERNEL_TABLE);

          //   printf("%d %d %d %f %f\n", cell[0], cell[1], cell[2], r/h, wk);


        

          grid[ll].Rho += Mass * wk;//will have to give right units eventually
          weight = Mass * wk;

          grid[ll].T += Temp * weight;

          for(int kvel=0.; kvel<3; kvel++)
            grid[ll].vel[kvel] += p->Vel[kvel]*weight;
        }	
    }
  }
   
  //   exit(-5);
}

void gridCIC_DM(struct gridStruct *grid, struct particle_data_DM *part_struct_DM, double BoxSize, int Nmesh, int NDM)
{
  int i, j, k, ii, jj, kk;
  float f1,f2,f3,f4,f5,f6,f7,f8, u,v,w, val;

  for(int iDM=0; iDM<NDM; iDM++)
    {
        struct particle_data_DM *p = &part_struct_DM[iDM];
      u = p->Pos[0] / BoxSize * Nmesh;
      v = p->Pos[1] / BoxSize * Nmesh;
      w = p->Pos[2] / BoxSize * Nmesh;
            
      i = (int) u;
      j = (int) v;
      k = (int) w;




      //  fprintf(stderr, "%le %le %le %d %d %d\n", u, v, w, i, j, k);
            
      if(i >= Nmesh)
        i = Nmesh - 1;
      if(j >= Nmesh)
        j = Nmesh - 1;
      if(k >= Nmesh)
        k = Nmesh - 1;
            
      u -= i;
      v -= j;
      w -= k;
            
      ii = i + 1;
      jj = j + 1;
      kk = k + 1;
            
      if(ii >= Nmesh)
        ii -= Nmesh;
      if(jj >= Nmesh)
        jj -= Nmesh;
      if(kk >= Nmesh)
        kk -= Nmesh;
            
      f1 = (1 - u) * (1 - v) * (1 - w);
      f2 = (1 - u) * (1 - v) * (w);
      f3 = (1 - u) * (v) * (1 - w);
      f4 = (1 - u) * (v) * (w);
      f5 = (u) * (1 - v) * (1 - w);
      f6 = (u) * (1 - v) * (w); 
      f7 = (u) * (v) * (1 - w);
      f8 = (u) * (v) * (w);
              
      val  = p->Mass; //assumes same number of DM particles
      grid[(i * Nmesh + j) * Nmesh + k].Rho += val*f1;
      grid[(i * Nmesh + j) * Nmesh + kk].Rho += val*f2;
      grid[(i * Nmesh + jj) * Nmesh + k].Rho += val*f3;
      grid[(i * Nmesh + jj) * Nmesh + kk].Rho += val*f4;
      grid[(ii * Nmesh + j) * Nmesh + k].Rho += val*f5;
      grid[(ii * Nmesh + j) * Nmesh + kk].Rho  += val*f6;
      grid[(ii * Nmesh + jj) * Nmesh + k].Rho  += val*f7;
      grid[(ii * Nmesh + jj) * Nmesh + kk].Rho += val*f8;


      //  fprintf(stdout, "sum f = %le\n", f1+f2+f3+f4+f5+f6+f7+f8);
    }
}



/************************************
old sph kenrels
***************************************/


/******************************************************
returns the separation between x and y assuming a periodic box
**************************************************************/
double getSeparation(double x, double y, double BoxSize)
{
  float lbox = BoxSize;
  float sep = x - y, sep2 = x +lbox-y, sep3=x-y-lbox;

  sep2 = (sep2*sep2 < sep3*sep3) ? sep2 : sep3;
  
  return (sep*sep < sep2*sep2) ? sep: sep2;
}



/******************************************************************************************
Gadget SPH kernel
*****************************************************************************************/
void set_sph_kernel(void)	/* `Kernel' with 2D normalization, `Kernel3D' with 3D normalization  */
{
  int i;

  for(i = 0; i <= KERNEL_TABLE; i++)
    KernelRad[i] = ((double) i) / KERNEL_TABLE;

  for(i = 0; i <= KERNEL_TABLE; i++)
    {
      if(KernelRad[i] <= 0.5)
	{
	  // Kernel[i] = 40 / (7.0 * PI) * (1 - 6 * KernelRad[i] * KernelRad[i] * (1 - KernelRad[i]));
	  Kernel3D[i] = 8 / M_PI * (1 - 6 * KernelRad[i] * KernelRad[i] * (1 - KernelRad[i]));
	}
      else
	{
	  //Kernel[i] = 40 / (7.0 * PI) * 2 * (1 - KernelRad[i]) * (1 - KernelRad[i]) * (1 - KernelRad[i]);
	  Kernel3D[i] = 8 / M_PI * 2 * (1 - KernelRad[i]) * (1 - KernelRad[i]) * (1 - KernelRad[i]);
	}
    }
}