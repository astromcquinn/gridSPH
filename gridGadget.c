
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "gridGadget.h"


#define KERNEL_TABLE 10000
double  Kernel[KERNEL_TABLE+1],
               Kernel3D[KERNEL_TABLE+1],
               KernelRad[KERNEL_TABLE+1];


/* //these are in header file now
struct particle_data 
{
  float  Pos[3];
  float  Vel[3];
  float  Mass;
  //int    Type;

  float  Rho, U, Temp, Ne, Nh, Hsml;
} *P;


struct gridStruct{
  float Rho;
  float T;
};

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
