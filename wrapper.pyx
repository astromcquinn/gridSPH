
cimport numpy as np
import numpy as np
from libc.stdlib cimport malloc, free


cdef extern from "/Users/matt/Dropbox/grad_students_projects/Hurum/gridGadget.h":
    struct gridStruct:
        float Rho
        float T
        float vel[3]
        # Ensure your C code can handle dynamically allocated arrays for Rho and T

    struct particle_data:
        float Pos[3];
        float Vel[3];
        float Mass;
        # int Type;
        float Rho, U, Temp, Ne, Nh, Hsml;

    void gridSPH(gridStruct *grid, particle_data *part_struct, double BoxSize, int Na, int Ngas);

def wrap_gridSPH(np.ndarray[np.float32_t, ndim=2] pos,
                 np.ndarray[np.float32_t, ndim=1] Hsml,
                 np.ndarray[np.float32_t, ndim=1] temp,
                 np.ndarray[np.float32_t, ndim=1] dens,
                 np.ndarray[np.float32_t, ndim=1] mass,
                 int Na, double BoxSize):
    cdef int Ngas = pos.shape[0]
    cdef int Ngrid = Na**3
    #cdef gridStruct c_grid

    cdef gridStruct *c_grid = <gridStruct*>malloc(sizeof(gridStruct)*Ngrid)

   

    cdef particle_data *c_particles = <particle_data*>malloc(Ngas * sizeof(particle_data))

    cdef int i, j
    # Initialize grid arrays to zero
    for i in range(Ngrid):
        c_grid[i].Rho = 0.0
        c_grid[i].T = 0.0

        for j  in range(3):
            c_grid[i].vel[j] = 0.0


    # Convert numpy arrays to particle_data array
    for i in range(Ngas):
        c_particles[i].Pos[0] = pos[i, 0]
        c_particles[i].Pos[1] = pos[i, 1]
        c_particles[i].Pos[2] = pos[i, 2]
        c_particles[i].Mass = mass[i]  # Assuming dens is used for Mass
        c_particles[i].Temp = temp[i]
        c_particles[i].Rho = dens[i]  # If Rho is also needed
        c_particles[i].Hsml = Hsml[i]  # If Rho is also needed

    # Call the C function
    gridSPH(c_grid, c_particles, BoxSize, Na, Ngas)



    # Create a numpy array to hold the grid data
    # Define a numpy dtype for the gridStruct struct
    grid_dtype = np.dtype([('Rho', np.float32), ('T', np.float32), ('vel', np.float32, (3,))])
    grid_array = np.empty(Ngrid, dtype=grid_dtype)

    # Copy the data from c_grid to grid_array (a python struct)
    for i in range(Ngrid):
        grid_array[i]['Rho'] = c_grid[i].Rho
        #print('ga ', i, grid_array[i]['Rho'])
        grid_array[i]['T'] = c_grid[i].T
        for j in range(3):
            grid_array[i]['vel'][j] = c_grid[i].vel[j]

    # Free allocated memory if not needed anymore
    free(c_particles)
    free(c_grid)

    # Return the numpy array
    return grid_array
