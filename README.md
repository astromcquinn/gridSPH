# gridSPH

#Code grids SPH particles and does simple analysis on them.


#to build C object file
#currently it does everything in double format
#To do: use openmp to make so can use many threads
gcc -c gridGadget.c -o gridGadget.o

#to compile cython wrapper
python setup.py build_ext --inplace
