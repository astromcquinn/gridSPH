# gridSPH

Code grids SPH particles and does simple analysis on them.

#example for how to run code is read_gadget.ipynb



To build C object file

To compile: gcc -c gridGadget.c -o gridGadget.o

To do: use openmp to make so can use many threads

to compile cython wrapper of C code that can be imported into .ipynb:
python setup.py build_ext --inplace


the power spectrum code is power_spectrum.py
