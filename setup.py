from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension("wrapper",
              sources=["wrapper.pyx"],
              include_dirs=[numpy.get_include()],  # Add NumPy headers to include dirs
              extra_objects=["/Users/matt/Dropbox/grad_students_projects/Hurum/gridGadget.o"],
              # Include any other necessary directories or libraries here
             )
]

setup(
    name="CythonWrapper",
    ext_modules=cythonize(extensions),
)

