# Copyright 2017 Fred Hutchinson Cancer Research Center################################################################################
### Build cython interface for BGEN code and package as a python module

import os

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ZSTD_INC = os.environ['ZSTD_INC']
ZSTD_LIB = os.environ['ZSTD_LIB']

# '-Izstd-1.1.0/lib', '-Lzstd-1.1.0/lib'],

setup(
  name = 'bgenpy',
  ext_modules=[
    Extension('bgenpy',
              # extra_compile_args=['-O3'],
              extra_compile_args= ['-Ibgen', '-I' + ZSTD_INC, '-L' + ZSTD_LIB],
              libraries = [ "zstd", "z"],
              library_dirs=[ZSTD_LIB],
              sources=["bgenpy.pyx", "bgen_wrapper.cpp", "bgen/zlib.cpp", "bgen/bgen.cpp", "bgen/MissingValue.cpp"],  #, "bgen/bgen_to_vcf.cpp"],    
              language='c++')
    ],
  cmdclass = {'build_ext': build_ext}
)
