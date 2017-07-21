# Copyright 2017 Fred Hutchinson Cancer Research Center
################################################################################
### Makefile for building BGEN to VCF example program


ZSTD_DIR =  # ... zstd-1.1.0/lib
ZSTD_INC = -I$(ZSTD_DIR)
ZSTD_LIB = -L$(ZSTD_DIR) -lzstd


INCS = $(ZSTD_INC)
LIBS = $(ZSTD_LIB) -lz

FLAGS = -g -std=c++11
bgen_to_vcf: bgen_to_vcf.cpp
	g++ $(FLAGS) $(INCS) -o bgen_to_vcf *.cpp  $(LIBS)
