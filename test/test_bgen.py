# Copyright 2017 Fred Hutchinson Cancer Research Center
################################################################################
### Functional test cases for using bgenpy

import unittest

import math
import time
import contextlib

import bgenpy


import h5py
import bgenpy as bgen

import numpy as np
from netCDF4 import Dataset

def convert_matrix_to_strings(M):
    """convert numpy character matrix to list of strings"""
    strings = [];
    (nrow, ncol) = M.shape;
    for i in range(nrow):
        v = M[i];
        s = v.tostring()
        s = s.strip(b"\x00")
        strings.append(s.decode(encoding="latin-1"))
    return strings;

@contextlib.contextmanager
def timeit(task):
    """session context generator thingy to get the run-time of a block of code"""
    start_t = time.time()
    yield
    end_t = time.time()
    print("Task {task} time {s} seconds".format(task=task,s=end_t-start_t))


def convert_file(output_bgen, input_netcdf):
    BFLAGS = 5; ## magic number for compressed bgen file, using layout_1, no sample IDs
    # hdf = h5py.File(input_hdf, 'r')
    ncvars = Dataset(input_netcdf, 'r', format='NETCDF4')
    n_samples, n_snps = ncvars['Prob_AA'].shape  # samples x snps

    b_out = bgen.Writer(output_bgen)
    b_out.write_attributes(n_samples, n_snps, BFLAGS)
    with timeit('load from NetCDF'):
        matrix_snp = ncvars['SNP_Name'][:,:]
        snp_names = convert_matrix_to_strings(matrix_snp)
        snps = {i:s for i,s in enumerate(snp_names) }
        count_allele = ncvars['Count_Allele'][:]
        baseline_allele = ncvars['Baseline_Allele'][:]
        chromosome = ncvars['Chromosome'][:]
        position = ncvars['Position'][:]

        P_AA = ncvars['Prob_AA'][:,:]
        P_AB = ncvars['Prob_AB'][:,:]

        #P_BB = np.ones((n_samples, delta), dtype='float64') - P_AA - P_AB
        P_BB = 1 - P_AA - P_AB  # BGEN doesn't like negative values ...
        P_BB = np.maximum(P_BB, np.zeros(shape=P_BB.shape))

        #print("P_AA", P_AA[0,:])
        #print("P_AB", P_AB[0,:])
        #print("P_BB", P_BB[0,:])
        with timeit('write to bgen'):
            for i in range(n_snps):
                snp_id = bytes(snps[i], "utf-8")
                rs_id = b""
                chrom = bytes(str(chromosome[i]), "utf-8")
                pos = position[i]
                alleles = [ bytes(count_allele[i], "utf-8"), bytes(baseline_allele[i], "utf-8")]
                probs = [ [P_AA[s,i], P_AB[s,i], P_BB[s,i]] for s in range(n_samples)]
                b_out.write_variant(snp_id, rs_id, chrom, pos, alleles, probs)

        ncvars.close()



class TestBGEN(unittest.TestCase):

    def test_bgen1(self):
        """Test whether a bgen v1.1 file can be read in, written, and read inagain"""
        # this bgen file does not have sample names
        B = bgen.Reader("test_files/test_data1.bgen")

        attributes =  B.attributes()
        #print attributes
        #print B.samples()
        # for info in B:
        #    print info


        b_out = bgen.Writer("test_files/tmp.bgen")
        b_out.write_attributes(attributes['number_of_samples'],
                               attributes['number_of_variants'],
                               attributes['flags'])

        example = []
        with timeit('load from BGEN'):
            B.seek_first()
            try:
                while True:
                    info = B.read_full_variant()
                    example.append(info)
                    b_out.write_variant(*info)
            except StopIteration:
                del b_out

        #for info in B:
        #    example.append(info)
        #    b_out.write_variant(*info)
        #del b_out


        b_in = bgen.Reader("test_files/tmp.bgen")
        test_attributes = b_in.attributes()
        test_samples = b_in.samples()
        test_data = [ info for info in b_in]



        self.assertTrue(test_attributes == attributes)
        self.assertTrue(example == test_data)



    def test_bgen2(self):

        convert_file("test_files/tmp.bgen", "test_files/test_data2.nc")

        ncvars = Dataset("test_files/test_data2.nc", 'r', format='NETCDF4')
        #start = 0
        #end = 

        matrix_snp = ncvars['SNP_Name'][:,:]

        snp_names = convert_matrix_to_strings(matrix_snp)
        snps = { bytes(s, "utf-8"):i for i,s in enumerate(snp_names) }

        count_allele = [ bytes(x, "utf-8") for x in ncvars['Count_Allele'][:]]
        baseline_allele = [ bytes(x, "utf-8") for x in ncvars['Baseline_Allele'][:]]
        chromosome = ncvars['Chromosome'][:]
        position = ncvars['Position'][:]

        P_AA = ncvars['Prob_AA'][:,:]
        P_AB = ncvars['Prob_AB'][:,:]


        #b_in = bgen.Reader("test_files/tmp.bgen")
        #attributes = b_in.attributes()
        #data = [ info for info in b_in]
        #del b_in
        
        data = []
        b_in = bgen.Reader("test_files/tmp.bgen")
        try:
            while True:
                info = b_in.read_full_variant()
                data.append(info)
        except StopIteration:
            pass



        for snp,rsname, chrom, pos, alleles, probs in data:
            i = snps[snp]
            a1, a2 = alleles
            self.assertTrue(a1 == count_allele[i])
            self.assertTrue(a2 == baseline_allele[i])
            self.assertTrue(int(chrom) == chromosome[i])
            self.assertTrue(pos == position[i])

            p_aa, p_ab, p_bb = zip(*probs)
            p_aa = np.array(p_aa)
            p_ab = np.array(p_ab)
            p_bb = np.array(p_bb)

            self.assertTrue(np.allclose(p_aa, P_AA[:,i], atol=1e-4))
            self.assertTrue(np.allclose(p_ab, P_AB[:,i], atol=1e-4))
            #print(i)
            #print(p_bb)
            # print(1 - P_AA[:,i] - P_AB[:,i])
            P_BB = 1 - P_AA[:,i] - P_AB[:,i]
            P_BB = np.maximum(P_BB, np.zeros(shape=P_BB.shape)) # BGEN does not like negative values
            self.assertTrue(np.allclose(p_bb, P_BB, atol=1e-4))


    def test_bgen3(self):
        """Test skipping over variants"""
        # this bgen file does not have sample names
        B = bgen.Reader("test_files/test_data1.bgen")

        attributes =  B.attributes()
        expected_variants = attributes['number_of_variants']

        example = []
        with timeit('load minimal variants from BGEN'):
            B.seek_first()
            try:
                while True:
                    info = B.read_minimal_variant()
                    example.append(info)
            except StopIteration:
                pass
        self.assertTrue(len(example) == expected_variants)
