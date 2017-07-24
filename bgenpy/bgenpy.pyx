# distutils: language = c++

# Copyright 2017 Fred Hutchinson Cancer Research Center

cimport cython

from libcpp.string cimport string as cpp_string
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "bgen.hpp" namespace "genfile::bgen":
    cdef cppclass Context:
        Context() except +
        int number_of_samples, number_of_variants, flags
        cpp_string magic
        cpp_string free_data



cdef extern from "bgen_wrapper.h":
    cdef cppclass BGENReader:
        BGENReader(cpp_string&) except +
        void load_samples()
        void seek_first_variant()
        void seek_to_variant(size_t)
        size_t offset()
        bool read_full_variant() except +
        bool read_minimal_variant() except +
        Context m_context
        vector[cpp_string] m_sample_ids

        cpp_string m_snp_id
        cpp_string m_rs_id
        cpp_string m_chromosome
        int m_position
        vector[cpp_string] m_alleles
        vector[vector[double]] m_probs
    cdef cppclass BGENWriter:
        BGENWriter(cpp_string&) except +
        void write_header(int, int, int)
        void write_variant(cpp_string&, cpp_string&, cpp_string&, int, vector[cpp_string], vector[vector[double]]) except +


# Ideas for improvement: 
#   * maybe it would be better to read samples into std::vector and return that.

cdef class Reader:
    cdef BGENReader* reader

    def __cinit__(self, filename):
        self.reader = new BGENReader(bytes(filename, "utf-8"))

    def __dealloc__(self):
        del self.reader

    def attributes(self):
        return { 'number_of_samples': self.reader.m_context.number_of_samples,
             'number_of_variants':self.reader.m_context.number_of_variants,
             'flags': self.reader.m_context.flags,
             'magic': self.reader.m_context.magic,
             'free_data': self.reader.m_context.free_data }
    def samples(self):
        # if flags & m_context.flags & genfile::bgen::e_SampleIdentifiers > 0: then there is block of sample IDs
        self.reader.load_samples()
        return self.reader.m_sample_ids
        pass

    def seek_first_variant(self):
        self.reader.seek_first_variant()

    def seek_to_variant_offset(self, offset):
        self.reader.seek_to_variant(offset)


    def read_full_variant(self):
        result = self.reader.read_full_variant()
        if not result:
            raise StopIteration
        return (self.reader.m_snp_id, 
                self.reader.m_rs_id, 
                self.reader.m_chromosome, 
                self.reader.m_position, 
                self.reader.m_alleles, 
                self.reader.m_probs)

    def read_minimal_variant(self):
        result = self.reader.read_minimal_variant()
        if not result:
            raise StopIteration
        return (self.reader.m_snp_id, 
                self.reader.m_rs_id, 
                self.reader.m_chromosome, 
                self.reader.m_position, 
                self.reader.m_alleles)

    def offset(self):
        return self.reader.offset()

    def __iter__(self):
        self.reader.seek_first_variant();
        return self

    def __next__(self):
        return self.read_full_variant()


cdef class Writer:
    cdef BGENWriter* writer
    def __cinit__(self, filename):
        self.writer = new BGENWriter(bytes(filename, "utf-8"))

    def __dealloc__(self):
        del self.writer

    def write_attributes(self, n_samples, n_variants, flags):
        self.writer.write_header(n_samples, n_variants, flags)

    def write_variant(self, snp_id, rs_id, chromosome, position, alleles, probs):
        self.writer.write_variant(snp_id, rs_id, chromosome, position, alleles, probs)



## Local Variables:
## mode: python
## End:
