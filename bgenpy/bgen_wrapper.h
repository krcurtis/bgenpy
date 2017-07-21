// Copyright 2017 Fred Hutchinson Cancer Research Center
////////////////////////////////////////////////////////////////////////////////
// BGEN wrapper

#include <stdint.h>
#include <string>

#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <memory>
#include "bgen.hpp"


class BGENReader {
public:
  BGENReader(std::string &filename);
  ~BGENReader();
  void load_samples();
  void seek_first_variant();
  bool read_full_variant();
  bool read_minimal_variant();

  std::istream *m_stream;
  uint32_t m_offset;
  genfile::bgen::Context m_context;
  std::vector< std::string > m_sample_ids;


  //current SNP
  std::string m_snp_id;
  std::string m_rs_id;
  std::string m_chromosome;
  uint32_t m_position;
  std::vector<std::string> m_alleles;
  std::vector< std::vector< double > > m_probs ;

  
  //mysterious buffers needed by bgen
  std::vector< genfile::byte_t > m_buffer1, m_buffer2 ;

};


class BGENWriter {
public:
  BGENWriter(std::string &filename);
  ~BGENWriter();
  void write_header(uint32_t n_samples, uint32_t n_variants, uint32_t flags);

  void write_variant(std::string& snp_id,
                               std::string& rs_id,
                               std::string& chromosome,
                               int position,
                               std::vector<std::string>& alleles,
                     std::vector<std::vector<double> > &probs);

  std::ostream *m_stream;
  uint32_t m_offset;
  genfile::bgen::Context m_context;
  std::vector< std::string > m_sample_ids;


  //current SNP
  std::string m_snp_id;
  std::string m_rs_id;
  std::string m_chromosome;
  uint32_t m_position;
  std::vector<std::string> m_alleles;
  std::vector< std::vector< double > > m_probs ;

  
  //mysterious buffers needed by bgen
  std::vector< genfile::byte_t > m_buffer1, m_buffer2 ;

};





// Local Variables:
// mode: c++
// End:
