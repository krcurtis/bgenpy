// Copyright 2017 Fred Hutchinson Cancer Research Center

#include <stdio.h>

// #define DEBUG_BGEN_FORMAT 1
#define HAVE_ZLIB 1


#include "bgen_wrapper.h"
#include "bgen/prob_setter.hpp"

class SampleWrap {
public:
  SampleWrap(BGENReader *B) {
    m_B = B;
  }

  void operator()(std::string &id) {
    m_B->m_sample_ids.push_back( id );
    //return true;
  }

private:
  BGENReader *m_B;
};



class SNPWrap {
public:
  SNPWrap( std::vector< std::string >* alleles) { m_alleles = alleles; }

  void operator()( std::size_t n ) { m_alleles->resize( n ) ; }
  void operator()( std::size_t i, std::string const& allele ) { m_alleles->at(i) = allele ; }

private:
  std::vector< std::string >* m_alleles;
};


BGENReader::BGENReader(std::string &filename)
{
  //printf("Here %s\n", filename.c_str());
  m_stream = new std::ifstream(filename.c_str(), std::ios_base::binary);
  genfile::bgen::read_offset( *m_stream, &m_offset ) ;
  //printf("offset %i\n", m_offset);
  genfile::bgen::read_header_block( *m_stream, &m_context ) ;
}

BGENReader::~BGENReader()
{
  //printf("BGENReader destructor called\n");
  delete m_stream;
}


void BGENReader::load_samples()
{
  m_sample_ids.clear();
  if (m_context.flags & genfile::bgen::e_SampleIdentifiers) {
    SampleWrap wrap(this);
    genfile::bgen::read_sample_identifier_block(*m_stream, m_context, wrap);
  }
}

void BGENReader::seek_first_variant()
{
  //printf("Here %s\n", __FUNCTION__);
  m_stream->seekg(m_offset + 4);  // why + 4?
}

bool BGENReader::read_full_variant()
{
  //printf("Here %s\n", __FUNCTION__);
  m_alleles.clear();
  SNPWrap wrap(&m_alleles);

  bool result = genfile::bgen::read_snp_identifying_data(*m_stream,
                                                         m_context,
                                                         &m_snp_id,
                                                         &m_rs_id,
                                                         &m_chromosome,
                                                         &m_position,
                                                         wrap, wrap);
  //have to read probs to advance file pointer to next variant, or call some special function ...
  if (!result)
    return result;

  m_probs.clear();
  ProbSetter setter(&m_probs);
  genfile::bgen::read_and_parse_genotype_data_block< ProbSetter >(
			*m_stream,
			m_context,
			setter,
			&m_buffer1,
			&m_buffer2
		) ;

  return result;
}

bool BGENReader::read_minimal_variant()
{
 //printf("Here %s\n", __FUNCTION__);
  m_probs.clear();
  m_alleles.clear();
  SNPWrap wrap(&m_alleles);

  bool result = genfile::bgen::read_snp_identifying_data(*m_stream,
                                                         m_context,
                                                         &m_snp_id,
                                                         &m_rs_id,
                                                         &m_chromosome,
                                                         &m_position,
                                                         wrap, wrap);
  if (!result)
    return result;

  //skip over probs data
  genfile::bgen::ignore_genotype_data_block(*m_stream, m_context);

  return result;

}


BGENWriter::BGENWriter(std::string &filename)
{
  //printf("Here %s\n", filename.c_str());
  m_stream = new std::ofstream(filename.c_str(), std::ios_base::binary);

}

BGENWriter::~BGENWriter()
{
  //printf("BGENWriter destructor called\n");
  delete m_stream;
}

void BGENWriter::write_header(uint32_t n_samples, uint32_t n_variants, uint32_t flags)
{
  m_context.number_of_samples = n_samples;
  m_context.number_of_variants = n_variants;
  m_context.flags = flags;
  uint32_t header_size = m_context.header_size();
  m_offset = header_size;
  genfile::bgen::write_offset(*m_stream, m_offset);
  genfile::bgen::write_header_block(*m_stream, m_context);
}


//TODO implement writing samples
// write_sample_identifier_block()


class AlleleProvider {
public:
  AlleleProvider( std::vector< std::string >* alleles) { m_alleles = alleles; }
  std::string operator()( std::size_t i) { return m_alleles->at(i); }

private:
  std::vector< std::string >* m_alleles;
};


void BGENWriter::write_variant(std::string& snp_id,
                               std::string& rs_id,
                               std::string& chromosome,
                               int position,
                               std::vector<std::string>& alleles,
                               std::vector<std::vector<double> > &probs)
{


  AlleleProvider ap(&alleles);

  genfile::byte_t* end = genfile::bgen::write_snp_identifying_data(&m_buffer1, 
                                         m_context, 
                                         snp_id, 
                                         rs_id, 
                                         chromosome, 
                                         position,
                                         alleles.size(),
                                                                 ap);

  m_stream->write( reinterpret_cast< char* >( &(m_buffer1[0]) ), end - &(m_buffer1[0]) ) ;

  //magic happens
  genfile::bgen::GenotypeDataBlockWriter bw(&m_buffer1, &m_buffer2, m_context, 16);
  //printf("before initialize %zu %zu\n", probs.size(), alleles.size());
  bw.initialise(probs.size(), alleles.size());

  //printf("Probs size %zu\n", probs.size());
  for (size_t i = 0; i < probs.size(); i++) {
    bw.set_sample(i);
    bw.set_number_of_entries(2, 3, genfile::ePerUnorderedGenotype, genfile::eProbability);
    for (size_t j = 0; j < probs[i].size(); j++) {
      bw.set_value(j, probs[i][j]);
    }
  }
  bw.finalise();
  //magic ends

  char const* start = reinterpret_cast< char const* >( bw.repr().first );
  size_t len = bw.repr().second  - bw.repr().first;
  //printf("%p %zu\n", start, len);
  m_stream->write( start, len);
  // reinterpret_cast< char const* >( bw.repr().first ), bw.repr().second  - bw.repr().first 

}




//byte_t* write_snp_identifying_data()
// how to write probability data ?????

// Local Variables:
// mode: c++
// End:
