/*  This file is part of jmers.

    jmers is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    jmers is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with jmers.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __JMERFISH_H__
#define __JMERFISH_H__

#include <iostream>
#include <string>
#include <jellyfish/file_header.hpp>
#include <jellyfish/mer_dna_bloom_counter.hpp>
#include <jellyfish/jellyfish.hpp>
#include "sequence_mers.hpp"

using jellyfish::mer_dna;
using jellyfish::mer_dna_bloom_counter;

namespace jmers {

class JellyfishDatabase {

    bool        is_empty;
    std::string description;
    std::ifstream in;
    jellyfish::file_header header;
    mer_dna_bloom_counter *bloomfilter;
    binary_query *binary;
    jellyfish::mapped_file binary_map;
    bool bloomcounter = false;
    
    template<typename Database>
    int query(const Database &db, std::string kmer)
    {
        sequence_mers       mers( this->header.canonical() );
        const sequence_mers mers_end( this->header.canonical() );
        
        mers = kmer;
        return db->check(*mers);
    }
    
  public:
    JellyfishDatabase(const char* filename) : is_empty(true) { this->open(filename); }
    JellyfishDatabase(std::string filename) : JellyfishDatabase( filename.c_str() ) {} 
    ~JellyfishDatabase() { in.close(); if (bloomcounter) { delete this->bloomfilter; } else { delete this->binary; } }
    
    std::string describe() const {
        if (is_empty)
            return "empty";
        return description;
    }
    
    int kmer = -1;
    int open(const char* filename)
    {
        this->in = std::ifstream( filename, std::ios::in|std::ios::binary );
        header = jellyfish::file_header( this->in );
        if( !this->in.good() )
        {
            // Raise exception
            std::cerr << "Failed to parse header of file '" << filename << "'" << std::endl;
            return 1;
        }
        // set kmer length - taken from the database
        this->kmer = this->header.key_len() / 2;
        mer_dna::k( this->kmer );
        
        /*
            Jellyfish databases can use a bloom counter to make a faster, less
            accurate kmer counter, or just count everything the old fashioned way.
            There's a slight difference in how the database is queried depending
            on whether it was created using the bloom filter.
        */
        if ( this->header.format() == "bloomcounter" )
        {
            this->bloomcounter = true;
            
            jellyfish::hash_pair<mer_dna> fns(this->header.matrix(1), this->header.matrix(2));
            this->bloomfilter = new mer_dna_bloom_counter(this->header.size(), this->header.nb_hashes(), this->in, fns);
            if ( !this->in.good() )
            {
                std::cerr << "Bloom filter file is truncated" << std::endl;
                return 1;
            }
        }
        else if ( this->header.format() == binary_dumper::format )
        {
            this->binary_map = jellyfish::mapped_file( filename );
            this->binary = new binary_query(this->binary_map.base() + this->header.offset(),
                            this->header.key_len(),
                            this->header.counter_len(),
                            this->header.matrix(),
                            this->header.size() - 1,
                            this->binary_map.length() - this->header.offset()
                           );
        } else
        {
            std::cerr << "Unsupported format '" << this->header.format() << "'. Must be a bloom counter or binary list." << std::endl;
            return 1;
        }
        return 0;
    }
    
    int query(std::string kmer)
    {
        if (this->bloomcounter)
            return this->query(this->bloomfilter, kmer);
        return this->query(this->binary, kmer);
    }
};

};

#endif