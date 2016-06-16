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


/*
    This file is heavily based on the Jellyfish "query_per_sequence.cc" example.
    If it breaks the GPL to use quite a bit of it as a template, please let me
    know!

    / Martin Norling
*/

//#include "Seq.h"
#include "Input.h"
#include <zlib.h>
#include <iostream>

#include <vector>

#include <jellyfish/file_header.hpp>
#include <jellyfish/mer_dna_bloom_counter.hpp>
#include <jellyfish/jellyfish.hpp>
#include "sequence_mers.hpp"

using jellyfish::mer_dna;
using jellyfish::mer_dna_bloom_counter;

int split_sequence_on_kmers(std::vector<int> &kmer_content, std::vector<int> &split_positions)
{
    std::cout << "kmers: ";
    bool last = false;
    int i = 0;
    for ( auto p : kmer_content )
    {
        if ( last != p )
            split_positions.push_back( i );
        last = p;
        i++;
    }
    if (last)
        split_positions.push_back( i );
    
    return 0;
}

template<typename Database>
int kmer_content(const Database& db, jmers::Seq s, bool canonical)
{
    sequence_mers                           mers(canonical);
    const sequence_mers                     mers_end(canonical);
    
    std::cout << ">" << s.name << "\n";
    mers = s.sequence;
    
    std::vector<int> kmer_content;
    std::vector<int> split_positions;
    for ( ; mers != mers_end; ++mers )
        kmer_content.push_back(db.check(*mers));
    
    split_sequence_on_kmers( kmer_content, split_positions );
    
    
    std::cout << "kmers: ";
    for ( auto p : kmer_content )
    {
        std::cout << p << " ";
    }
    std::cout << "\n";
    
    std::cout << "split positions: ";
    for ( auto p : split_positions )
    {
        std::cout << p << " ";
    }
    std::cout << "\n";
    
    return 0;
}

int main(int argc, char *argv[])
{
    if ( argc < 3 )
    {
        std::cerr << "USAGE: " << argv[0] << " jellyfish_db.jf seq_file [seq_file [...]]" << std::endl;
        return 0;
    }
    
    // fastq-parsing variables
    jmers::Seq s;
    
    // Parse jellyfish database
    std::ifstream in(argv[1], std::ios::in|std::ios::binary);
    jellyfish::file_header header(in);
    if(!in.good())
    {
        std::cerr << "Failed to parse header of file '" << argv[1] << "'" << std::endl;
        return 1;
    }

    // set kmer length - taken from the database
    mer_dna::k(header.key_len() / 2);

    /*
        Jellyfish databases can use a bloom counter to make a faster, less
        accurate kmer counter, or just count everything the old fashioned way.
        There's a slight difference in how the database is queried depending
        on whether it was created using the bloom filter.
    */
    if ( header.format() == "bloomcounter" )
    {
        jellyfish::hash_pair<mer_dna> fns(header.matrix(1), header.matrix(2));
        mer_dna_bloom_counter filter(header.size(), header.nb_hashes(), in, fns);
        if ( !in.good() )
        {
            std::cerr << "Bloom filter file is truncated" << std::endl;
            return 1;
        }
        in.close();

        // Parse through everything with a bloom counter database
        for (int i = 2; i < argc; i++)
        {
            try
            {
                jmers::Input seq_file = jmers::Input(argv[i]);
                while ( seq_file.read(s) )
                {
                    s.dump(std::cout);
                    kmer_content(filter, s, header.canonical());
                }
            }
            catch (std::exception &e)
            {
                std::cerr << e.what() << std::endl;
            }
        }
    }
    else if ( header.format() == binary_dumper::format )
    {
        jellyfish::mapped_file binary_map(argv[1]);
        binary_query bq(binary_map.base() + header.offset(),
                        header.key_len(),
                        header.counter_len(),
                        header.matrix(),
                        header.size() - 1,
                        binary_map.length() - header.offset()
                       );

        // Parse through everything with a 'regular' database
        for (int i = 2; i < argc; i++)
        {
            try
            {
                jmers::Input seq_file = jmers::Input(argv[i]);
                while ( seq_file.read(s) )
                {
                    s.dump(std::cout);
                    kmer_content(bq, s, header.canonical());
                }
            }
            catch (std::exception &e)
            {
                std::cerr << e.what() << std::endl;
            }
        }
    } else
    {
        std::cerr << "Unsupported format '" << header.format() << "'. Must be a bloom counter or binary list." << std::endl;
        return 1;
    }
    
    return 0;
}
